#include "mpi.h"
#include <stdio.h>
#include <vector>
#include <list>
#include <iostream>
#include "../utility/debug.cpp"
#include "../utility/asserts.cpp"
#include "../utility/errors.cpp"
#include "../utility/MPIdata.cpp"
using namespace std;

//int MPIrank;
//int MPInumprocs;

#define debugout cout << "(" << MPIdata::get_rank() << ")" \
  << "DEBUG " << __func__ << ", " << __FILE__ << ":" << __LINE__ << ": "

// hacked mechanisms to print output one process at a time
//
// for some reason I still get threads overwriting one another,
// although not as badly as without this.
//
void output_barrier()
{
  fflush(stdout);
  cout << flush;
  MPI_Barrier(MPI_COMM_WORLD);
}
int barrier_ret_nprocs()
{
  output_barrier();
  return MPIdata::get_nprocs();
}
int barrier_ret_1()
{
  output_barrier();
  return 1;
}
#define criticalout \
  for(int rank=0;rank<barrier_ret_nprocs();rank++) \
      if(rank==MPIdata::get_rank()) 
#define masterout \
  for(int i=0;i<barrier_ret_1();i++) \
    if(!MPIdata::get_rank())
// Creating MPI barriers to atomize output is a poor solution;
#define dout criticalout debugout
// a better solution is to write output to separate files:
//ofstream outfile;
//#define dout outfile;

// particle ("element")
struct Particle
{
  double u[8];
};

std::ostream& operator<<(std::ostream& os, const Particle& pcl)
{
  //os << showpos;
  os << "[";
  for(int i=0; i<7; i++)
    os << pcl.u[i] << " ";
  os << pcl.u[7] << "]"; // << "\n";
  return os;
}

namespace Comm
{
  enum Enum
  {
    FINISHED = 0,
    UNFINISHED = 1
  };
};

// The combination of group and neighbor
// should be unique for each connection
//
// In a Cartesian topology, the challenges to this uniqueness
// occur for periodic boundary conditions when there is a
// dimension only 1 or 2 processes thick.
//
// In the 1-process-thick case, we can forgo MPI communication
// in that direction altogether.  But in the 2-process case,
// we need to make a distinction between upward and downward
// channels of communication.  Distinguishing upward and
// downward directions is sufficient for any topology based
// on convex MPI subdomains.
//
class Connection
{
 private:
  MPI_Comm _comm;
  int _rank; // rank within group of neighbor we're connecting to
 public:
  Connection(const MPI_Comm& comm_,int rank_):
    _comm(comm_),
    _rank(rank_)
  {}
  MPI_Comm comm()const{return _comm;}
  int rank()const{return _rank;}
};

// block of elements
template <class type>
struct Block
{
 private: // initialized in constructor
  // change the block to use an aligned allocator
  vector<type> block;
  int capacity;
  int id;
  // used for MPI communication
  MPI_Request request;
  int flag;
 private: // initialized at compile time
  // assumes using MPI_DOUBLE
  static const int NUMBERS_PER_ELEMENT = sizeof(type)/sizeof(double);
 public:
  Block(int capacity_, int id_):
    capacity(capacity_),
    id(id_),
    request(MPI_REQUEST_NULL),
    flag(0)
  {
    block.reserve(capacity);
  }
 public: // accessors
  MPI_Request& fetch_request(){return request;}
  vector<type>& fetch_block(){return block;}
  const vector<type>& get_block()const{return block;}
  int get_block_id(){ return id; }
 public: // sending
  bool isfull()
  {
    bool retval = (block.size()>=capacity);
    return retval;
  }
  void push_back(const type& in)
  {
    block.push_back(in);
  }
  void send(const Connection& dest)
  {
    assert_le(block.size(), capacity);
    MPI_Isend(&block[0], NUMBERS_PER_ELEMENT*block.size(), MPI_DOUBLE,
      dest.rank(), 0, dest.comm(), &request);
  }
 public: // receiving
  // post a receive
  void recv(const Connection& source)
  {
    // make sure that space exists to receive
    block.resize(capacity);
    int tag=0;
    MPI_Irecv(&block[0], NUMBERS_PER_ELEMENT*block.size(), MPI_DOUBLE,
      source.rank(), tag, source.comm(), &request);
  }
  // processing received data
  //
  // returns true if shrink was necessary.
  bool shrink_received_block(MPI_Status& status)
  {
    int count;
    MPI_Get_count(&status, MPI_DOUBLE, &count);
    const int num_elements_received = count / NUMBERS_PER_ELEMENT;
    if(num_elements_received == capacity)
      return false;
    assert_lt(num_elements_received,capacity);
    // shrink to eliminate the trailing garbage
    block.resize(num_elements_received);
    return true;
  }
  void waitfor_recv()
  {
    MPI_Status status;
    MPI_Wait(&request, &status);
    shrink_received_block(status);
  }
  // returns true if message has been received.
  bool test_recv(MPI_Status& status)
  {
    MPI_Test(&request, &flag, &status);
    if(!flag)
      return false;
    // MPI_Test man page says this should now be true
    assert(request==MPI_REQUEST_NULL);
    return true;
  }
  // if communication is complete then shrinks block
  // if necessary and returns true
  //bool test(MPI_Status& status)
  //{
  //  if(flag) // communication complete
  //  {
  //    // MPI_Test man page says this should now be true
  //    assert(request==MPI_REQUEST_NULL);
  //    shrink_received_block(status);
  //    return true;
  //  }
  //  return false;
  //}
};

std::ostream& operator<<(std::ostream& os, const Block<Particle>& block_)
{
  const vector<Particle>& block = block_.get_block();
  for(int k=0; k<block.size();k++)
  {
    os << "\n  block[" << k << "] = " << block[k];
  }
  return os;
}

// We use the same class for sending and receiving, because
// in the case where a dimension is only one process thick,
// the sender and receiver are the same process (in which case
// the send and receive methods do not need to be called at
// all).  A more correct way to handle this would be to inherit
// from BlockSender_interface and BlockRecver_interface abstract
// base classes.
//
template <class type>
class BlockCommunicator
{
 private: // data
  list<void*>::iterator curr_block;
  list<void*> blockList;
  //list<Block<type>*>::iterator curr_block;
  //list<Block<type>*> blockList;
  Connection connection;
  int blocksize;
  int nextid;
 // methods
 public: // construction
  BlockCommunicator(int blocksize_, int numblocks, Connection connection_):
    blocksize(blocksize_),
    connection(connection_)
  {
    assert(numblocks>0);
    for(nextid=0;nextid<numblocks;nextid++)
    {
      Block<type>* newBlock = new Block<type>(blocksize, nextid);
      blockList.push_back(newBlock);
    }
    curr_block = blockList.begin();
  }
 //protected: // if destructor is public then it should be virtual
  ~BlockCommunicator()
  {
    list<void*>::iterator i;
    for(i=blockList.begin(); i != blockList.end(); ++i)
      delete &fetch_block(i);
  }
 // access
 private:
  Block<type>& fetch_block(list<void*>::iterator block)
  {
    return *((Block<type>*)*block);
  }
 public: // information
  // access current block
  Block<type>& fetch_curr_block()
  {
    assert(curr_block != blockList.end());
    return fetch_block(curr_block);
    //return *((Block<type>*)*curr_block);
  }
  bool at_end() const
  {
    return curr_block == blockList.end();
  }
 //private: // access
 // int get_current_block_id()
 // { return at_end() ? -1 : fetch_curr_block().get_block_id(); }

 // operations
 //
 public:
  void rewind()
  {
    curr_block = blockList.begin();
  }
  void advance_block()
  {
    assert(!at_end());
    curr_block++;
  }
 public: // receiving operations
  // this assumes that the list of blocks is not empty
  void allocate_more_blocks(int numblocks=1)
  {
    // allocate more blocks and insert them
    //
    // to make sure that curr_block points to a valid block at
    // the conclusion, we decrement it prior to the push_backs...
    curr_block--;
    for(int i=0;i<numblocks;i++)
    {
      Block<type>* newBlock = new Block<type>(blocksize, nextid++);
      blockList.push_back(newBlock);
    }
    // ...and then increment it again:
    curr_block++;
  }
  void allocate_more_recv_blocks(int numblocks=1)
  {
    // allocate blocks
    allocate_more_blocks(numblocks);

    // post receives on all blocks allocated
    list<void*>::iterator i;
    for(i=curr_block;i!=blockList.end();i++)
    {
      fetch_block(i).recv(connection);
    }
  }
 public: // receiving operations
  bool test_recv_curr_block(MPI_Status& status)
  {
    return fetch_curr_block().test_recv(status);
  }
  // handle received block and advance curr_block
  //
  // assumes that curr_block was just received with argument status
  void handle_received_block(MPI_Status& status)
  {
    bool islast = fetch_curr_block().shrink_received_block(status);
    if(islast)
    {
      // this was the last block, so cancel and free
      // the remaining receive requests
      for(++curr_block; curr_block!=blockList.end();curr_block++)
      {
        MPI_Request& pending_request = fetch_curr_block().fetch_request();
        MPI_Cancel(&pending_request);
        MPI_Request_free(&pending_request);
      }
      // curr_block = blockList.end();
    }
    else
    {
      // this was not the last block, so advance curr_block,
      // allocating more blocks to receive if necessary
      //
      advance_to_next_recv_block();
    }
  }
  void advance_to_next_recv_block()
  {
    assert(!at_end());
    curr_block++;
    if(at_end())
    {
      allocate_more_recv_blocks();
    }
  }
  MPI_Request get_curr_request()
  {
    if(at_end())
      return MPI_REQUEST_NULL;
    return fetch_curr_block().fetch_request();
  }

 // sending routines
 //
 private: // methods
  // assumes that there is a block to send
  void send_block()
  {
    // send the block
    fetch_curr_block().send(connection);

    // proceed to the next block
    curr_block++;

    // alternatively, we could recycle blocks
    // for which MPI_Test(receive) returns true.
    if(curr_block == blockList.end())
    {
      allocate_more_blocks();
    }
  }
 public: // interface
  // returns true iff block was sent
  // (whereupon the caller might want to check for incoming messages).
  bool send(const type& in)
  {
    // append the particle to the block.
    fetch_curr_block().push_back(in);

    // if the block is full, send it
    if(fetch_curr_block().isfull())
    {
      send_block();
      return true;
    }
    return false;
  }
  // send the remaining particles, sending an empty message if
  // necessary to signal that no more particles will be sent.
  //
  // this assumes that curr_block exists and has all unsent particles
  void send_complete()
  {
    assert(curr_block != blockList.end());
    assert(!fetch_curr_block().isfull());
    fetch_curr_block().send(connection);
    curr_block = blockList.end();
  }
 
 public: // receiving methods
  // post receive on all blocks starting from starting_block
  bool recv_blocks(list<void*>::iterator starting_block)
  {
    for(; starting_block != blockList.end(); ++starting_block)
      fetch_block(starting_block).recv(connection);
  }
  bool recv_blocks()
  {
    // post receives on all blocks
    recv_blocks(blockList.begin());
    // reset curr_block to initial block
    curr_block = blockList.begin();
  }
  void waitfor_recv_curr_block(MPI_Status& status)
  {
    fetch_curr_block().waitfor_recv(status);
  }
};

//std::ostream& operator<<(std::ostream& os, const BlockCommunicator<Particle>& bc)
//{
//  for(bc.rewind();!bc.at_end();bc.advance_to_next_block())
//  {
//    Block<Particle>& curr_block = bc.fetch_curr_block();
//    os << curr_block;
//  }
//  return os;
//}


//template <class type>
//struct BlockSender : public BlockCommunicator<type>
//{
//};

//template <class type>
//struct BlockRecver : public BlockCommunicator<type>
//{
//};

int main(int argc, char **argv)
{
  MPIdata::init(&argc, &argv);

  //char outfilename[30];
  //snprintf(outfilename,"testcomm.%d",MPIdata::get_rank());
  //outfile.open(outfilename);

  //MPI_Init(&argc, &argv);
  //MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  //MPI_Comm_size(MPI_COMM_WORLD, &MPInumprocs);

  // define separate communicators for opposite
  // directions of information flow
  MPI_Comm up_comm;
  MPI_Comm dn_comm;
  const int reorder = 1;
  const int ndims = 1;
  /*const*/ int dims[ndims] = {2};
  /*const*/ int periods[ndims] = {1};
  const int XDIR = 0;
  const int UPWARD = 1;
  const int DNWARD = -1;
  // create a communicator for upward communication
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &up_comm);
  // create a communicator for downward communication
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 0, &dn_comm);

  int up_cart_rank;
  int dn_cart_rank;
  MPI_Comm_rank(up_comm, &up_cart_rank);
  MPI_Comm_rank(dn_comm, &dn_cart_rank);
  criticalout dprint(up_cart_rank);
  criticalout dprint(dn_cart_rank);

  int up_coords[ndims];
  int dn_coords[ndims];
  MPI_Cart_coords(up_comm, up_cart_rank, ndims, up_coords);
  MPI_Cart_coords(dn_comm, dn_cart_rank, ndims, dn_coords);
  criticalout dprint(up_coords[0]);
  criticalout dprint(dn_coords[0]);

  //MPI_Comm_dup(up_comm, &dn_comm);
  //MPI_Comm_dup(MPI_COMM_WORLD, &up_comm);
  //MPI_Comm_dup(MPI_COMM_WORLD, &dn_comm);

  // This does not actually do a shift; rather, it returns the
  // ranks of the neighbors that would be used to do this shift.
  // We use this as a mechanism to identify neighbors.
  // Shifting upward means that rank_source will be
  // lower_neighbor and rank_dest will be upper_neighbor.
  int up_src; // from lower
  int up_dst; // to upper
  int dn_src; // from upper
  int dn_dst; // to lower
  MPI_Cart_shift(up_comm, XDIR, UPWARD, &up_src, &up_dst);
  //MPI_Cart_shift(dn_comm, XDIR, DNWARD, &lower_neighbor, &upper_neighbor);
  MPI_Cart_shift(dn_comm, XDIR, DNWARD, &dn_src, &dn_dst);
  criticalout dprint(up_src); // from lower
  criticalout dprint(up_dst); // to upper
  criticalout dprint(dn_src); // from upper
  criticalout dprint(dn_dst); // to lower
  Connection lowXrecvConn(up_comm,up_src);
  Connection hghXrecvConn(dn_comm,dn_src);
  Connection lowXsendConn(dn_comm,dn_dst);
  Connection hghXsendConn(up_comm,up_dst);

  // showing that we can propagate a message upward
  {
    const int count=16;
    char recvbuff[count];
    char sendbuff[count];
    const int sendcount = 1+ // because of terminating '\0'
      snprintf(sendbuff,count,"hello from %d", MPIdata::get_rank());
    criticalout dprint(sendbuff);
    // receive message from lower
    MPI_Request recv_request;
    MPI_Irecv(recvbuff, count, MPI_CHAR, up_src, 0, up_comm, &recv_request);
    // send a message to upper
    MPI_Request send_request;
    MPI_Isend(sendbuff, sendcount, MPI_CHAR, up_dst, 0, up_comm, &send_request);
    // wait for message to arrive
    MPI_Status recv_status;
    MPI_Wait(&recv_request, &recv_status);
    // print the message
    criticalout dprint(recvbuff);
  }

  // showing that upward and downward messages are kept straight
  {
    const int count=16;
    char up_recvbuff[count];
    char up_sendbuff[count];
    char dn_recvbuff[count];
    char dn_sendbuff[count];
    const int up_sendcount = 1+ // because of terminating '\0'
      snprintf(up_sendbuff,count,"up_hello from %d", MPIdata::get_rank());
    const int dn_sendcount = 1+ // because of terminating '\0'
      snprintf(dn_sendbuff,count,"dn_hello from %d", MPIdata::get_rank());
    criticalout dprint(up_sendbuff);
    criticalout dprint(dn_sendbuff);
    // receive message from lower
    MPI_Request up_recv_request;
    MPI_Request dn_recv_request;
    MPI_Irecv(up_recvbuff, count, MPI_CHAR, up_src, 0, up_comm, &up_recv_request);
    MPI_Irecv(dn_recvbuff, count, MPI_CHAR, dn_src, 0, dn_comm, &dn_recv_request);
    // send a message to upper
    MPI_Request up_send_request;
    MPI_Request dn_send_request;
    MPI_Isend(up_sendbuff, up_sendcount, MPI_CHAR, up_dst, 0, up_comm, &up_send_request);
    MPI_Isend(dn_sendbuff, dn_sendcount, MPI_CHAR, dn_dst, 0, dn_comm, &dn_send_request);
    // wait for message to arrive
    MPI_Status up_recv_status;
    MPI_Status dn_recv_status;
    MPI_Wait(&up_recv_request, &up_recv_status);
    MPI_Wait(&dn_recv_request, &dn_recv_status);
    // print the message
    criticalout dprint(up_recvbuff);
    criticalout dprint(dn_recvbuff);
  }

  // showing that we can propagate a particle upward
  {
    const int count=8;
    Particle send_pcl;
    for(int i=0;i<count;i++)
      send_pcl.u[i]=i+.1*MPIdata::get_rank()+.09;
    criticalout debugout << "send_pcl = " << send_pcl << endl;
    // receive message from lower
    Particle recv_pcl;
    MPI_Request recv_request;
    MPI_Irecv(recv_pcl.u, count, MPI_DOUBLE, up_src, 0, up_comm, &recv_request);
    // send a message to upper
    MPI_Request send_request;
    MPI_Isend(send_pcl.u, count, MPI_DOUBLE, up_dst, 0, up_comm, &send_request);
    // wait for message to arrive
    MPI_Status recv_status;
    MPI_Wait(&recv_request, &recv_status);
    // print the message
    criticalout debugout << "recv_pcl = " << recv_pcl << endl;
  }

  // showing that we can propagate a block of particles upward
  {
    masterout debugout << "=== propagating particles upward ===" << endl;

    const int blocksize=8; //
    const int count=8*blocksize;
    //
    // receive message from lower
    //
    Block<Particle> recv_pcls(blocksize, MPIdata::get_rank());
    //vector<Particle>& recv_block = recv_pcls.fetch_block();
    //recv_block.resize(blocksize);
    //MPI_Request recv_request;
    //MPI_Irecv(recv_block[0].u, count, MPI_DOUBLE, up_src, 0, up_comm, &recv_request);
    recv_pcls.recv(lowXrecvConn);
    //
    // send particles
    //
    Block<Particle> send_pcls(blocksize, MPIdata::get_rank());
    vector<Particle>& send_block = send_pcls.fetch_block();
    // initialize particles
    for(int p=0;p<blocksize;p++)
    {
      Particle pcl;
      for(int i=0;i<8;i++)
        pcl.u[i]=p + .1*i+.01*MPIdata::get_rank()+.009;
      send_block.push_back(pcl);
    }
    criticalout debugout << "send_pcls = " << send_pcls << endl;
    // send a message to upper
    //MPI_Request send_request;
    //MPI_Isend(send_block[0].u, count, MPI_DOUBLE, up_dst, 0, up_comm, &send_request);
    send_pcls.send(hghXsendConn);
    // wait for message to arrive
    recv_pcls.waitfor_recv();
    ////MPI_Wait(&recv_request, &recv_status);
    //MPI_Wait(&recv_pcls.fetch_request(),&recv_status);
    // print the message one process at a time
    criticalout debugout << "recv_pcls = " << recv_pcls << endl;
  }

  // showing that we can propagate a list of blocks of particles upward
  if(1)
  {
    masterout debugout << "=== propagating particles upward ===" << endl;

    const int blocksize=4;
    const int numblocks = 2;
    const int numpcls = blocksize*numblocks-1;
    //
    // receive message from lower
    //
    BlockCommunicator<Particle> recv_pcls(blocksize, numblocks, lowXrecvConn);
    recv_pcls.recv_blocks();
    //
    // send particles
    //
    BlockCommunicator<Particle> send_pcls(blocksize, numblocks, hghXsendConn);
    // create and send particles
    int num_blocks_sent = 0;
    for(int p=0;p<numpcls;p++)
    {
      Particle pcl;
      for(int i=0;i<8;i++)
        pcl.u[i]=p + .1*i+.01*MPIdata::get_rank()+.009;
      num_blocks_sent += send_pcls.send(pcl);
    }
    send_pcls.send_complete();
    //assert_eq(num_blocks_sent, numblocks);
    // print all the blocks that were sent.
    for(send_pcls.rewind();!send_pcls.at_end();send_pcls.advance_block())
    {
      Block<Particle>& curr_block = send_pcls.fetch_curr_block();
      dout << curr_block << endl;
    }

    // read and print blocks as they arrive
    {
      recv_pcls.rewind();
      while(!recv_pcls.at_end())
      {
        MPI_Status status;
        MPI_Wait(&recv_pcls.fetch_curr_block().fetch_request(),&status);
        //recv_pcls.waitfor_recv_curr_block(status);
        Block<Particle>& curr_block = recv_pcls.fetch_curr_block();
        recv_pcls.handle_received_block(status);
        dout << curr_block << endl;
      }
    }
  }

  // communicating particles in blocks
  if(1)
  {
    masterout debugout << "=== propagating particles both ways ===" << endl;
    const int blocksize = 8;
    const int numblocks = 2;

    // for each neighbor, create a receive communicator
    //
    BlockCommunicator<Particle> lowXrecv(blocksize, numblocks, lowXrecvConn);
    BlockCommunicator<Particle> hghXrecv(blocksize, numblocks, hghXrecvConn);

    // post receives 
    {
      lowXrecv.recv_blocks();
      hghXrecv.recv_blocks();
    }

    // for each neighbor, create a send communicator
    //
    BlockCommunicator<Particle> lowXsend(blocksize, numblocks, lowXsendConn);
    BlockCommunicator<Particle> hghXsend(blocksize, numblocks, hghXsendConn);

    // create, send, and receive particles
    for(int i=0; i<blocksize*numblocks;i++)
    {
      Particle lowPcl, hghPcl;
      for(int j=0;j<8;j++)
      {
        const double num = 100+i+.1*j;
        lowPcl.u[j] = num+.09;
        hghPcl.u[j] = num+.01;
      }
      hghPcl.u[0] = MPIdata::get_rank();
      lowPcl.u[0] = MPIdata::get_rank();
      //dprint(i)
      bool lowXblockSent = lowXsend.send(lowPcl);
      bool hghXblockSent = hghXsend.send(hghPcl);
      if(lowXblockSent || hghXblockSent)
      {
        // check if particles have arrived,
        // and if so deal with them
        MPI_Status status;
        if(lowXrecv.test_recv_curr_block(status))
        {
          Block<Particle>& curr_block = lowXrecv.fetch_curr_block();
          lowXrecv.handle_received_block(status);
          dout << ": upon sending particle " << i
            << " received low particle block " << curr_block.get_block_id() << ":"
            << curr_block << endl;
        }
        if(hghXrecv.test_recv_curr_block(status))
        {
          Block<Particle>& curr_block = hghXrecv.fetch_curr_block();
          hghXrecv.handle_received_block(status);
          dout << ": upon sending particle " << i
            << " received hgh particle block " << curr_block.get_block_id() << ":"
            << curr_block << endl;
        }
      }
    }

    // send any remaining unsent particles
    //
    lowXsend.send_complete();
    hghXsend.send_complete();

    const int incount=2;
    MPI_Request recv_requests[incount] = 
    {
      lowXrecv.get_curr_request(),
      hghXrecv.get_curr_request()
    };
    // wait on and deal with remaining incoming blocks
    while(!lowXrecv.at_end() || !hghXrecv.at_end())
    {
      // wait for some blocks to be received
      //
      //int recv_indices[incount]; // which requests completed
      //MPI_Status recv_statuses[incount]; // status of completed requests
      //int outcount; // number of requests that returned true
      //MPI_Waitsome(incount, recv_requests, &outcount, recv_indices, recv_statuses);

      //#if 0
      int recv_index;
      MPI_Status recv_status;
      MPI_Waitany(incount, recv_requests, &recv_index, &recv_status);
      switch(recv_index)
      {
        default:
          invalid_value_error(recv_index);
        case MPI_UNDEFINED:
          eprintf("recv_requests contains no active handles");
          break;
        case 0: // lowXrecv
         {
          Block<Particle>& curr_block = lowXrecv.fetch_curr_block();
          lowXrecv.handle_received_block(recv_status);
          dout << ": received lowXrecv."
            << curr_block.get_block_id()
            << curr_block << endl;
          recv_requests[0] = lowXrecv.get_curr_request();
         }
          break;
        case 1: // hghXrecv
         {
          Block<Particle>& curr_block = hghXrecv.fetch_curr_block();
          hghXrecv.handle_received_block(recv_status);
          dout << ": received hghXrecv."
            << curr_block.get_block_id()
            << curr_block << endl;
          recv_requests[1] = hghXrecv.get_curr_request();
         }
          break;
      }
      //#endif

      MPI_Status status;
      if(!lowXrecv.at_end() && lowXrecv.test_recv_curr_block(status))
      {
        Block<Particle>& curr_block = lowXrecv.fetch_curr_block();
        lowXrecv.handle_received_block(status);
        dout << "(" << MPIdata::get_rank() << ") line "
          << __LINE__ << ": received lowXrecv."
          << curr_block.get_block_id() << ":"
          << curr_block << endl;
      }
      if(!hghXrecv.at_end() && hghXrecv.test_recv_curr_block(status))
      {
        Block<Particle>& curr_block = hghXrecv.fetch_curr_block();
        hghXrecv.handle_received_block(status);
        dout << "(" << MPIdata::get_rank() << ") line "
          << __LINE__ << ": received hghXrecv."
          << curr_block.get_block_id() << ":"
          << curr_block << endl;
      }
    }
  }

  // probably MPI_Finalize takes care of this anyway...
  MPI_Comm_free(&up_comm);
  MPI_Comm_free(&dn_comm);
  MPIdata::finalize_mpi();
  //MPI_Finalize();
}
