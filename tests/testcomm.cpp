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

// particle ("element")
struct Particle
{
  double u[8];
};

std::ostream& operator<<(std::ostream& os, const Particle& pcl)
{
  //os << showpos;
  for(int i=0; i<7; i++)
    os << showpos << pcl.u[i] << " ";
  os << pcl.u[8] << "\n";
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
  MPI_Request request;
  int flag;
  int id;
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
  MPI_Request& fetch_request(){return request;}
  //int size()
  //{
  //  return block.size();
  //}
  vector<type>& fetch_block(){return block;}
  const vector<type>& get_block()const{return block;}
  int get_block_id(){ return id; }
  bool isfull()
  {
    return (block.size()>=capacity);
  }
  void push_back(const type& in)
  {
    block.push_back(in);
  }
  void send(const Connection& dest, int tag)
  {
    assert(block.size() <= capacity);
    MPI_Isend(&block[0], NUMBERS_PER_ELEMENT*block.size(), MPI_DOUBLE,
      dest.rank(), tag, dest.comm(), &request);
  }
  void recv(const Connection& source)
  {
    // make sure that space exists to receive
    block.resize(capacity);
    int tag=0; // what on earth is this used for?
    MPI_Irecv(&block[0], NUMBERS_PER_ELEMENT*block.size(), MPI_DOUBLE,
      source.rank(), tag, source.comm(), &request);
  }
  // returns true if message has been received.
  bool test_recv(int& count)
  {
    MPI_Status status;
    bool ret = test(status);
    if(!ret)
      return false;
    MPI_Get_count(&status, MPI_DOUBLE, &count);
    assert(count<=capacity);
    // shrink to eliminate the trailing garbage
    block.resize(count);
    return true;
  }
  bool test_recv()
  {
    int count;
    return test_recv(count);
  }
  // returns true if communication is complete
  bool test(MPI_Status& status)
  {
    MPI_Test(&request, &flag, &status);
    int count;
    MPI_Get_count(&status, MPI_DOUBLE, &count);
    if(flag) // communication complete
    {
      assert(request==MPI_REQUEST_NULL); // based on MPI_Test man page
      return true;
    }
    return false;
  }
};

std::ostream& operator<<(std::ostream& os, const Block<Particle>& block_)
{
  const vector<Particle>& block = block_.get_block();
  for(int k=0; k<block.size();k++)
  {
    os << "  block[" << k << "] = " << block[k] << "\n";
  }
  return os;
}

list<Block<Particle>*>::iterator blah;

// We use the same class for sending and receiving, because
// in the case where a dimension is only one process thick,
// the sender and receiver are the same process (in which case
// the send and receive methods do not need to be called at
// all).  A more correct way to handle this would be to inherit
// from BlockSender_interface and BlockRecver_interface abstract
// base classes.
//
template <class type>
struct BlockCommunicator
{
 private: // data
  list<void*>::iterator curr_block;
  list<void*> blockList;
  //list<Block<type>*>::iterator curr_block;
  //list<Block<type>*> blockList;
  Connection connection;
  int blocksize;
  int nextid;
 protected: // methods
  void allocate_more_blocks()
  {
    // allocate more blocks and insert them
    //
    Block<type>* newBlock = new Block<type>(blocksize, nextid++);
    blockList.push_back(newBlock);
  }
 public: // interface
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
    //list<Block<type>*>::iterator i;
    for(curr_block=blockList.begin(); curr_block != blockList.end(); ++curr_block)
      delete &fetch_curr_block();
  }

  Block<type>& fetch_curr_block()
  {
    assert(curr_block != blockList.end());
    return *((Block<type>*)*curr_block);
  }
  bool atend()
  {
    return curr_block == blockList.end();
  }
  void advance_to_next_block()
  {
    assert(curr_block != blockList.end());
    curr_block++;
  }
  bool test_recv_curr_block()
  {
    return fetch_curr_block().test_recv();
  }
  bool test_curr_block()
  {
    return fetch_curr_block().test();
  }
  MPI_Request get_curr_request()
  {
    if(atend())
      return MPI_REQUEST_NULL;
    return fetch_curr_block().fetch_request();
  }

 // sending routines
 //
 private: // methods
  void send_block(int block_tag)
  {
    // send the block
    fetch_curr_block().send(connection, block_tag);

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
    if(fetch_curr_block().isfull());
    {
      int block_tag = Comm::UNFINISHED;
      send_block(block_tag);
      return true;
    }
    return false;
  }
 
 public: // receiving methods
  // post receive on all blocks
  // and reset 
  bool recv_blocks()
  {
    // post receives on all blocks
    //list<Block<type>*>::iterator i;
    for(curr_block=blockList.begin(); curr_block != blockList.end(); ++curr_block)
      fetch_curr_block().recv(connection);
    // reset curr_block to initial block
    curr_block = blockList.begin();
  }
};

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
  //MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 0, &dn_comm);
  MPI_Comm_dup(up_comm, &dn_comm);
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
  Connection lowXrecvConn(up_comm,up_src);
  Connection hghXrecvConn(dn_comm,dn_src);
  Connection lowXsendConn(dn_comm,dn_dst);
  Connection hghXsendConn(up_comm,up_dst);

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
      lowPcl.u[j] = i+.1*(j+1);
      hghPcl.u[j] = -lowPcl.u[j];
    }
    hghPcl.u[0] = MPIdata::get_rank();
    lowPcl.u[0] = MPIdata::get_rank();
    hghPcl.u[1] = MPIdata::get_rank();
    lowPcl.u[1] = MPIdata::get_rank();
    bool lowXblockSent = lowXsend.send(lowPcl);
    bool hghXblockSent = hghXsend.send(hghPcl);
    if(lowXblockSent || hghXblockSent)
    {
      // check if particles have arrived,
      // and if so deal with them
      if(lowXrecv.test_recv_curr_block())
      {
        Block<Particle>& curr_block = lowXrecv.fetch_curr_block();
        cout << "(" << MPIdata::get_rank() << ") line "
          << __LINE__ << ": while sending block " << i
          << " received particle block" << curr_block.get_block_id() << "\n"
          << curr_block << endl;
        lowXrecv.advance_to_next_block();
      }
      if(hghXrecv.test_recv_curr_block())
      {
        Block<Particle>& curr_block = hghXrecv.fetch_curr_block();
        cout << "(" << MPIdata::get_rank() << ") line "
          << __LINE__ << ": while sending block " << i
          << " received particle block" << curr_block.get_block_id() << "\n"
          << curr_block << endl;
      }
    }
  }
  const int incount=2;
  MPI_Request recv_requests[incount] = 
  {
    lowXrecv.get_curr_request(),
    hghXrecv.get_curr_request()
  };
  // wait on and deal with remaining incoming blocks
  while(!lowXrecv.atend() || !hghXrecv.atend())
  {
    // wait for some blocks to be received
    //
    //int recv_indices[incount]; // which requests completed
    //MPI_Status recv_statuses[incount]; // status of completed requests
    //int outcount; // number of requests that returned true
    //MPI_Waitsome(incount, recv_requests, &outcount, recv_indices, recv_statuses);

    int recv_index;
    MPI_Status recv_status;
    MPI_Waitany(incount, recv_requests, &recv_index, &recv_status);
    switch(recv_index)
    {
      default:
        invalid_value_error(recv_index)
      case MPI_UNDEFINED:
        eprintf("recv_requests contains no active handles");
        break;
      case 0: // lowXrecv
       {
        Block<Particle>& curr_block = lowXrecv.fetch_curr_block();
        bool test = curr_block.test_recv(); // resize the list
        assert(test);
        cout << "(" << MPIdata::get_rank() << ") line "
          << __LINE__ << ": received lowXrecv."
          << curr_block.get_block_id() << "\n"
          << curr_block << endl;
        lowXrecv.advance_to_next_block();
        recv_requests[0] = lowXrecv.get_curr_request();
       }
        break;
      case 1: // hghXrecv
       {
        Block<Particle>& curr_block = hghXrecv.fetch_curr_block();
        bool test = curr_block.test_recv(); // resize the list
        assert(test);
        cout << "(" << MPIdata::get_rank() << ") line "
          << __LINE__ << ": received hghXrecv."
          << curr_block.get_block_id() << "\n"
          << curr_block << endl;
        hghXrecv.advance_to_next_block();
        recv_requests[1] = hghXrecv.get_curr_request();
       }
        break;
    }

    #if 0
    if(!lowXrecv.atend() && lowXrecv.test_curr_block())
    {
      Block<Particle>& curr_block = lowXrecv.fetch_curr_block();
      cout << "(" << MPIdata::get_rank() << ") line "
        << __LINE__ << ": received lowXrecv."
        << block.get_block_id() << "\n"
        << block << endl;
      lowXrecv.advance_to_next_block();
    }
    if(!hghXrecv.atend() && hghXrecv.test_curr_block())
    {
      Block<Particle>& curr_block = hghXrecv.fetch_curr_block();
      cout << "(" << MPIdata::get_rank() << ") line "
        << __LINE__ << ": received hghXrecv."
        << block.get_block_id() << "\n"
        << block << endl;
      hghXrecv.advance_to_next_block();
    }
    #endif
  }

  // probably MPI_Finalize takes care of this anyway...
  MPI_Comm_free(&up_comm);
  MPI_Comm_free(&dn_comm);
  MPIdata::finalize_mpi();
  //MPI_Finalize();
}
