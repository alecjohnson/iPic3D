#ifndef BlockCommunicator_h
#define BlockCommunicator_h
#include "aligned_allocator.h"
#include "asserts.h"
#include "debug.h" // temporary
#include <vector> // needed for aligned_vector
#include <list>
#include <assert.h>

// The combination of group (comm), tag, and neighbor
// should be unique for each connection.
//
// In a Cartesian topology, the challenges to this uniqueness
// occur for periodic boundary conditions when there is a
// dimension only 1 or 2 processes thick.
//
// In the 1-process-thick case, we can forgo MPI communication in
// that direction altogether.  But in the 2-process-thick case,
// we need to make a distinction between upward and downward
// channels of communication.  Distinguishing upward and downward
// directions is sufficient for any topology based on convex MPI
// subdomains, since in this case two subdomains can share at
// most two faces.  Any mechanism for distinguishing upward and
// downward directions for a pair of subdomains is sufficient;
// one way is to compare the positions of the centroids, first in
// the x direction, then in the y direction, and then in the z
// direction.  One might use tag=0 for downward communication and
// tag=1 for upward communication, and use MPI_COMM_WORLD for the
// group.
//
class Connection
{
 public: // constants
  enum Tag
  {
    DEFAULT = 0,
    PARTICLE_DN, // downward communication of particles
    PARTICLE_UP // upward communication of particles
  };
 private: // data
  // In MPI a message envelope includes the following
  // information plus the rank of this process.
  int _rank; // rank within group of neighbors we're connecting to
  int _tag; // tag to attach to messages
  MPI_Comm _comm; // communicator group
 public: // init
  Connection(int rank_, int tag_, const MPI_Comm& comm_):
    _rank(rank_),
    _tag(tag_),
    _comm(comm_)
  {}
 public: // accessors
  int rank()const{return _rank;}
  int tag()const{return _tag;}
  MPI_Comm comm()const{return _comm;}
};

bool signal_hack()
{
  return false;
}

// block of elements
template <class type>
struct Block
{
 private:
 private: // initialized in constructor
  // change the block to use an aligned allocator
  aligned_vector(type) block;
  int capacity;
  int id;
  // used for MPI communication
  MPI_Request request;
  int flag;
  // hack to piggy-back information onto message
  int signal; // char
 private: // initialized at compile time
  // assumes using MPI_DOUBLE
  static const int NUMBERS_PER_ELEMENT = sizeof(type)/sizeof(double);
 public:
  Block(int capacity_, int id_):
    capacity(capacity_),
    id(id_),
    request(MPI_REQUEST_NULL),
    flag(0),
    signal(0)
  {
    block.reserve(capacity);
  }
 public: // accessors
  MPI_Request& fetch_request(){return request;}
  aligned_vector(type)& fetch_block(){return block;}
  const aligned_vector(type)& get_block()const{return block;}
  int get_block_id(){ return id; }
  void set_finished_flag() { signal |= 2; }
  bool finished_flag_is_set() { return signal & 2; }
  void set_insert_flag() { signal |= 1; }
  bool insert_flag_is_set() { return signal & 1; }
  int size(){return block.size();}

 public: // operations
  bool test_comm(MPI_Status& status)
  {
    MPI_Test(&request, &flag, &status);
    if(!flag)
      return false;
    // MPI_Test man page says this should now be true
    assert(request==MPI_REQUEST_NULL);
    return true;
  }
  bool test_comm()
  {
    MPI_Status status;
    return test_comm(status);
  }

 // sending
 //
 public: // sending
  type& operator[](int i){return block[i];}
  bool is_active()
  {
    return request!=MPI_REQUEST_NULL;
  }
  bool is_finished()
  {
    return request==MPI_REQUEST_NULL;
  }
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
    // hack: put signal in an extra "particle"
    if(signal_hack())
    {
      type signal_element;
      // convert signal to double and put in signal_element
      *((double*) &signal_element) = double(signal);
      push_back(signal_element);
    }

    dprintf("sending block number %d", id);
    MPI_Isend(&block[0], NUMBERS_PER_ELEMENT*block.size(), MPI_DOUBLE,
      dest.rank(), dest.tag(), dest.comm(), &request);
  }
  void waitfor_send()
  {
    MPI_Status status;
    MPI_Wait(&request, &status);
  }
  void set_request(MPI_Request request_)
  {
    request = request_;
  }
  void clear()
  {
    block.clear();
    request = MPI_REQUEST_NULL;
  }
 public: // receiving
  // post a receive
  void recv(const Connection& source)
  {
    // make sure that space exists to receive
    int newsize = signal_hack() ? capacity+1 : capacity;
    block.resize(newsize);
    assert(request==MPI_REQUEST_NULL);
    MPI_Irecv(&block[0], NUMBERS_PER_ELEMENT*block.size(), MPI_DOUBLE,
      source.rank(), source.tag(), source.comm(), &request);
  }
  // processing received data
  //
  // returns true if number of elements received
  // was smaller than capacity
  bool shrink_received_block(MPI_Status& status)
  {
    int count;
    MPI_Get_count(&status, MPI_DOUBLE, &count);
    int num_elements_received = count / NUMBERS_PER_ELEMENT;
    // hack: handle the signal particle
    if(signal_hack())
    {
      // discount signal particle
      num_elements_received--;
      // extract signal from signal element
      signal = int(*((double*)&block[num_elements_received]));
    }
    else
    {
      if(num_elements_received < capacity)
      {
        set_finished_flag();
      }
    }
    bool retval = (num_elements_received < capacity);
    assert_le(num_elements_received,capacity);
    // shrink to eliminate the trailing garbage
    block.resize(num_elements_received);
    return retval;
  }
  void waitfor_recv()
  {
    MPI_Status status;
    MPI_Wait(&request, &status);
    shrink_received_block(status);
  }
  // returns true if message has been received.
};

// We use the same class for sending and receiving, because
// in the case where a dimension is only one process thick,
// the sender and receiver are the same process (in which case
// the send and receive methods do not need to be called at
// all).  A more correct way to handle this would be to inherit
// from BlockSender_interface and BlockRecver_interface abstract
// base classes.
//
// The implementation should probably be changed to use
// persistent requests.  My analysis of the decision 
// regarding message size is as follows:
//
#if 0
   Communication time is a linear combination of total data
   to be sent and the number of messages it is broken into:

     time = data_time + message_overhead, i.e.,
     time = data_size*data_coef + num_msgs*msg_coef.

   Assuming no blocking issues,
   data_time is presumably the same for both persistent
   and non-persistent communication.

   message_overhead is presumably much smaller for persistent
   than non-persistent communication (or persistent would not add
   value).

   data_time must dominate message_overhead for relevant cases or
   persistent communication adds no value.

   Assuming persistent communication can add value,
   messages should be large enough so that data_time dominates
   for persistent communication but small enough so that
   message_overhead dominates for non-persistent communication.
   This determines a window of message sizes appropriate
   for persistent communication.
#endif
//
template <class type>
class BlockCommunicator
{
 enum CommState
 {
   INITIAL=0,
   ACTIVE,
   FINISHED
 };
 private: // data
  std::list<void*>::iterator curr_block;
  // list of pointers avoids need for placement new
  std::list<void*> blockList;
  //std::list<Block<type>*>::iterator curr_block;
  //std::list<Block<type>*> blockList;
  Connection connection;
  int blocksize;
  int nextid;
  int commState;
 // methods
 public: // construction
  BlockCommunicator(int blocksize_, int numblocks, Connection connection_):
    blocksize(blocksize_),
    connection(connection_),
    commState(INITIAL)
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
    std::list<void*>::iterator i;
    for(i=blockList.begin(); i != blockList.end(); ++i)
      delete &fetch_block(i);
  }
 // access
 private:
  Block<type>& fetch_block(std::list<void*>::iterator block)
  {
    return *((Block<type>*)*block);
  }
 public: // information
  // access current block
  Block<type>& fetch_curr_block()
  {
    assert(curr_block != blockList.end());
    return fetch_block(curr_block);
  }
  bool comm_finished() const
  {
    return commState==FINISHED;
  }

 // operations
 //
 private:
  // mimic behavior of a ring
  void increment_block(list<void*>iterator& block_iter)
  {
    block_iter++;
    if(block_iter == blockList.end())
    block_iter = blockList.begin();
  }
  void increment_curr_block()
  {
    increment_block(curr_block)
  }
  //void rewind()
  //{
  //  curr_block = blockList.begin();
  //}
 public: // receiving operations

  // insert blocks prior to the current block
  // and call recv on them.
  //
  // If there is nowhere to write an incoming message, then I
  // suppose that a good MPI implementation would buffer it
  // somewhere until we post a receive that can accept it, rather
  // than blocking on the completing the send.  If MPI provides a
  // mechanism to query how much such intermediate buffering is
  // happening then we could use it as a trigger to insert more
  // receive blocks, but I don't know of such a mechanism.
  //
  // Alternatively, I think that we can take the number of send
  // blocks needed to keep from having to wait for sends to
  // complete as a reasonable indicator of the number of receive
  // blocks needed to avoid having to buffer.  So when the sender
  // adds more blocks, it can signal the receiver to do likewise.
  //
  void insert_more_recv_blocks(int numblocks=1)
  {
    for(int i=0;i<numblocks;i++)
    {
      Block<type>* newBlock = new Block<type>(blocksize, nextid++);
      blockList.insert(curr_block, newBlock);
      newBlock->recv(connection);
    }
  }
 public: // receiving operations
  //
  // if this returns true, then the user must follow up by
  // calling fetch_received_block(status);
  //
  bool test_recv_curr_block(MPI_Status& status)
  {
    return fetch_curr_block().test_comm(status);
  }
  // After the user has processed the block returned by this
  // call, he should call release_received_block() to
  // advance to the next block and free this block to
  // be used for receiving another message.
  //
  Block<type>& fetch_received_block(MPI_Status& status)
  {
    bool shrunk = fetch_curr_block().shrink_received_block(status);
    fetch_curr_block().set_request(MPI_REQUEST_NULL);
    if(signal_hack())
    {
      if(fetch_curr_block().insert_flag_is_set())
      {
        // insert another block to receive particles
        // prior to the block that we just received
        dprintf("as requested, I am inserting a block.");
        insert_more_recv_blocks();
      }
    }
    if(fetch_curr_block().finished_flag_is_set())
    {
      dprintf("looks like that was the last block.");
      commState=FINISHED;
    }
    return fetch_curr_block();
  }
  // reactivate recv request and advance to next block
  void release_received_block()
  {
    // post receive again and proceed to next block
    fetch_curr_block().recv(connection);
    increment_curr_block();
  }

  // cancel and free any open receive requests
  void cancel_recvs()
  {
    std::list<void*>::iterator b;
    for(b=blockList.begin(); b!=blockList.end();b++)
    {
      MPI_Request& pending_request = fetch_curr_block().fetch_request();
      if(pending_request!=MPI_REQUEST_NULL)
      {
        MPI_Cancel(&pending_request);
        MPI_Request_free(&pending_request);
      }
    }
  }

  MPI_Request get_curr_request()
  {
    return fetch_curr_block().fetch_request();
  }
 public:
  bool recv_blocks()
  {
    // post receives on all blocks
    std::list<void*>::iterator b;
    for(b=blockList.begin(); b != blockList.end(); ++b)
      fetch_block(b).recv(connection);
    // reset curr_block to initial block
    curr_block = blockList.begin();
  }
  //void waitfor_recv_curr_block(MPI_Status& status)
  //{
  //  fetch_curr_block().waitfor_recv(status);
  //}

 // sending routines
 //
 private:
  // assumes that there is a block to send
  void send_block()
  {
    dprint(fetch_curr_block().size());
    std::list<void*>::iterator block_to_send = curr_block;

    std::list<void*>::iterator next_block = curr_block;
    increment_block(next_block);

    // proceed to the next block
    if(!signal_hack())
    {
      // if next block is still sending, wait for it to finish
      // sending so that we can reuse it.
      if(fetch_curr_block().is_active())
      {
        dprintf("waiting for send of block %d to complete.", fetch_curr_block().get_id());
        fetch_curr_block().waitfor_send();
        dprintf("done waiting for send to complete.");
        assert(!fetch_curr_block().is_active());
        assert(fetch_curr_block().test_comm());
      }
    }

    // if next block is still sending, insert another block for
    // sending, use it instead, and tell the receiver to insert
    // another block for receiving.
    //
    if(signal_hack() && (!fetch_curr_block().test_comm()))
    {
      Block<type>* newBlock = new Block<type>(blocksize, nextid);
      blockList.insert(curr_block,newBlock);
      curr_block = block_to_send;
      increment_curr_block();
      // the new block is not "still sending"
      assert(fetch_curr_block().test_comm());

      // hack: set flag in block to tell receiver
      // to insert another block for receiving
      fetch_block(block_to_send).set_insert_flag();
    }

    // send the block
    dprint(fetch_block(block_to_send).size());
    fetch_block(block_to_send).send(connection);
    // clear the new block
    fetch_curr_block().clear();
  }
 public: // interface
  // send a single particle.
  //
  // returns true iff block was sent
  // (whereupon the caller might want to check for incoming messages).
  bool send(const type& in)
  {
    // append the particle to the block.
    fetch_curr_block().push_back(in);
    dprint(fetch_curr_block().size());
    dprint(fetch_curr_block().isfull());

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
  // this assumes that curr_block exists and contains all unsent particles
  void send_complete()
  {
    assert(!fetch_curr_block().isfull());
    fetch_curr_block().set_finished_flag();
    send_block();
    commState=FINISHED;
  }
  //void clear_send()
  //{
  //  // make sure that pending sends have cleared
  //  //
  //  std::list<void*>::iterator b = blockList.begin();
  //  for(;b!=blockList.end();b++)
  //  {
  //    // might be more efficient to make
  //    // a single call to MPI_Waitall
  //    fetch_block(b).waitfor_send();
  //    // clear the block
  //    fetch_block(b).clear();
  //  }
  //  // reset curr_block
  //  rewind();
  //}
 public:
};

//#include <iosfwd> // for ostream
//
//inline std::ostream& operator<<(std::ostream& os, const Block<Particle>& block_)
//{
//  const aligned_vector(Particle)& block = block_.get_block();
//  for(int k=0; k<block.size();k++)
//  {
//    os << "\n  block[" << k << "] = " << block[k];
//  }
//  return os;
//}

#endif
