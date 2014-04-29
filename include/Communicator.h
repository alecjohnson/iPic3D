
template<class type>
class Ring
{
  std::list<Larray<type> >;
}

// class to communicate equal-sized messages
// in batches between MPI processes

template<class type>
class BatchCommunicator
{
 private:
   // data structure to handle messages
  //Ring<Larray<type> > ring;
  stl::list<vector<type> > bufferList;

 public:
  // append msg to send buffer
  void send(const type& msg)
  {
    
  }
};
