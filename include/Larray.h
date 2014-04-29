#ifndef Larray_h
#define Larray_h

#include "Alloc.h" // for ALIGNED
#include "ipicmath.h" // for pow2roundup
#include "asserts.h"
#include "string.h" // for memcpy

// linear array (e.g. of particles)
// 
// this is pretty much the same thing as stl::vector
//
template<class type>
class Larray
{
 private: // members
  type* list;
  int _size; // number of particles in list
  int _maxsize; // maximum number of particles
 private:
  void check_index(int i)
  {
    #ifdef CHECK_BOUNDS
      assert_le(0, i);
      assert_lt(i, _size);
    #endif
  }
 public: // access
  int size()const { return _size; }
  int maxsize()const { return _maxsize; }
  void pop()
  {
    assert(_size>0);
    _size--;
  }
  void clear()
  {
    _size = 0;
  }
  void push_back(const type& element)
  {
    if(_size>=_maxsize)
    {
      int newmaxsize = pow2roundup(_size+1);
      realloc(newmaxsize);
    }
    list[_size] = element;
    _size++;
  }
  inline type& operator[](int i)
  {
    check_index(i);
    ALIGNED(list);
    return list[i];
  }
  void delete_element(int i)
  {
    // works even for last particle,
    // though pop would be faster in that case.
    // 
    // why doesn't this compile?
    //this->operator[i] = list[--_size];
    check_index(i);
    ALIGNED(list);
    list[i] = list[--_size];
  }
 public: // memory
  ~Larray()
  {
    AlignedFree(list);
  }
  Larray():
    list(0),
    _size(0),
    _maxsize(0)
  {}
  Larray(int requested_size)
    list(0),
    _size(0),
    _maxsize(0)
  { if(requested_size > 0) realloc(requested_size); }
  // reset maximum size without deleting elements
  void realloc(int newmaxsize)
  {
    // ignore request if requested size is too small
    if(_size > newmaxsize) return;

    // round up size to a multiple of num_elem_in_block
    const int num_elem_in_block = 8;
    newmaxsize = ((newmaxsize-1)/num_elem_in_block+1)*num_elem_in_block;
    if(newmaxsize != _maxsize)
    {
      _maxsize = newmaxsize;
      type* oldList = list;
      // the next two lines assume that type has no indirection
      list = AlignedAlloc(type,_maxsize);
      memcpy(list,oldList,sizeof(type)*_size);
      // for(int i=0;i<_size;i++) list[i] = oldList[i];
      AlignedFree(oldList);
    }
  }
  void realloc_if_smaller_than(int required_max_size)
  {
    if(_size < required_max_size)
      realloc(required_max_size);
  }
  void shrink()
  {
    // shrink _maxsize by a factor of two if elements will fit.
    int proposed_size = pow2rounddown(_maxsize/2);
    if( _size <= proposed_size && proposed_size < _maxsize)
    {
      realloc(proposed_size);
    }
  }
};

#endif
