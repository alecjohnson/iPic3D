#ifndef ipic_arrays_h
#define ipic_arrays_h
#include "Alloc.h"
#include <new> // needed for placement new

// should make a const version of this class
class vector_array3_double
{
    const int numel;
    array3_double*list;
  public:
    ~vector_array3_double()
    {
      for(int i=0;i<numel;i++)
      {
        list[i].~array3_double();
      }
      free(list);
    }
    vector_array3_double(int numel_, size_t nx, size_t ny, size_t nz):
      numel(numel_),
      list(0)
    {
      list = (array3_double*)malloc(numel*sizeof(array3_double));
      for(int i=0;i<numel;i++)
      {
        // placement new
        new(list+i) array3_double(nx,ny,nz);
      }
    }
    array3_double& operator[](int i)
    {
      #ifdef CHECK_BOUNDS
      assert_le(0,i);
      assert_lt(i,numel);
      #endif // CHECK_BOUNDS
      return list[i];
    }
};
#endif // ipic_arrays_h
