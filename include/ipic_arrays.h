#ifndef ipic_arrays_h
#define ipic_arrays_h
#include "Alloc.h"
#include <new> // needed for placement new

// with a properly written array class this class would not
// be needed.  Could generalize this to make arr3_double a
// template parameter.

// Inheritance chain:
// const_vector_arr3_double
// -> ref_vector_arr3_double
// -> vector_arr3_double
class const_vector_arr3_double
{
  protected:
    const int numel;
    arr3_double*list;
  public:
    void free()
    {
      for(int i=0;i<numel;i++)
      {
        list[i].free();
      }
      ::free(list);
    }
    const_vector_arr3_double(int numel_, size_t nx, size_t ny, size_t nz):
      numel(numel_),
      list(0)
    {
      list = (arr3_double*)malloc(numel*sizeof(arr3_double));
      for(int i=0;i<numel;i++)
      {
        // placement new
        new(list+i) arr3_double(nx,ny,nz);
      }
    }
    const_arr3_double operator[](int i)
    {
      #ifdef CHECK_BOUNDS
      assert_le(0,i);
      assert_lt(i,numel);
      #endif // CHECK_BOUNDS
      return list[i];
    }
};
class ref_vector_arr3_double : public const_vector_arr3_double
{
  public:
    ref_vector_arr3_double(int numel_, size_t nx, size_t ny, size_t nz):
      const_vector_arr3_double(numel_, nx, ny, nz)
    {}
    arr3_double operator[](int i)
    {
      #ifdef CHECK_BOUNDS
      assert_le(0,i);
      assert_lt(i,numel);
      #endif // CHECK_BOUNDS
      return list[i];
    }
};
class vector_arr3_double : public ref_vector_arr3_double
{
  public:
    vector_arr3_double(int numel_, size_t nx, size_t ny, size_t nz):
      ref_vector_arr3_double(numel_, nx, ny, nz)
    {}
    ~vector_arr3_double(){free();}
};
#endif // ipic_arrays_h
