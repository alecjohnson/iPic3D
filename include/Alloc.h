#ifndef IPIC_ALLOC_H
#define IPIC_ALLOC_H
#include <cstddef> // for alignment stuff
#include "ipicdefs.h" // for CHECK_BOUNDS
#include "asserts.h" // for assert_le, assert_lt

/*
    Array class devloped by Alec Johnson
      consolidating arrays developed by 
    Reger Ferrer, Vicenç Beltran, and Florentino Sainz
      and earlier arrays written by
    Jorge Amaya and Stefano Markidis.
*/
#define ALIGNMENT (64)
#ifdef __INTEL_COMPILER
    #define ALIGNED(X) __assume_aligned(X, ALIGNMENT)
    #define AlignedAlloc(T, NUM) \
        (T *const __restrict__)(_mm_malloc(sizeof(T)*NUM, ALIGNMENT))
    #define AlignedFree(S) (_mm_free(S))
#else
    #define ALIGNED(X)
    #define AlignedFree(S) (delete[] S)
    #define AlignedAlloc(T, NUM) (new T[NUM]) 
#endif

/****** begin Alec Johnson's code ******/

// methods to allocate arrays
//
template < class type > inline type * newArray1(int sz1) {
  type *arr = AlignedAlloc(type, sz1); // new type [sz1];
  return arr;
}
template < class type > inline type ** newArray2(int sz1, int sz2) {
  type **arr = AlignedAlloc(type*, sz1); // new type *[sz1];
  type *ptr = newArray1<type>(sz1*sz2);
  for (int i = 0; i < sz1; i++)
  {
    arr[i] = ptr;
    ptr += sz2;
  }
  return arr;
}
template < class type > inline type *** newArray3(int sz1, int sz2, int sz3) {
  // type ***arr = AlignedAlloc(type**, sz1); // new type **[sz1];
  type ***arr = new type **[sz1];
  type **ptr = newArray2<type>(sz1*sz2, sz3);
  for (int i = 0; i < sz1; i++)
  {
    arr[i] = ptr;
    ptr += sz2;
  }
  return arr;
}
template < class type > inline type **** newArray4(int sz1, int sz2, int sz3, int sz4) {
  type ****arr = AlignedAlloc(type***, sz1); //(new type ***[sz1]);
  type ***ptr = newArray3<type>(sz1*sz2, sz3, sz4);
  for (int i = 0; i < sz1; i++) {
    arr[i] = ptr;
    ptr += sz2;
  }
  return arr;
}

// methods to deallocate arrays
//
template < class type > inline void delArr1(type * arr)
{
  AlignedFree(arr);
}
template < class type > inline void delArr2(type ** arr, int sz1)
{
  delArr1(arr[0]);
  AlignedFree(arr);
}
template < class type > inline void delArr3(type *** arr, int sz1, int sz2)
{
  delArr2(arr[0],sz2);
  AlignedFree(arr);
}
template < class type > inline void delArr4(type **** arr,
  int sz1, int sz2, int sz3)
{
  delArr3(arr[0],sz2,sz3);
  AlignedFree(arr);
}


// classes to dereference arrays.
//
// The purpose of this class is to allow elements of
// multidimensional arrays to be accessed with a calculated
// one-dimensional index while using the same syntax as is used
// for a nested array.  This gives correct results, but is slow
// unless optimization is turned on, allowing the compiler to
// figure out that the whole chain of calls to the operator[]
// methods and to the RefN constructors reduces to computing
// a one-dimensional subscript used to access a one-dimensional
// array.
//
template <class type>
class Ref1
{
  type* __restrict__ const arr;
  const int* const sizes;
  int k;
 public:
  inline Ref1(type *const arr_in, int const* const sizes_in, int k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline type& operator[](int idx){
    #ifdef CHECK_BOUNDS
    assert_le(0, idx); assert_lt(idx, sizes[0]);
    #endif
    k *= sizes[0]; k += idx;
    ALIGNED(arr);
    return arr[k];
  }
};

template <class type>
class Ref2
{
  type* __restrict__ const arr;
  const int* const sizes;
  int k;
 public:
  inline Ref2(type *const arr_in, int const*const sizes_in, int k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline Ref1<type> operator[](int idx){
    #ifdef CHECK_BOUNDS
    assert_le(0, idx); assert_lt(idx, sizes[0]);
    #endif
    k *= sizes[0]; k += idx;
    return Ref1<type>(arr, sizes+1, k);
  }
};

template <class type>
class Ref3
{
  type* __restrict__ const arr;
  const int* const sizes;
  int k;
 public:
  inline Ref3(type*const arr_in, int const*const sizes_in, int k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline Ref2<type> operator[](int idx){
    #ifdef CHECK_BOUNDS
    assert_le(0, idx); assert_lt(idx, sizes[0]);
    #endif
    k *= sizes[0]; k += idx;
    return Ref2<type>(arr, sizes+1, k);
  }
};

// ArrN adopts an array allocated by newArrN
//
// The purpose of these classes is to provide more efficient
// and more regulated access to array elements.  The idea is to
// maintain backward compatibility while allowing us to move
// toward a proper array abstraction.
//
// proposed improvements:
// - alternative constructor that takes a list of array dimensions
// - reference counting and destructor that deallocates if appropriate
// - methods that use parallel arithmetic for omp and vectorized code

template <class type>
class Arr1
{
  private: // data
    type* __restrict__ arr;
    int sizes[3];
    bool owner;
  private: // methods
    void set_size(int s1)
    {
      sizes[0] = 1;
      sizes[1] = s1;
      sizes[2] = 0;
    }
  public:
    ~Arr1()
    {
      if(owner) AlignedFree(arr);
    }
    Arr1(int s1) :
      owner(true)
    {
      set_size(s1);
      arr = AlignedAlloc(type, s1);
    }
    inline Arr1(type* in, int s1) :
      arr(in),
      owner(false)
    {
      set_size(s1);
    }
    inline type operator[](int idx){
      #ifdef CHECK_BOUNDS
      assert_le(0, idx); assert_lt(idx, sizes[1]);
      #endif
      ALIGNED(arr);
      return arr[idx];
    }
};

template <class type>
class Arr2
{
  private: // data
    type* __restrict__ arr;
    int sizes[4];
    bool owner;
  private: // methods
    void set_sizes(int s1, int s2)
    {
      sizes[0] = 2;
      sizes[1] = s1;
      sizes[2] = s2;
      sizes[3] = 0;
    }
  public:
    ~Arr2()
    {
      if(owner) AlignedFree(arr);
    }
    Arr2(int s1, int s2) :
      owner(true)
    {
      set_sizes(s1,s2);
      arr = AlignedAlloc(type, s1*s2);
    }
    Arr2(type*const* in, int s1, int s2) :
      arr(*in),
      owner(false)
    {
      set_sizes(s1,s2);
    }
    inline Ref1<type> operator[](int idx){
      #ifdef CHECK_BOUNDS
      assert_le(0, idx); assert_lt(idx, sizes[1]);
      #endif
      return Ref1<type>(arr, sizes+2, idx);
    }
    int getidx(int n1, int n2) const
    {
        #ifdef CHECK_BOUNDS
        assert_le(0, n1); assert_lt(n1, sizes[1]);
        assert_le(0, n2); assert_lt(n2, sizes[2]);
        #endif
        int k = n1;
        k *= sizes[2]; k += n2;
        return k;
    }
    const type& get(int n1,int n2) const
      { ALIGNED(arr); return arr[getidx(n1,n2)]; }
    type& fetch(int n1,int n2)
      { ALIGNED(arr); return arr[getidx(n1,n2)]; }
    void set(int n1,int n2, type value)
      { ALIGNED(arr); arr[getidx(n1,n2)] = value; }
};

template <class type>
class Arr3
{
  private: // data
    type* __restrict__ arr;
    int sizes[5];
    bool owner;
  private: // methods
    void set_sizes(int s1, int s2, int s3)
    {
      sizes[0] = 3;
      sizes[1] = s1;
      sizes[2] = s2;
      sizes[3] = s3;
      sizes[4] = 0;
    }
  public:
    ~Arr3()
    {
      if(owner) AlignedFree(arr);
    }
    Arr3(int s1, int s2, int s3) :
      owner(true)
    {
      set_sizes(s1,s2,s3);
      arr = AlignedAlloc(type, s1*s2*s3);
    }
    Arr3(type*const*const* in, int s1, int s2, int s3) :
      arr(**in),
      owner(false)
    {
      set_sizes(s1,s2,s3);
    }
    inline Ref2<type> operator[](int idx){
      #ifdef CHECK_BOUNDS
      assert_le(0, idx); assert_lt(idx, sizes[1]);
      #endif
      return Ref2<type>(arr, sizes+2, idx);
    }
    int getidx(int n1, int n2, int n3) const
    {
        #ifdef CHECK_BOUNDS
        assert_le(0, n1); assert_lt(n1, sizes[1]);
        assert_le(0, n2); assert_lt(n2, sizes[2]);
        assert_le(0, n3); assert_lt(n3, sizes[3]);
        #endif
        int k = n1;
        k *= sizes[2]; k += n2;
        k *= sizes[3]; k += n3;
        return k;
    }
    const type& get(int n1,int n2,int n3) const
      { ALIGNED(arr); return arr[getidx(n1,n2,n3)]; }
    type& fetch(int n1,int n2,int n3)
      { ALIGNED(arr); return arr[getidx(n1,n2,n3)]; }
    void set(int n1,int n2,int n3, type value)
      { ALIGNED(arr); arr[getidx(n1,n2,n3)] = value; }
};

template <class type>
class Arr4
{
  private: // data
    type* __restrict__ arr;
    int sizes[6];
    bool owner;
  private: // methods
    void set_sizes(int s1, int s2, int s3, int s4)
    {
      sizes[0] = 4;
      sizes[1] = s1;
      sizes[2] = s2;
      sizes[3] = s3;
      sizes[4] = s4;
      sizes[5] = 0;
    }
  public:
    ~Arr4()
    {
      if(owner) AlignedFree(arr);
    }
    Arr4(int s1, int s2, int s3, int s4) :
      owner(true)
    {
      set_sizes(s1,s2,s3,s4);
      arr = AlignedAlloc(type, s1*s2*s3*s4);
    }
    Arr4(type*const*const*const* in,
      int s1, int s2, int s3, int s4) :
      arr(***in),
      owner(false)
    {
      set_sizes(s1,s2,s3,s4);
    }
    inline Ref3<type> operator[](int idx){
      #ifdef CHECK_BOUNDS
      assert_le(0, idx); assert_lt(idx, sizes[1]);
      #endif
      return Ref3<type>(arr, sizes+2, idx);
    }
    int getidx(int n1, int n2, int n3, int n4) const
    {
        #ifdef CHECK_BOUNDS
        assert_le(0, n1); assert_lt(n1, sizes[1]);
        assert_le(0, n2); assert_lt(n2, sizes[2]);
        assert_le(0, n3); assert_lt(n3, sizes[3]);
        assert_le(0, n4); assert_lt(n4, sizes[4]);
        #endif
        int k = n1;
        k *= sizes[2]; k += n2;
        k *= sizes[3]; k += n3;
        k *= sizes[4]; k += n4;
        return k;
    }
    const type& get(int n1,int n2,int n3,int n4) const
      { return arr[getidx(n1,n2,n3,n4)]; }
    type& fetch(int n1,int n2,int n3,int n4)
      { return arr[getidx(n1,n2,n3,n4)]; }
    void set(int n1,int n2,int n3,int n4, type value)
      { arr[getidx(n1,n2,n3,n4)] = value; }
};

// aliases to avoid filling the code with template brackets
//
typedef Arr1<int> intArr1;
typedef Arr2<int> intArr2;
typedef Arr3<int> intArr3;
typedef Arr4<int> intArr4;
typedef Ref1<int> intRef1;
typedef Ref2<int> intRef2;
typedef Ref3<int> intRef3;
typedef Arr1<double> doubleArr1;
typedef Arr2<double> doubleArr2;
typedef Arr3<double> doubleArr3;
typedef Arr4<double> doubleArr4;
typedef Ref1<double> doubleRef1;
typedef Ref2<double> doubleRef2;
typedef Ref3<double> doubleRef3;

#define newArr4(type,sz1,sz2,sz3,sz4) newArray4<type>((sz1),(sz2),(sz3),(sz4))
#define newArr3(type,sz1,sz2,sz3) newArray3<type>((sz1),(sz2),(sz3))
#define newArr2(type,sz1,sz2) newArray2<type>((sz1),(sz2))
/****** end Alec Johnson's code ******/
#endif
