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

#undef CHECK_BOUNDS

/****** begin Alec Johnson's code ******/

// methods to allocate arrays.
// These are a succinct equivalent of Jorge's earler methods,
// except for the use of AlignedAlloc in place of new.
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
// RefN is essentially a dumbed-down version of ArrN with an
// index shift applied to the underlying array.  The purpose of
// RefN is to allow elements of multidimensional arrays to be
// accessed with a calculated one-dimensional index while using
// the same syntax as is used for a nested array.  This gives
// correct results, but is slow unless optimization is turned on,
// allowing the compiler to figure out that the whole chain of
// calls to the operator[] methods and to the RefN constructors
// reduces to computing a one-dimensional subscript used to
// access a one-dimensional array.
//
template <class type>
class Ref1
{
  type* __restrict__ const arr;
  int shift;
  int S1;
 public:
  inline Ref1(type*const arr_, int k, int s1) :
    arr(arr_), shift(k), S1(s1)
  {}
  inline type& operator[](int n1){
    #ifdef CHECK_BOUNDS
    assert_le(0, idx); assert_lt(n1, S1);
    #endif
    return arr[shift+n1];
  }
};

template <class type>
class Ref2
{
  type* __restrict__ const arr;
  int shift;
  int S2, S1;
 public:
  inline Ref2(type*const arr_, int k, int s2, int s1) :
    arr(arr_), shift(k), S2(s2), S1(s1)
  {}
  inline Ref1<type> operator[](int n2){
    #ifdef CHECK_BOUNDS
    assert_le(0, idx); assert_lt(n2, S2);
    #endif
    return Ref1<type>(arr, (shift+n2)*S1, S1);
  }
};

template <class type>
class Ref3
{
  type* __restrict__ const arr;
  int shift;
  int S3, S2, S1;
 public:
  inline Ref3(type*const arr_, int k, int s3, int s2, int s1) :
    arr(arr_), shift(k), S3(s3), S2(s2), S1(s1)
  {}
  inline Ref2<type> operator[](int n3){
    #ifdef CHECK_BOUNDS
    assert_le(0, idx); assert_lt(n3, S3);
    #endif
    return Ref2<type>(arr, (shift+n3)*S2, S2, S1);
  }
};

// ArrN can adopt an array allocated by newArrN
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
    int S1;
    bool owner;
  public:
    ~Arr1()
    {
      if(owner) AlignedFree(arr);
    }
    Arr1(int s1) :
      S1(s1),
      owner(true)
    {
      arr = AlignedAlloc(type, s1);
    }
    Arr1(type* in,
      int s1) :
      S1(s1),
      arr(in),
      owner(false)
    { }
    inline type& operator[](int n1){
      #ifdef CHECK_BOUNDS
      assert_le(0, n1); assert_lt(n1, S1);
      #endif
      return arr[n1];
    }
    inline int getidx(int n1) const
    {
      #ifdef CHECK_BOUNDS
      assert_le(0, n1); assert_lt(n1, S1);
      #endif
      return n1;
    }
    const type& get(int n1) const
      { return arr[getidx(n1)]; }
    type& fetch(int n2,int n1)
      { return arr[getidx(n1)]; }
    void set(int n1, type value)
      { arr[getidx(n1)] = value; }
};

template <class type>
class Arr2
{
  private: // data
    type* __restrict__ arr;
    int S2,S1;
    bool owner;
  public:
    ~Arr2()
    {
      if(owner) AlignedFree(arr);
    }
    Arr2(int s2, int s1) :
      S2(s2), S1(s1),
      owner(true)
    {
      arr = AlignedAlloc(type, ss2*s1);
    }
    Arr2(type*const* in,
      int s2, int s1) :
      S2(s2), S1(s1),
      arr(*in),
      owner(false)
    { }
    inline Ref1<type> operator[](int n2){
      #ifdef CHECK_BOUNDS
      assert_le(0, n2); assert_lt(n2, S2);
      #endif
      return Ref1<type>(arr, n2*S1, S1);
    }
    inline int getidx(int n2, int n1) const
    {
      #ifdef CHECK_BOUNDS
      assert_le(0, n2); assert_lt(n2, S2);
      assert_le(0, n1); assert_lt(n1, S1);
      #endif
      return n2*S1+n1;
    }
    const type& get(int n2,int n1) const
      { return arr[getidx(n2,n1)]; }
    type& fetch(int n2,int n1)
      { return arr[getidx(n2,n1)]; }
    void set(int n2,int n1, type value)
      { arr[getidx(n2,n1)] = value; }
};

template <class type>
class Arr3
{
  private: // data
    type* __restrict__ arr;
    int S3,S2,S1;
    bool owner;
  public:
    ~Arr3()
    {
      if(owner) AlignedFree(arr);
    }
    Arr3(int s3, int s2, int s1) :
      S3(s3), S2(s2), S1(s1),
      owner(true)
    {
      arr = AlignedAlloc(type, s3*s2*s1);
    }
    Arr3(type*const*const* in,
      int s3, int s2, int s1) :
      S3(s3), S2(s2), S1(s1),
      arr(**in),
      owner(false)
    { }
    inline Ref2<type> operator[](int n3){
      #ifdef CHECK_BOUNDS
      assert_le(0, n3); assert_lt(n3, S3);
      #endif
      return Ref2<type>(arr, n3*S2, S2, S1);
    }
    inline int getidx(int n3, int n2, int n1) const
    {
      #ifdef CHECK_BOUNDS
      assert_le(0, n3); assert_lt(n3, S3);
      assert_le(0, n2); assert_lt(n2, S2);
      assert_le(0, n1); assert_lt(n1, S1);
      #endif
      return (n3*S2+n2)*S1+n1;
    }
    const type& get(int n3,int n2,int n1) const
      { return arr[getidx(n3,n2,n1)]; }
    type& fetch(int n3,int n2,int n1)
      { return arr[getidx(n3,n2,n1)]; }
    void set(int n3,int n2,int n1, type value)
      { arr[getidx(n3,n2,n1)] = value; }
};

template <class type>
class Arr4
{
  private: // data
    type* __restrict__ arr;
    int S4,S3,S2,S1;
    bool owner;
  public:
    ~Arr4()
    {
      if(owner) AlignedFree(arr);
    }
    Arr4(int s4, int s3, int s2, int s1) :
      S4(s4), S3(s3), S2(s2), S1(s1),
      owner(true)
    {
      arr = AlignedAlloc(type, s4*s3*s2*s1);
    }
    Arr4(type*const*const*const* in,
      int s4, int s3, int s2, int s1) :
      S4(s4), S3(s3), S2(s2), S1(s1),
      arr(***in),
      owner(false)
    { }
    inline Ref3<type> operator[](int n4){
      #ifdef CHECK_BOUNDS
      assert_le(0, n4); assert_lt(n4, S4);
      #endif
      return Ref3<type>(arr, n4*S3, S3, S2, S1);
      //return Arr3<type>(arr, n4*S3, S3, S2, S1); // could do this instead...
    }
    inline int getidx(int n4, int n3, int n2, int n1) const
    {
      #ifdef CHECK_BOUNDS
      assert_le(0, n4); assert_lt(n4, S4);
      assert_le(0, n3); assert_lt(n3, S3);
      assert_le(0, n2); assert_lt(n2, S2);
      assert_le(0, n1); assert_lt(n1, S1);
      #endif
      return ((n4*S3+n3)*S2+n2)*S1+n1;
    }
    const type& get(int n4,int n3,int n2,int n1) const
      { return arr[getidx(n4,n3,n2,n1)]; }
    type& fetch(int n4,int n3,int n2,int n1)
      { return arr[getidx(n4,n3,n2,n1)]; }
    void set(int n4,int n3,int n2,int n1, type value)
      { arr[getidx(n4,n3,n2,n1)] = value; }
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
