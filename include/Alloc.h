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
#ifdef CHECK_BOUNDS
  #define check_bounds(n,S) {assert_le(0, n); assert_lt(n, S);}
#else
  #define check_bounds(n,S)
#endif

/****** begin Alec Johnson's code ******/

// methods to allocate arrays.
// These are a succinct equivalent of Jorge's earler methods,
// except for the use of AlignedAlloc in place of new.
//
template < class type > inline type * newArray1(size_t sz1) {
  type *arr = AlignedAlloc(type, sz1); // new type [sz1];
  return arr;
}
template < class type > inline type ** newArray2(size_t sz1, size_t sz2) {
  type **arr = AlignedAlloc(type*, sz1); // new type *[sz1];
  type *ptr = newArray1<type>(sz1*sz2);
  for (size_t i = 0; i < sz1; i++)
  {
    arr[i] = ptr;
    ptr += sz2;
  }
  return arr;
}
template < class type > inline type *** newArray3(size_t sz1, size_t sz2, size_t sz3) {
  // type ***arr = AlignedAlloc(type**, sz1); // new type **[sz1];
  type ***arr = new type **[sz1];
  type **ptr = newArray2<type>(sz1*sz2, sz3);
  for (size_t i = 0; i < sz1; i++)
  {
    arr[i] = ptr;
    ptr += sz2;
  }
  return arr;
}
template < class type > inline type **** newArray4(size_t sz1, size_t sz2, size_t sz3, size_t sz4) {
  type ****arr = AlignedAlloc(type***, sz1); //(new type ***[sz1]);
  type ***ptr = newArray3<type>(sz1*sz2, sz3, sz4);
  for (size_t i = 0; i < sz1; i++) {
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
template < class type > inline void delArr2(type ** arr, size_t sz1)
{
  delArr1(arr[0]);
  AlignedFree(arr);
}
template < class type > inline void delArr3(type *** arr, size_t sz1, size_t sz2)
{
  delArr2(arr[0],sz2);
  AlignedFree(arr);
}
template < class type > inline void delArr4(type **** arr,
  size_t sz1, size_t sz2, size_t sz3)
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
  type* const __restrict__ arr;
  const size_t S1;
  const size_t shift;
 public:
  inline Ref1(type*const arr_, size_t k, size_t s1) :
    arr(arr_), shift(k), S1(s1)
  {}
  inline type& operator[](size_t n1){
    check_bounds(n1, S1);
    ALIGNED(arr);
    return arr[shift+n1];
  }
};

template <class type>
class Ref2
{
  type* const __restrict__ arr;
  const size_t shift;
  const size_t S2, S1;
 public:
  inline Ref2(type*const arr_, size_t k, size_t s2, size_t s1) :
    arr(arr_), shift(k), S2(s2), S1(s1)
  {}
  inline Ref1<type> operator[](size_t n2){
    check_bounds(n2,S2);
    return Ref1<type>(arr, (shift+n2)*S1, S1);
  }
};

template <class type>
class Ref3
{
  type* const __restrict__ arr;
  const size_t shift;
  const size_t S3, S2, S1;
 public:
  inline Ref3(type*const arr_, size_t k, size_t s3, size_t s2, size_t s1) :
    arr(arr_), shift(k), S3(s3), S2(s2), S1(s1)
  {}
  inline Ref2<type> operator[](size_t n3){
    check_bounds(n3, S3);
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
// - inheriting wrapper with automatic destructor
//   that returns this as a "fast accessor"
// - methods that use parallel arithmetic for omp and vectorized code

template <class type>
class Arr1
{
  private: // data
    type* const __restrict__ arr;
    const size_t S1;
  public:
    ~Arr1() { }
    void free() { AlignedFree(arr); }
    Arr1(size_t s1) :
      S1(s1),
      arr(AlignedAlloc(type, s1))
    { }
    Arr1(type* in,
      size_t s1) :
      S1(s1),
      arr(in)
    { }
    inline type& operator[](size_t n1){
      check_bounds(n1, S1);
      ALIGNED(arr);
      return arr[n1];
    }
    inline size_t getidx(size_t n1) const
    {
      check_bounds(n1, S1);
      return n1;
    }
    const type& get(size_t n1) const
      { ALIGNED(arr); return arr[getidx(n1)]; }
    type& fetch(size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[getidx(n1)]; }
    void set(size_t n1, type value)
      { ALIGNED(arr); arr[getidx(n1)] = value; }
};

template <class type>
class Arr2
{
  private: // data
    const size_t S2,S1;
    type* const __restrict__ arr;
  public:
    ~Arr2(){}
    void free() { AlignedFree(arr); }
    Arr2(size_t s2, size_t s1) :
      S2(s2), S1(s1),
      arr(AlignedAlloc(type, s2*s1))
    {
    }
    Arr2(type*const* in,
      size_t s2, size_t s1) :
      S2(s2), S1(s1),
      arr(*in)
    { }
    inline Ref1<type> operator[](size_t n2){
      check_bounds(n2, S2);
      return Ref1<type>(arr, n2*S1, S1);
    }
    inline size_t getidx(size_t n2, size_t n1) const
    {
      check_bounds(n2, S2);
      check_bounds(n1, S1);
      return n2*S1+n1;
    }
    type& operator()(size_t n2, size_t n1) const
      { ALIGNED(arr); return arr[n1+S1*n2]; }
    type& fetch(size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[n1+S1*n2]; } //arr[n2*S1+n1]; }
    //{ ALIGNED(arr); return arr[getidx(n2,n1)]; }
    const type& get(size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[n1+S1*n2]; } // arr[n2*S1+n1]; 
    //{ ALIGNED(arr); return arr[getidx(n2,n1)]; }
    void set(size_t n2,size_t n1, type value)
      { ALIGNED(arr); arr[[n2*S1+n1]] = value; }
    //{ ALIGNED(arr); arr[getidx(n2,n1)] = value; }
};

template <class type>
class Arr3
{
  private: // data
    type* const __restrict__ arr;
    const size_t S3,S2,S1;
  public:
    ~Arr3(){} // nonempty destructor would kill performance
    void free() { AlignedFree(arr); }
    Arr3(size_t s3, size_t s2, size_t s1) :
      S3(s3), S2(s2), S1(s1),
      arr(AlignedAlloc(type, s3*s2*s1))
    { }
    Arr3(type*const*const* in,
      size_t s3, size_t s2, size_t s1) :
      S3(s3), S2(s2), S1(s1),
      arr(**in)
    { }
    inline Ref2<type> operator[](size_t n3){
      check_bounds(n3, S3);
      return Ref2<type>(arr, n3*S2, S2, S1);
    }
    inline size_t getidx(size_t n3, size_t n2, size_t n1) const
    {
      check_bounds(n3, S3);
      check_bounds(n2, S2);
      check_bounds(n1, S1);
      return (n3*S2+n2)*S1+n1;
    }
    type& fetch(size_t n3,size_t n2,size_t n1) const
    {
       ALIGNED(arr);
       return arr[n1+S1*(n2+S2*n3)];
    }
    //{ ALIGNED(arr); return arr[(n3*S2+n2)*S1+n1]; }
    //{ ALIGNED(arr); return arr[getidx(n3,n2,n1)]; }
    const type& get(size_t n3,size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[(n3*S2+n2)*S1+n1]; }
    //{ ALIGNED(arr); return arr[getidx(n3,n2,n1)]; }
    void set(size_t n3,size_t n2,size_t n1, type value)
      { ALIGNED(arr); arr[(n3*S2+n2)*S1+n1] = value; }
    //{ ALIGNED(arr); arr[getidx(n3,n2,n1)] = value; }
    type& operator()(size_t n3, size_t n2, size_t n1) const
    {
       ALIGNED(arr);
       return arr[n1+S1*(n2+S2*n3)];
    }
};

template <class type>
class Arr4
{
  private: // data
    const size_t S4,S3,S2,S1;
    type* const __restrict__ arr;
  public:
    ~Arr4(){} // nonempty destructor would kill performance
    void free() { AlignedFree(arr); }
    Arr4(size_t s4, size_t s3, size_t s2, size_t s1) :
      arr(AlignedAlloc(type, s4*s3*s2*s1)),
      S4(s4), S3(s3), S2(s2), S1(s1)
    { }
    Arr4(type*const*const*const* in,
      size_t s4, size_t s3, size_t s2, size_t s1) :
      S4(s4), S3(s3), S2(s2), S1(s1),
      arr(***in)
    { }
    inline Ref3<type> operator[](size_t n4){
      check_bounds(n4, S4);
      return Ref3<type>(arr, n4*S3, S3, S2, S1);
    }
    inline size_t getidx(size_t n4, size_t n3, size_t n2, size_t n1) const
    {
      check_bounds(n4, S4);
      check_bounds(n3, S3);
      check_bounds(n2, S2);
      check_bounds(n1, S1);
      return ((n4*S3+n3)*S2+n2)*S1+n1;
    }
    const type& get(size_t n4,size_t n3,size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[getidx(n4,n3,n2,n1)]; }
    type& fetch(size_t n4,size_t n3,size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[getidx(n4,n3,n2,n1)]; }
    void set(size_t n4,size_t n3,size_t n2,size_t n1, type value)
      { ALIGNED(arr); arr[getidx(n4,n3,n2,n1)] = value; }
};

/**** begin assimilated array classes for debugging ***/

template <class type>
class ArrRank3
{
  private: // data
    const size_t S3,S2,S1;
    type* const __restrict__ arr;
  public:
    // nonempty destructor kills performance
    ~ArrRank3()
    { }
    ArrRank3(size_t s3, size_t s2, size_t s1) :
      S3(s3), S2(s2), S1(s1),
      arr(AlignedAlloc(type, s3*s2*s1))
    {
    }
    ArrRank3(type*const*const* in,
      size_t s3, size_t s2, size_t s1) :
      S3(s3), S2(s2), S1(s1)
      arr(**in)
    { }
    inline Ref2<type> operator[](size_t n3){
      check_bounds(n3, S3);
      return Ref2<type>(arr, n3*S2, S2, S1);
    }
    inline size_t getidx(size_t n3, size_t n2, size_t n1) const
    {
      check_bounds(n3, S3);
      check_bounds(n2, S2);
      check_bounds(n1, S1);
      return (n3*S2+n2)*S1+n1;
    }
    type& fetch(size_t n3,size_t n2,size_t n1) const
    { ALIGNED(arr); return arr[(n3*S2+n2)*S1+n1]; }
    //{ ALIGNED(arr); return arr[getidx(n3,n2,n1)]; }
    const type& get(size_t n3,size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[(n3*S2+n2)*S1+n1]; }
    //{ ALIGNED(arr); return arr[getidx(n3,n2,n1)]; }
    void set(size_t n3,size_t n2,size_t n1, type value)
      { ALIGNED(arr); arr[(n3*S2+n2)*S1+n1] = value; }
    //{ ALIGNED(arr); arr[getidx(n3,n2,n1)] = value; }
    type& operator()(size_t n3, size_t n2, size_t n1) const
    {
       ALIGNED(arr);
       return arr[n1+S1*(n2+S2*n3)];
    }
    void free()
    {
      AlignedFree(arr);
    };
};

/**** end assimilated array classes for debugging ***/

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
