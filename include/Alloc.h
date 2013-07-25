#ifndef IPIC_ALLOC_H
#define IPIC_ALLOC_H
#include <cstddef>

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
    #define AlignedAlloc(T, NUM) (T *const __restrict__)(_mm_malloc(sizeof(T)*NUM, ALIGNMENT))
    #define AlignedFree(S) (_mm_free(S))
#else
    #define ALIGNED(X)
    #define AlignedFree(S) (delete[] S)
    #define AlignedAlloc(T, NUM) (new T[NUM]) 
#endif

/*** Array classes if dimensions are known at compile time. ***/

template <class type, size_t N, size_t M>
class Array2D
{
public:
   type data [N][M];
};

template <class type, size_t N1, size_t N2, size_t N3>
class Array3D
{
public:
   type data [N1][N2][N3];
};

/******** end [i][j] arrays from Reger Ferrer and Vicenç Beltran ******/

/****** begin Alec Johnson's code ******/

// methods to allocate arrays
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
  type ***arr = AlignedAlloc(type**, sz1); // new type **[sz1];
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
// The purpose of this class is to allow elements of multidimensional
// arrays to be accessed with a calculated one-dimensional index while
// using the same syntax as is used for a nested array.  This gives
// correct results, but unfortunately is too slow; evidently the
// compiler is not intelligent enough to figure out that the whole chain
// of calls to the operator[] methods and to the DerefN constructors
// reduces to computing a one-dimensional subscript used to access a
// one-dimensional array.
//
template <class type>
class Deref1
{
  type* __restrict__ const arr;
  const size_t* const sizes;
  size_t k;
 public:
  inline Deref1(type *const arr_in, size_t const* const sizes_in, size_t k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline type& operator[](size_t idx){
    k *= sizes[0]; k += idx;
    ALIGNED(arr);
    return arr[k];
  }
};

template <class type>
class Deref2
{
  type* __restrict__ const arr;
  const size_t* const sizes;
  size_t k;
  inline Deref2(type *const arr_in, size_t const*const sizes_in, size_t k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline Deref1<type> operator[](size_t idx){
    k *= sizes[0]; k += idx;
    return Deref1<type>(arr, sizes+1, k);
  }
};

template <class type>
class Deref3
{
  type* __restrict__ const arr;
  const size_t* const sizes;
  size_t k;
  inline Deref3(type*const arr_in, size_t const*const sizes_in, size_t k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline Deref2<type> operator[](size_t idx){
    k *= sizes[0]; k += idx;
    return Deref2<type>(arr, sizes+1, k);
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

// class that can adopt array allocated by newArr1
template <class type>
class Arr1
{
  private:
    type* __restrict__ arr;
    size_t sizes[3];
    bool owner;
  private:
    void set_size(size_t s1)
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
    Arr1(size_t s1) :
      owner(true)
    {
      set_size(s1);
      arr = AlignedAlloc(type, s1);
    }
    inline Arr1(type* in, size_t s1) :
      arr(in),
      owner(false)
    {
      set_size(s1);
    }
    inline type operator[](size_t idx){
      ALIGNED(arr);
      return arr[idx];
    }
};

// class that can adopt array allocated by newArr2
template <class type>
class Arr2
{
  private:
    type* __restrict__ arr;
    size_t sizes[4];
    bool owner;
  private:
    void set_sizes(size_t s1, size_t s2)
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
    Arr2(size_t s1, size_t s2) :
      owner(true)
    {
      set_sizes(s1,s2);
      arr = AlignedAlloc(type, s1*s2);
    }
    Arr2(type*const* in, size_t s1, size_t s2) :
      arr(*in),
      owner(false)
    {
      set_sizes(s1,s2);
    }
    inline Deref1<type> operator[](size_t idx){
      return Deref1<type>(arr, sizes+2, idx);
    }
    size_t getidx(size_t n1, size_t n2) const
    {
        size_t k = n1;
        k *= sizes[2]; k += n2;
        return k;
    }
    const type& get(size_t n1,size_t n2) const
      { ALIGNED(arr); return arr[getidx(n1,n2)]; }
    type& fetch(size_t n1,size_t n2)
      { ALIGNED(arr); return arr[getidx(n1,n2)]; }
    void set(size_t n1,size_t n2, type value)
      { ALIGNED(arr); arr[getidx(n1,n2)] = value; }
};

// class that can adopt array allocated by newArr3
template <class type>
class Arr3
{
  private:
    type* __restrict__ arr;
    size_t sizes[5];
    bool owner;
  private:
    void set_sizes(size_t s1, size_t s2, size_t s3)
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
    Arr3(size_t s1, size_t s2, size_t s3) :
      owner(true)
    {
      set_sizes(s1,s2,s3);
      arr = AlignedAlloc(type, s1*s2*s3);
    }
    Arr3(type*const*const* in, size_t s1, size_t s2, size_t s3) :
      arr(**in),
      owner(false)
    {
      set_sizes(s1,s2,s3);
    }
    inline Deref2<type> operator[](size_t idx){
      return Deref2<type>(arr, sizes+2, idx);
    }
    size_t getidx(size_t n1, size_t n2, size_t n3) const
    {
        size_t k = n1;
        k *= sizes[2]; k += n2;
        k *= sizes[3]; k += n3;
        return k;
    }
    const type& get(size_t n1,size_t n2,size_t n3) const
      { ALIGNED(arr); return arr[getidx(n1,n2,n3)]; }
    type& fetch(size_t n1,size_t n2,size_t n3)
      { ALIGNED(arr); return arr[getidx(n1,n2,n3)]; }
    void set(size_t n1,size_t n2,size_t n3, type value)
      { ALIGNED(arr); arr[getidx(n1,n2,n3)] = value; }
};

// class that can adopt array allocated by newArr4
template <class type>
class Arr4
{
  private: // data
    type* __restrict__ arr;
    size_t sizes[6];
    bool owner;
  private: // methods
    void set_sizes(size_t s1, size_t s2, size_t s3, size_t s4)
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
    Arr4(size_t s1, size_t s2, size_t s3, size_t s4) :
      owner(true)
    {
      set_sizes(s1,s2,s3,s4);
      arr = AlignedAlloc(type, s1*s2*s3*s4);
    }
    Arr4(type*const*const*const* in,
      size_t s1, size_t s2, size_t s3, size_t s4) :
      arr(***in),
      owner(false)
    {
      set_sizes(s1,s2,s3,s4);
    }
    inline Deref3<type> operator[](size_t idx){
      return Deref3<type>(arr, sizes+2, idx);
    }
    size_t getidx(size_t n1, size_t n2, size_t n3, size_t n4) const
    {
        size_t k = n1;
        k *= sizes[2]; k += n2;
        k *= sizes[3]; k += n3;
        k *= sizes[4]; k += n4;
        return k;
    }
    const type& get(size_t n1,size_t n2,size_t n3,size_t n4) const
      { return arr[getidx(n1,n2,n3,n4)]; }
    type& fetch(size_t n1,size_t n2,size_t n3,size_t n4)
      { return arr[getidx(n1,n2,n3,n4)]; }
    void set(size_t n1,size_t n2,size_t n3,size_t n4, type value)
      { arr[getidx(n1,n2,n3,n4)] = value; }
};

#define newArr4(type,sz1,sz2,sz3,sz4) newArray4<type>((sz1),(sz2),(sz3),(sz4))
#define newArr3(type,sz1,sz2,sz3) newArray3<type>((sz1),(sz2),(sz3))
#define newArr2(type,sz1,sz2) newArray2<type>((sz1),(sz2))
/****** end Alec Johnson's code ******/
#endif
