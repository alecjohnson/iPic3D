/****** begin (i,j) arrays from Reger Ferrer and Vicenç Beltran ******/

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

template <class ElementType>
class Rank1
{
    const size_t _N1;
    ElementType  * __restrict__ const  _storage;

public:

    Rank1(size_t N1) : _N1(N1), _storage(AlignedAlloc(ElementType, N1)) {}

    //Rank1( const Rank1& other ) : _N1( other._N1 ), _storage( other._storage ) {}

    ElementType& operator()(size_t n1) const
    {
       ALIGNED(_storage);
       return _storage[n1];
    }

    size_t dim1() const { return _N1; }

    ~Rank1() { };
};

template <class ElementType>
class Rank2
{
    const size_t  _N1, _N2;
    ElementType * __restrict__ const _storage;


public:
    Rank2(size_t N1, size_t N2) : _N1(N1), _N2(N2), _storage(AlignedAlloc(ElementType, N1*N2)) {}

    //Rank2( const Rank2& other ) : _N1( other._N1 ), _N2( other._N2 ), _storage( other._storage ) {}

    ElementType& operator()(size_t n1, size_t n2) const
    {
       ALIGNED(_storage);
       return _storage[n2+_N2*n1];
    }

    size_t dim1() const { return _N1; }
    size_t dim2() const { return _N2; }
    
    void deleteArr() {
       AlignedFree(_storage);
    }

    ~Rank2() { };
};
    
template <class ElementType>
class Rank3
{
    const size_t _N1, _N2, _N3;
    ElementType *    const __restrict__ _storage;


public:

    Rank3(size_t N1, size_t N2, size_t N3) : _N1(N1), _N2(N2), _N3(N3),
    _storage(AlignedAlloc(ElementType, N1*N2*N3)) {}

    //Rank3( const Rank3& other ) : _N1( other._N1 ), _N2( other._N2 ), _N3( other._N3 ),
    //_storage( other._storage ) {}

    ElementType& operator()(size_t n1, size_t n2, size_t n3) const
    {
       ALIGNED(_storage);
       return _storage[n3+_N3*(n2+_N2*n1)];
    }

    ~Rank3() { }

    size_t dim1() const { return _N1; }
    size_t dim2() const { return _N2; }
    size_t dim3() const { return _N3; }

    void deleteArr() {
       AlignedFree(_storage);
    }
};

template <class ElementType>
class Rank4
{
    const size_t _N1, _N2, _N3, _N4;
    ElementType* __restrict__ const _storage;

public:

    Rank4(size_t N1, size_t N2, size_t N3, size_t N4) : _N1(N1), _N2(N2), _N3(N3), _N4(N4),
    _storage(AlignedAlloc(ElementType, N1*N2*N3*N4)) {}

    //Rank4( const Rank4& other ) : _N1( other._N1 ), _N2( other._N2 ), _N3( other._N3 ), _N4( other._N4 ),
    //_storage( other._storage ) {}

    ElementType& operator()(size_t n1, size_t n2, size_t n3, size_t n4) const
    {
       ALIGNED(_storage);
       return _storage[n4+_N4*(n3+_N3*(n2+_N2*n1))];
    }

    ~Rank4() { }

    size_t dim1() const { return _N1; } 
    size_t dim2() const { return _N2; }
    size_t dim3() const { return _N3; }
    size_t dim4() const { return _N4; }
    
    void deleteArr() { AlignedFree(_storage); }

};

/******** end (i,j) arrays from Reger Ferrer and Vicenç Beltran ******/

/****** begin [i][j] arrays from Reger Ferrer and Vicenç Beltran ******/

template <class ElementType>
class Array1
{
    const size_t _N1;
    ElementType  * __restrict__ const  _storage;

public:
    Array1(size_t N1, void * __restrict__ const storage) : _N1(N1),
         _storage(reinterpret_cast<ElementType * __restrict__ const>(storage)){}

    Array1(size_t N1) : _N1(N1), _storage(new ElementType[N1]){}

    ElementType& operator[](size_t i) const
    {
        return _storage[i];
    }

};

template <class ElementType>
class Array2
{
    const size_t _N1, _N2;
    ElementType * __restrict__ const _storage;

public:
    Array2(size_t N1, size_t N2, void *storage) : _N1(N1), _N2(N2),
         _storage(reinterpret_cast<ElementType * __restrict__ const>(storage)){}

    Array2(size_t N1, size_t N2) : _N1(N1), _N2(N2),
         _storage(new ElementType[N1*N2]) {}

    Array1<ElementType> operator[](size_t i) const
    {
        return Array1<ElementType>(_N2, _storage + i * _N2);
    }
    ElementType& operator()(size_t n1, size_t n2) const
    {
       ALIGNED(_storage);
       return _storage[n2+_N2*n1];
    }
};

template <class ElementType>
class Array3
{
    const size_t _N1, _N2, _N3;
    ElementType* __restrict__ const _storage;

public:
    Array3(size_t N1, size_t N2, size_t N3, void *storage) : _N1(N1), _N2(N2), _N3(N3),
         _storage(reinterpret_cast<ElementType * __restrict__ const>(storage)) {}


    Array3(size_t N1, size_t N2, size_t N3) : _N1(N1), _N2(N2), _N3(N3),
         _storage(new ElementType(N1*N2*N3)) {}

    Array2<ElementType> operator[](size_t i) const
    {
        return Array2<ElementType>(_N2, _N3, _storage + i * _N2 * _N3);
    }
    ElementType& operator()(size_t n1, size_t n2, size_t n3) const
    {
       ALIGNED(_storage);
       return _storage[n3+_N3*(n2+_N2*n1)];
    }
};

template <class ElementType>
class Array4
{
    const size_t _N1, _N2, _N3, _N4;
    ElementType* __restrict__ const _storage;

public:
    Array4(size_t N1, size_t N2, size_t N3, size_t N4, void *storage) : _N1(N1), _N2(N2), _N3(N3), _N4(N4),
         _storage(reinterpret_cast<ElementType * __restrict__ const>(storage)) {}

    Array4(size_t N1, size_t N2, size_t N3, size_t N4) : _N1(N1), _N2(N2), _N3(N3), _N4(N4),
         _storage(new ElementType(N1*N2*N3*N4)) {}

    Array3<ElementType> operator[](size_t i) const
    {
        return Array3<ElementType>(_N2, _N3, _N4, _storage + i * _N2 * _N3 * _N4);
    }
    ElementType& operator()(size_t n1, size_t n2, size_t n3, size_t n4) const
    {
       ALIGNED(_storage);
       return _storage[n4+_N4*(n3+_N3*(n2+_N2*n1))];
    }

};

/*** Array classes if dimensions are known at compile time. ***/

template <size_t N, size_t M>
class IntArray2D
{
public:
   int data [N][M];
};

template <size_t N1, size_t N2, size_t N3>
class IntArray3D
{
public:
   int data [N1][N2][N3];
};

template <size_t N, size_t M>
class DoubleArray2D
{
public:
   double data [N][M];
};

template <size_t N1, size_t N2, size_t N3>
class DoubleArray3D
{
public:
   double data [N1][N2][N3];
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
template < class type > inline void delArr1(type * arr) {
  AlignedFree(arr);
}
template < class type > inline void delArr2(type ** arr, size_t sz1) {
  delArr1(arr[0]);
  AlignedFree(arr);
}
template < class type > inline void delArr3(type *** arr, size_t sz1, size_t sz2) {
  delArr2(arr[0],sz2);
  AlignedFree(arr);
}
template < class type > inline void delArr4(type **** arr, size_t sz1, size_t sz2, size_t sz3) {
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
struct Deref1
{
  double* __restrict__ const arr;
  const size_t* const sizes;
  size_t k;
  inline Deref1(double *const arr_in, size_t const* const sizes_in, size_t k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline double& operator[](size_t idx){
    k *= sizes[0]; k += idx;
    ALIGNED(arr);
    return arr[k];
  }
};

struct Deref2
{
  double* __restrict__ const arr;
  const size_t* const sizes;
  size_t k;
  inline Deref2(double *const arr_in, size_t const*const sizes_in, size_t k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline Deref1 operator[](size_t idx){
    k *= sizes[0]; k += idx;
    return Deref1(arr, sizes+1, k);
  }
};

struct Deref3
{
  double* __restrict__ const arr;
  const size_t* const sizes;
  size_t k;
  inline Deref3(double*const arr_in, size_t const*const sizes_in, size_t k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline Deref2 operator[](size_t idx){
    k *= sizes[0]; k += idx;
    return Deref2(arr, sizes+1, k);
  }
};

// doubleArrN adopts an array allocated by newArrN
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
class doubleArr1
{
  private:
    double* __restrict__ const arr;
    size_t sizes[3];
  public:
    inline doubleArr1(double* in, size_t s1) :
      arr(in)
    {
      sizes[0] = 2; // not used
      sizes[1] = s1; //get_size(in);
      sizes[2] = 0;
    }
    inline double operator[](size_t idx){
      ALIGNED(arr);
      return arr[idx];
    }
};

// class that can adopt array allocated by newArr2
class doubleArr2
{
  private:
    double* __restrict__ const arr;
    size_t sizes[4];
  public:
    inline doubleArr2(double*const* in, size_t s1, size_t s2) :
      arr(*in)
    {
      sizes[0] = 2; // not used
      sizes[1] = s1; //get_size(in);
      sizes[2] = s2; //get_size(*in);
      sizes[3] = 0;
    }
    inline Deref1 operator[](size_t idx){
      return Deref1(arr, sizes+2, idx);
    }
    size_t getidx(size_t n1, size_t n2) const
    {
        size_t k = n1;
        k *= sizes[2]; k += n2;
        return k;
    }
    const double& get(size_t n1,size_t n2) const
      { ALIGNED(arr); return arr[getidx(n1,n2)]; }
    double& fetch(size_t n1,size_t n2)
      { ALIGNED(arr); return arr[getidx(n1,n2)]; }
    void set(size_t n1,size_t n2, double value)
      { ALIGNED(arr); arr[getidx(n1,n2)] = value; }
};

// class that can adopt array allocated by newArr3
class doubleArr3
{
  private:
    double* arr;
    size_t sizes[5];
  public:
    inline doubleArr3(double*const*const* in, size_t s1, size_t s2, size_t s3) :
      arr(**in)
    {
      sizes[0] = 3; // not used
      sizes[1] = s1; // get_size(in);
      sizes[2] = s2; // get_size(*in);
      sizes[3] = s3; // get_size(**in);
      sizes[4] = 0;
    }
    inline Deref2 operator[](size_t idx){
      return Deref2(arr, sizes+2, idx);
    }
    size_t getidx(size_t n1, size_t n2, size_t n3) const
    {
        size_t k = n1;
        k *= sizes[2]; k += n2;
        k *= sizes[3]; k += n3;
        return k;
    }
    const double& get(size_t n1,size_t n2,size_t n3) const
      { ALIGNED(arr); return arr[getidx(n1,n2,n3)]; }
    double& fetch(size_t n1,size_t n2,size_t n3)
      { ALIGNED(arr); return arr[getidx(n1,n2,n3)]; }
    void set(size_t n1,size_t n2,size_t n3, double value)
      { ALIGNED(arr); arr[getidx(n1,n2,n3)] = value; }
};

// class that can adopt array allocated by newArr4
class doubleArr4
{
  private:
    double* arr;
    size_t sizes[6];
  public:
    doubleArr4(double*const*const*const* in,
      size_t s1, size_t s2, size_t s3, size_t s4) :
      arr(***in)
    {
      sizes[0] = 4; // not used
      sizes[1] = s1; //get_size(in);
      sizes[2] = s2; //get_size(*in);
      sizes[3] = s3; //get_size(**in);
      sizes[4] = s4; //get_size(***in);
      sizes[5] = 0;
    }
    inline Deref3 operator[](size_t idx){
      return Deref3(arr, sizes+2, idx);
    }
    size_t getidx(size_t n1, size_t n2, size_t n3, size_t n4) const
    {
        size_t k = n1;
        k *= sizes[2]; k += n2;
        k *= sizes[3]; k += n3;
        k *= sizes[4]; k += n4;
        return k;
    }
    const double& get(size_t n1,size_t n2,size_t n3,size_t n4) const
      { return arr[getidx(n1,n2,n3,n4)]; }
    double& fetch(size_t n1,size_t n2,size_t n3,size_t n4)
      { return arr[getidx(n1,n2,n3,n4)]; }
    void set(size_t n1,size_t n2,size_t n3,size_t n4, double value)
      { arr[getidx(n1,n2,n3,n4)] = value; }
};

#define newArr4(type,sz1,sz2,sz3,sz4) newArray4<type>((sz1),(sz2),(sz3),(sz4))
#define newArr3(type,sz1,sz2,sz3) newArray3<type>((sz1),(sz2),(sz3))
#define newArr2(type,sz1,sz2) newArray2<type>((sz1),(sz2))
/****** end Alec Johnson's code ******/
