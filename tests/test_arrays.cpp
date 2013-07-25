/*
   Reger Ferrer
   Vicenç Beltran
   Alec Johnson
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "stopwatch.h"
#include "../include/Alloc.h"
#include "assert.h"

/****** begin (i,j) arrays from Reger Ferrer and Vicenç Beltran ******/

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

/******** end [i][j] arrays from Reger Ferrer and Vicenç Beltran ******/

using namespace std;

const int ITERS = 10000;
const size_t DIM_X = 64;
const size_t DIM_Y = 64;

template <class type>
void test_arrays()
{
   const size_t n = DIM_X;
   const size_t m = DIM_Y;

   Array2<type> Abra(n, m);
   Array2<type> Bbra(n, m);
   Array2<type> Cbra(n, m);

   Rank2<type> Apar(n, m);
   Rank2<type> Bpar(n, m);
   Rank2<type> Cpar(n, m);

   Array2D<type, DIM_X, DIM_Y> Afix ;
   Array2D<type, DIM_X, DIM_Y> Bfix ;
   Array2D<type, DIM_X, DIM_Y> Cfix ;

   type** Aold = newArr2(type, DIM_X, DIM_Y);
   type** Bold = newArr2(type, DIM_X, DIM_Y);
   type** Cold = newArr2(type, DIM_X, DIM_Y);

   Arr2<type> Aeaj(Aold, DIM_X, DIM_Y);
   Arr2<type> Beaj(Bold, DIM_X, DIM_Y);
   Arr2<type> Ceaj(Cold, DIM_X, DIM_Y);

   printf("Initializing data ...\n");
   for(size_t i=0; i<DIM_X; i++){
      for(size_t j=i; j<DIM_Y; j++){
         Bbra[i][j] = rand();
         Cbra[i][j] = rand();
         Bpar(i,j) = Bbra[i][j];
         Cpar(i,j) = Cbra[i][j];
         Bfix.data[i][j] = Bbra[i][j];
         Cfix.data[i][j] = Cbra[i][j];
         Bold[i][j] = Bbra[i][j];
         Cold[i][j] = Cbra[i][j];
      }
   }

   stopwatch(START);
   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<DIM_X; i++){
      for(size_t j=i; j<DIM_Y; j++){
         Abra[i][j] = Bbra[i][j] * Cbra[i][j];
      }
   }
   printf("Total time bracket array: %d ms\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<DIM_X; i++){
      for(size_t j=i; j<DIM_Y; j++){
         Apar(i,j) = Bpar(i,j) * Cpar(i,j);
      }
   }
   printf("Total time parentheses array: %d ms\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<DIM_X; i++){
      for(size_t j=i; j<DIM_Y; j++){
         Afix.data[i][j] = Bfix.data[i][j] * Cfix.data[i][j];
      }
   }
   printf("Total time fixed array: %d ms\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<DIM_X; i++){
      for(size_t j=i; j<DIM_Y; j++){
         Aold[i][j] = Bold[i][j] * Cold[i][j];
      }
   }
   printf("Total time indirect array: %d ms\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<DIM_X; i++){
      for(size_t j=i; j<DIM_Y; j++){
         Aeaj[i][j] = Beaj[i][j] * Ceaj[i][j];
      }
   }
   printf("Total time Alec [][] array: %d ms\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<DIM_X; i++){
      for(size_t j=i; j<DIM_Y; j++){
         Aeaj.fetch(i,j) = Beaj.get(i,j) * Ceaj.get(i,j);
      }
   }
   printf("Total time Alec (,) array: %d ms\n", tv_to_ms(stopwatch(LAP)));

   for(size_t i=0; i<DIM_X; i++){
      for(size_t j=i; j<DIM_Y; j++){
         assert(Aold[i][j] == Abra[i][j]);
         assert(Afix.data[i][j] == Abra[i][j]);
      }
   }

   printf("Verification done!\n");
   stopwatch(STOP);
}

int main()
{
  printf("=== testing integer arrays ===\n");
  test_arrays<int>();
  printf("=== testing double arrays ===\n");
  test_arrays<double>();
}
