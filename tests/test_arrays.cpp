/*
   Reger Ferrer
   Vicenç Beltran
   Alec Johnson
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "stopwatch.h"
#include "Alloc.h"
#include "asserts.h"
#include "debug.h"

/*! The allocator for 4D array */
template < class type > type **** newArray4_Amaya(int sz1, int sz2, int sz3, int sz4) {

  type ****all_x;
  type ***all_y;
  type **all_z;
  type *all_r;

  all_x = new type ***[sz1];
  all_y = new type **[sz1 * sz2];
  all_z = new type *[sz1 * sz2 * sz3];
  all_r = new type[sz1 * sz2 * sz3 * sz4];

  type ****result = all_x;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
    for (int j = 0; j < sz2; j++, all_z += sz3) {
      result[i][j] = all_z;
      for (int k = 0; k < sz3; k++, all_r += sz4) {
        result[i][j][k] = all_r;
      }
    }
  }

  return result;
}

/*! Deallocator for 4D arrays */
template < class type > void delArr4_Amaya(type **** arr, int dummyx, int dummyy, int dummyz) {
  delete[]arr[0][0][0];
  delete[]arr[0][0];
  delete[]arr[0];
  delete[]arr;
}

/*! The allocator for 3D array */
template < class type > type *** newArray3_Amaya(int sz1, int sz2, int sz3) {

  type ***all_x;
  type **all_y;
  type *all_z;

  all_x = new type **[sz1];
  all_y = new type *[sz1 * sz2];
  all_z = new type[sz1 * sz2 * sz3];

  type ***result = all_x;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
    for (int j = 0; j < sz2; j++, all_z += sz3) {
      result[i][j] = all_z;
    }
  }

  return result;

}

/*! Deallocator for 3D arrays */
template < class type > void delArr3_Amaya(type *** arr, int dummyx, int dummyy) {
  delete[]arr[0][0];
  delete[]arr[0];
  delete[]arr;
}

/*! The allocator for 2D array */
template < class type > type ** newArr2_Amaya(int sz1, int sz2) {

  type **all_x;
  type *all_y;

  all_x = new type *[sz1];
  all_y = new type[sz1 * sz2];

  type **result = all_x;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
  }

  return result;

}

/*! Deallocator for 2D arrays */
template < class type > void delArr2_Amaya(type ** arr, int dummyx) {
  delete[]arr[0];
  delete[]arr;
}

#define newArr4_Amaya(type,sz1,sz2,sz3,sz4) newArray4_Amaya<type>((sz1),(sz2),(sz3),(sz4))
#define newArr3_Amaya(type,sz1,sz2,sz3) newArray3_Amaya<type>((sz1),(sz2),(sz3))
#define newArr2_Amaya(type,sz1,sz2) newArray2_Amaya<type>((sz1),(sz2))

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

/*** Array classes if dimensions are known at compile time. ***/

template <class type, size_t N1, size_t N2>
class Array2D
{
public:
   type data [N1][N2];
};

template <class type, size_t N1, size_t N2, size_t N3>
class Array3D
{
public:
   type data [N1][N2][N3];
};

/******** end [i][j] arrays from Reger Ferrer and Vicenç Beltran ******/

using namespace std;

template <class type>
void testArr2_diagonal()
{
   const int ITERS = 10000;
   const size_t dim1 = 64;
   const size_t dim2 = 64;

   Array2<type> Abra(dim1, dim2);
   Array2<type> Bbra(dim1, dim2);
   Array2<type> Cbra(dim1, dim2);

   Rank2<type> Apar(dim1, dim2);
   Rank2<type> Bpar(dim1, dim2);
   Rank2<type> Cpar(dim1, dim2);

   Array2D<type, dim1, dim2> Afix ;
   Array2D<type, dim1, dim2> Bfix ;
   Array2D<type, dim1, dim2> Cfix ;

   type** Aold = newArr2(type, dim1, dim2);
   type** Bold = newArr2(type, dim1, dim2);
   type** Cold = newArr2(type, dim1, dim2);

   Arr2<type> Aeaj(Aold, dim1, dim2);
   Arr2<type> Beaj(Bold, dim1, dim2);
   Arr2<type> Ceaj(Cold, dim1, dim2);

   printf("Initializing data ...\n");
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Bbra[i][j] = rand();
      Cbra[i][j] = rand();
      Bpar(i,j) = Bbra[i][j];
      Cpar(i,j) = Cbra[i][j];
      Bfix.data[i][j] = Bbra[i][j];
      Cfix.data[i][j] = Cbra[i][j];
      Bold[i][j] = Bbra[i][j];
      Cold[i][j] = Cbra[i][j];
      Beaj.fetch(i,j) = Bbra[i][j];
      Ceaj.fetch(i,j) = Cbra[i][j];
   }

   stopwatch(START);
   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Abra[i][j] = Bbra[i][j] * Cbra[i][j];
   }
   printf("%d ms = Total time [i][j] Vincenc array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Apar(i,j) = Bpar(i,j) * Cpar(i,j);
   }
   printf("%d ms = Total time (i,j) Vincenc array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Afix.data[i][j] = Bfix.data[i][j] * Cfix.data[i][j];
   }
   printf("%d ms = Total time [i][j] fixed-dimension array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Aold[i][j] = Bold[i][j] * Cold[i][j];
   }
   printf("%d ms = Total time [i][j] chained-pointer array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Aeaj[i][j] = Beaj[i][j] * Ceaj[i][j];
   }
   printf("%d ms = Total time [i][j] access of Arr2\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Aeaj.fetch(i,j) = Beaj.get(i,j) * Ceaj.get(i,j);
   }
   printf("%d ms = Total time (i,j) access of Arr2\n", tv_to_ms(stopwatch(LAP)));

   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      assert(Afix.data[i][j] == Abra[i][j]);
      assert(Aold[i][j] == Abra[i][j]);
      assert(Aeaj.get(i,j) == Abra[i][j]);
   }

   printf("Verification done!\n");
   stopwatch(STOP);

   delArr2(Aold,dim1);
   delArr2(Bold,dim1);
   delArr2(Cold,dim1);
}

template <class type>
void testArr2()
{
   const int ITERS = 10000;
   const size_t dim1 = 64;
   const size_t dim2 = 64;

   Array2<type> Abra(dim1, dim2);
   Array2<type> Bbra(dim1, dim2);
   Array2<type> Cbra(dim1, dim2);

   Rank2<type> Apar(dim1, dim2);
   Rank2<type> Bpar(dim1, dim2);
   Rank2<type> Cpar(dim1, dim2);

   Array2D<type, dim1, dim2> Afix ;
   Array2D<type, dim1, dim2> Bfix ;
   Array2D<type, dim1, dim2> Cfix ;

   type** Aold = newArr2(type, dim1, dim2);
   type** Bold = newArr2(type, dim1, dim2);
   type** Cold = newArr2(type, dim1, dim2);

   Arr2<type> Aeaj(Aold, dim1, dim2);
   Arr2<type> Beaj(Bold, dim1, dim2);
   Arr2<type> Ceaj(Cold, dim1, dim2);

   printf("Initializing data ...\n");
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Bbra[i][j] = rand();
      Cbra[i][j] = rand();
      Bpar(i,j) = Bbra[i][j];
      Cpar(i,j) = Cbra[i][j];
      Bfix.data[i][j] = Bbra[i][j];
      Cfix.data[i][j] = Cbra[i][j];
      Bold[i][j] = Bbra[i][j];
      Cold[i][j] = Cbra[i][j];
      Beaj.fetch(i,j) = Bbra[i][j];
      Ceaj.fetch(i,j) = Cbra[i][j];
   }

   stopwatch(START);
   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Abra[i][j] = Bbra[i][j] * Cbra[i][j];
   }
   printf("%d ms = Total time [i][j] Vincenc array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Apar(i,j) = Bpar(i,j) * Cpar(i,j);
   }
   printf("%d ms = Total time (i,j) Vincenc array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Afix.data[i][j] = Bfix.data[i][j] * Cfix.data[i][j];
   }
   printf("%d ms = Total time [i][j] fixed-dimension array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Aold[i][j] = Bold[i][j] * Cold[i][j];
   }
   printf("%d ms = Total time [i][j] chained-pointer array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Aeaj[i][j] = Beaj[i][j] * Ceaj[i][j];
   }
   printf("%d ms = Total time [i][j] access of Arr2\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Aeaj.fetch(i,j) = Beaj.get(i,j) * Ceaj.get(i,j);
   }
   printf("%d ms = Total time (i,j) access of Arr2\n", tv_to_ms(stopwatch(LAP)));

   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      assert(Afix.data[i][j] == Abra[i][j]);
      assert(Aold[i][j] == Abra[i][j]);
      assert(Aeaj.get(i,j) == Abra[i][j]);
   }

   printf("Verification done!\n");
   stopwatch(STOP);

   delArr2(Aold,dim1);
   delArr2(Bold,dim1);
   delArr2(Cold,dim1);
}

template <class type>
void testArr3()
{
   const int ITERS = 100;
   const size_t dim1 = 64;
   const size_t dim2 = 64;
   const size_t dim3 = 64;

   // Using this class causes a segmentation fault (why?)
   //Array3<type> Abra(dim1, dim2, dim3);
   //Array3<type> Bbra(dim1, dim2, dim3);
   //Array3<type> Cbra(dim1, dim2, dim3);

   Rank3<type> Apar(dim1, dim2, dim3);
   Rank3<type> Bpar(dim1, dim2, dim3);
   Rank3<type> Cpar(dim1, dim2, dim3);

   Array3D<type, dim1, dim2, dim3> Afix ;
   Array3D<type, dim1, dim2, dim3> Bfix ;
   Array3D<type, dim1, dim2, dim3> Cfix ;

   type*** Aold = newArr3(type, dim1, dim2, dim3);
   type*** Bold = newArr3(type, dim1, dim2, dim3);
   type*** Cold = newArr3(type, dim1, dim2, dim3);

   //Arr3<type> Aeaj(Aold, dim1, dim2, dim3);
   //Arr3<type> Beaj(Bold, dim1, dim2, dim3);
   //Arr3<type> Ceaj(Cold, dim1, dim2, dim3);
   Arr3<type> Aeaj(dim1, dim2, dim3);
   Arr3<type> Beaj(dim1, dim2, dim3);
   Arr3<type> Ceaj(dim1, dim2, dim3);

   printf("Initializing data ...\n");
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Beaj.fetch(i,j,k) = rand();
      Ceaj.fetch(i,j,k) = rand();
      //Bbra[i][j][k] = Beaj.get(i,j,k);
      //Cbra[i][j][k] = Ceaj.get(i,j,k);
      Bpar(i,j,k) = Beaj.get(i,j,k);
      Cpar(i,j,k) = Ceaj.get(i,j,k);
      Bfix.data[i][j][k] = Beaj.get(i,j,k);
      Cfix.data[i][j][k] = Ceaj.get(i,j,k);
      Bold[i][j][k] = Beaj.get(i,j,k);
      Cold[i][j][k] = Ceaj.get(i,j,k);
   }

   stopwatch(START);

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Apar(i,j,k) = Bpar(i,j,k) * Cpar(i,j,k);
   }
   printf("%d ms = Total time (i,j,k) Vincenc array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Afix.data[i][j][k] = Bfix.data[i][j][k] * Cfix.data[i][j][k];
   }
   printf("%d ms = Total time [i][j][k] fixed-dimension array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Aold[i][j][k] = Bold[i][j][k] * Cold[i][j][k];
   }
   printf("%d ms = Total time [i][j][k] chained-pointer array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Aeaj[i][j][k] = Beaj[i][j][k] * Ceaj[i][j][k];
   }
   printf("%d ms = Total time [i][j][k] access of Arr3\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Aeaj.fetch(i,j,k) = Beaj.get(i,j,k) * Ceaj.get(i,j,k);
   }
   printf("%d ms = Total time (i,j,k) access of Arr3\n", tv_to_ms(stopwatch(LAP)));

   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      assert_eq(Aold[i][j][k], Aeaj.get(i,j,k));
      assert_eq(Apar(i,j,k), Aeaj.get(i,j,k));
      assert_eq(Afix.data[i][j][k], Aeaj.get(i,j,k));
   }

   printf("Verification done!\n");
   stopwatch(STOP);
}

int main()
{
  printf("=== testing Arr2<int> (diagonal) ===\n");
  testArr2_diagonal<int>();
  printf("=== testing Arr2<double> (diagonal) ===\n");
  testArr2_diagonal<double>();
  printf("=== testing Arr2<int> ===\n");
  testArr2<int>();
  printf("=== testing Arr2<double> ===\n");
  testArr2<double>();
  printf("=== testing Arr3<int> ===\n");
  testArr3<int>();
  printf("=== testing Arr3<double> ===\n");
  testArr3<double>();
}
