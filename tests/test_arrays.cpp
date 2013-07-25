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

using namespace std;

const int ITERS = 10000;
const size_t DIM_X = 64;
const size_t DIM_Y = 64;

void test_int_arrays()
{
   const size_t n = DIM_X;
   const size_t m = DIM_Y;

   Array2<int> Abra(n, m);
   Array2<int> Bbra(n, m);
   Array2<int> Cbra(n, m);

   Rank2<int> Apar(n, m);
   Rank2<int> Bpar(n, m);
   Rank2<int> Cpar(n, m);

   IntArray2D<DIM_X, DIM_Y> Afix ;
   IntArray2D<DIM_X, DIM_Y> Bfix ;
   IntArray2D<DIM_X, DIM_Y> Cfix ;

   int** Aold = newArr2(int, DIM_X, DIM_Y);
   int** Bold = newArr2(int, DIM_X, DIM_Y);
   int** Cold = newArr2(int, DIM_X, DIM_Y);

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

   for(size_t i=0; i<DIM_X; i++){
      for(size_t j=i; j<DIM_Y; j++){
         assert(Aold[i][j] == Abra[i][j]);
         assert(Afix.data[i][j] == Abra[i][j]);
      }
   }

   printf("Verification done!\n");
   stopwatch(STOP);

}

void test_double_arrays()
{
   printf("=============================\n");
   const size_t n = DIM_X;
   const size_t m = DIM_Y;

   Array2<double> Abra(n, m);
   Array2<double> Bbra(n, m);
   Array2<double> Cbra(n, m);

   Rank2<double> Apar(n, m);
   Rank2<double> Bpar(n, m);
   Rank2<double> Cpar(n, m);

   DoubleArray2D<DIM_X, DIM_Y> Afix ;
   DoubleArray2D<DIM_X, DIM_Y> Bfix ;
   DoubleArray2D<DIM_X, DIM_Y> Cfix ;

   double** Aold = newArr2(double, DIM_X, DIM_Y);
   double** Bold = newArr2(double, DIM_X, DIM_Y);
   double** Cold = newArr2(double, DIM_X, DIM_Y);

   doubleArr2 Aeaj(Aold, DIM_X, DIM_Y);
   doubleArr2 Beaj(Bold, DIM_X, DIM_Y);
   doubleArr2 Ceaj(Cold, DIM_X, DIM_Y);

   printf("Initializing double data ...\n");
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

}

int main()
{
  test_int_arrays();
  test_double_arrays();
}
