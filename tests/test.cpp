#include <stdio.h>
#include <malloc.h>
#define NO_NEW
#include "Alloc.h"

using namespace iPic3D;

int main()
{
  printf("hello world\n");
  int* arr = (int*)malloc(sizeof(int)*5);
  for(int i=0;i<5;i++)
    arr[i] = i;
  for(int i=0;i<5;i++)
    printf("arr[%d] = %d\n", i, arr[i]);

  array1<int> arr1(5);
  for(int i=0;i<5;i++)
    arr1[i] = i;
  for(int i=0;i<5;i++)
    printf("arr1[%d] = %d\n", i, arr1[i]);

  array2<int> arr2(2,3);
  for(int i=0;i<2;i++)
  for(int j=0;j<3;j++)
    arr2[i][j]=i+j;
  for(int i=0;i<2;i++)
  for(int j=0;j<3;j++)
    printf("arr2[%d][%d] = %d\n", i,j, arr2[i][j]);

  for(int i=0;i<2*3;i++)
    printf("arr2.fetch(%d) = %d\n", i, arr2.fetch(i));
    
  array3<int> arr3(2,3,4);
  //for(int i=0;i<2;i++)
  //for(int j=0;j<3;j++)
  //for(int k=0;k<4;k++)
  //  arr2[i][j]=i+j+k;
  //for(int i=0;i<2;i++)
  //for(int j=0;j<3;j++)
  //for(int k=0;k<4;k++)
  //  printf("arr3[%d][%d][%d] = %d\n", i,j,k, arr3[i][j][k]);
  //for(int i=0;i<2*3*4;i++)
  //  printf("arr3.fetch(%d) = %d\n", i, arr3.fetch(i));
}
