#ifndef ipicarrays_h
#define ipicarrays_h
#include "Alloc.h"

// should make a const version of this class
class vector_array3_double
{
  int numel;
  array3_double*list;
  ~vector_array3_double()
  {
    for(int i=0;i<numel;i++)
    {
      list[i].~array3_double();
    }
    free(list);
  }
  vector_array3_double(int numel_, int nx, int ny, int nz):
    numel(numel_),
    list(0)
  {
    list = malloc(numel*sizeof(array3_double*));
    for(int i=0;i<numel;i++)
    {
      new(&list[i]) array3_double(nx,ny,nz);
    }
  }
   array3_double& operator[](int i)
   {
     #ifdef CHECK_BOUNDS
     assert_le(0,i);
     assert_lt(i,numel);
     #endif // CHECK_BOUNDS
     return list[i];
   }
};
#endif // ipicarrays_h
