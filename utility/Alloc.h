/***************************************************************************

  -------------------

developers           : Stefano Markidis, Giovanni Lapenta

This file was rewritten by eaj (Alec Johnson)

 ***************************************************************************/
#ifndef Alloc_H
#define Alloc_H
#include "stddef.h"
#include "assert.h"

// macros to allocate arrays
//
#define newArr1(type,sz1) _new_1d_array((sz1),(type *) NULL)
#define newArr2(type,sz1,sz2) _new_2d_array((sz1),(sz2),(type *) NULL)
#define newArr3(type,sz1,sz2,sz3) _new_3d_array((sz1),(sz2),(sz3),(type *) NULL)
#define newArr4(type,sz1,sz2,sz3,sz4) \
    _new_4d_array((sz1),(sz2),(sz3),(sz4),(type *) NULL);

// methods to allocate arrays
//
template < class type > inline type * _new_1d_array(int sz1, type * dummy) {
  type *arr = new type [sz1];
  //set_sizes(arr,sz1);
  return arr;
}
template < class type > inline type ** _new_2d_array(int sz1, int sz2, type * dummy) {
  type **arr = new type *[sz1];
  type *ptr = newArr1(type, sz1*sz2);
  for (int i = 0; i < sz1; i++)
  {
    arr[i] = ptr;
    ptr += sz2;
  }
  //set_sizes(arr,sz1);
  //set_size(*arr,sz2);
  return arr;
}
template < class type > inline type *** _new_3d_array(int sz1, int sz2, int sz3, type * dummy) {
  type ***arr = new type **[sz1];
  type **ptr = newArr2(type, sz1*sz2, sz3);
  for (int i = 0; i < sz1; i++)
  {
    arr[i] = ptr;
    ptr += sz2;
  }
  //set_sizes(arr,sz1);
  //set_size(*arr,sz2);
  return arr;
}
template < class type > inline type **** _new_4d_array(int sz1, int sz2, int sz3, int sz4, type * dummy) {
  type ****arr = (new type ***[sz1]);
  type ***ptr = newArr3(type, sz1*sz2, sz3, sz4);
  for (int i = 0; i < sz1; i++) {
    arr[i] = ptr;
    ptr += sz2;
  }
  //set_sizes(arr,sz1);
  //set_size(*arr,sz2);
  return arr;
}

// methods to deallocate arrays
//
template < class type > inline void delArr1(type * arr) {
  delete[](arr);
}
template < class type > inline void delArr2(type ** arr, int sz1) {
  delArr1(arr[0]);
  delete[](arr);
}
template < class type > inline void delArr3(type *** arr, int sz1, int sz2) {
  delArr2(arr[0],sz2);
  delete[](arr);
}
template < class type > inline void delArr4(type **** arr, int sz1, int sz2, int sz3) {
  delArr3(arr[0],sz2,sz3);
  delete[](arr);
}


class doubleArr3
{
  private:
    bool owned;
    double* arr1d;
    double*** arr;
    int sizes[5];
  private:
    void set_sizes(int s1, int s2, int s3)
    {
      sizes[0] = 3; // not used
      sizes[1] = s1;
      sizes[2] = s2;
      sizes[3] = s3;
      sizes[4] = 0;
    }
  public:
    doubleArr3(int s1, int s2, int s3)
    {
      owned = true;
      arr = newArr3(double, s1,s2,s3);
      set_sizes(s1,s2,s3);
    }
    doubleArr3(double*** in, int s1, int s2, int s3)
    {
      owned = false;
      arr = in;
      set_sizes(s1,s2,s3);
    }
    ~doubleArr3()
    {
      if(!owned) return;
      delArr3(arr, sizes[1],sizes[2]);
    }
    operator double***(){return arr;}
    int getidx(int n1, int n2, int n3) const
    {
        int k = n1;
        k *= sizes[2]; k += n2;
        k *= sizes[3]; k += n3;
        return k;
    }
    const double& get(int n1,int n2,int n3) const
      { return arr1d[getidx(n1,n2,n3)]; }
    double& fetch(int n1,int n2,int n3)
      { return arr1d[getidx(n1,n2,n3)]; }
    void set(int n1,int n2,int n3, double value)
      { arr1d[getidx(n1,n2,n3)] = value; }
};

class doubleArr4
{
  private:
    bool owned;
    double* arr1d;
    double**** arr;
    int sizes[6];
  private:
    void set_sizes(int s1, int s2, int s3, int s4)
    {
      sizes[0] = 4; // not used
      sizes[1] = s1;
      sizes[2] = s2;
      sizes[3] = s3;
      sizes[4] = s4;
      sizes[5] = 0;
    }
  public:
    doubleArr4(int s1, int s2, int s3, int s4)
    {
      owned = true;
      arr = newArr4(double, s1,s2,s3,s4);
      set_sizes(s1,s2,s3,s4);
    }
    doubleArr4(double**** in, int s1, int s2, int s3, int s4)
    {
      owned = false;
      arr = in;
      set_sizes(s1,s2,s3,s4);
    }
    ~doubleArr4()
    {
      if(!owned) return;
      delArr4(arr, sizes[1],sizes[2],sizes[3]);
    }
    operator double****(){return arr;}
    int getidx(int n1, int n2, int n3, int n4) const
    {
        int k = n1;
        k *= sizes[2]; k += n2;
        k *= sizes[3]; k += n3;
        k *= sizes[4]; k += n4;
        return k;
    }
    const double& get(int n1,int n2,int n3,int n4) const
      { return arr1d[getidx(n1,n2,n3,n4)]; }
    double& fetch(int n1,int n2,int n3,int n4)
      { return arr1d[getidx(n1,n2,n3,n4)]; }
    void set(int n1,int n2,int n3,int n4, double value)
      { arr1d[getidx(n1,n2,n3,n4)] = value; }
};

#endif
