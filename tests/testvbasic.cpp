
// I started with Reza Rahman's
// "Intel Xeon Phi Coprocessor Architecture and Tools:
//   The Guide for Application Developers".
// I referred to:
// * shuffle intrinsics documentation:
//     http://software.intel.com/en-us/node/460872
// * vector mask intrinsics:
//     http://goo.gl/SLFmT7
// * intel xeon phi instruction set architecture reference manual:
//    http://software.intel.com/sites/default/files/forum/278102/327364001en.pdf
// * for double precision the mask bits are the bottom 8 bits:
//    http://stackoverflow.com/questions/18742863/manipulating-masks-for-doubles-on-xeon-phi
// * intrinsics for type casting:
//    http://software.intel.com/en-us/node/485381
// * masked swizzle parameters:
//    http://goo.gl/IXNDdP
// * SIMD MIC classes, e.g. F64vec8, M512 (__m512):
//    http://goo.gl/vMKWo5
// * SIMD Xeon classes (seem to be better documented):
//    http://vhdplus.sourceforge.net/doc/class_f64vec4.html

#define MICVEC_DEFINE_OUTPUT_OPERATORS
#include <iostream>
#include <micvec.h>
//#include <immintrin.h> // for __mm512_int2mask(0xf);
#include <stdio.h>
#include <stdint.h> // for int64_t

#define printexpr(var) cout << "line " << __LINE__ << ": " \
  << #var << " = " << var << endl;

#define print_double_arr(var) \
  { \
  printf("line %d: %s = (", __LINE__, #var); \
  for(int i=0; i<7; i++) printf("%g, ", var[i]); \
  printf("%g)\n", var[7]); \
  }

#define ASSUME_ALIGNED(X) __assume_aligned(X, 64)

using namespace std;

void testI8()
{
  _MM_PERM_ENUM p32;
  __declspec(align(64)) I64vec8 inputData(7,6,5,4,3,2,1,0);
  __declspec(align(64)) I64vec8 outputData;

  int64_t* idata = (int64_t*) &inputData;
  for(int i=0;i<8;i++) printf("idata[%d]=%d\n",i,idata[i]);

  std::cout << "input = " << inputData << endl;
  printexpr(inputData);

  // swizzle input data and print
  //std::cout << "swizzle data for pattern 'cdab' \n" << inputData.cdab() << endl;
  printexpr(inputData.cdab()); // use for 2x2 matrix transpose
  printexpr(inputData.badc()); // use for 4x4 matrix transpose of 2x2 blocks
  printexpr(inputData.dacb()); // use for cross product
  printexpr(inputData.aaaa());
  printexpr(inputData.bbbb());
  printexpr(inputData.cccc());
  printexpr(inputData.dddd());
}

void testF8()
{
  _MM_PERM_ENUM p32;
  __declspec(align(64)) F64vec8 data1(1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.8);
  __declspec(align(64)) F64vec8 data2(2.7,2.6,2.5,2.4,2.3,2.2,2.1,2.8);
  __declspec(align(64)) F64vec8 outp1;
  __declspec(align(64)) F64vec8 outp2;
  __declspec(align(64)) F64vec8 outputData;

  double* idata = (double*) &data1;
  for(int i=0;i<8;i++) printf("idata[%d]=%g\n",i,idata[i]);

  std::cout << "input = " << data1 << endl;

  // swizzle input data and print
  //std::cout << "swizzle data for pattern 'cdab' \n" << inputData.cdab() << endl;
  printexpr(data1.cdab());
  printexpr(data1.badc()); // use for 2x2 transpose
  printexpr(data1.dacb());
  printexpr(data1.aaaa());
  printexpr(data1.bbbb());
  printexpr(data1.cccc());
  printexpr(data1.dddd());

  // mask values
  //
  // rmask  hex   mask  dec   oct
  // 0000 : 0x0 = 0000 =  0 =   0
  // 1000 : 0x1 = 0001 =  1 =   1
  // 0100 : 0x2 = 0010 =  2 =   2
  // 1100 : 0x3 = 0011 =  3 =   3
  // 0010 : 0x4 = 0100 =  4 =   4
  // 1010 : 0x5 = 0101 =  5 =   5
  // 0110 : 0x6 = 0110 =  6 =   6
  // 1110 : 0x7 = 0111 =  7 =   7
  // 0001 : 0x8 = 1000 =  8 =  01
  // 1001 : 0x9 = 1001 =  9 =  02
  // 0101 : 0xa = 1010 = 10 =  03
  // 1101 : 0xb = 1011 = 11 =  04
  // 0011 : 0xc = 1100 = 12 =  05
  // 1011 : 0xd = 1101 = 13 =  06
  // 0111 : 0xe = 1110 = 14 =  07
  // 1111 : 0xf = 1111 = 15 = 010
  
  // rmask_11001100 (reversed mask) would equal mask_00110011
  // rmask bits go from lowest to highest
  // I designate the bits from lowest to highest,
  // but intel goes from highest to lowest
}

// transpose each 2x2 block in the two rows of a 2x8 matrix
// used to transpose the 2x2 blocks of the 8x8 matrix
inline void trans2x2(double in1[8], double in2[8])
{
  F64vec8& data1 = (F64vec8&)in1[0];
  F64vec8& data2 = (F64vec8&)in2[0];
  // copy the data that we will first overwrite
  const F64vec8 buff1 = data1;

  // data1: copy odds from evens of data2,
  //   i.e.,
  // for(int i=1; i<8; i+=2) data1[i] = data2[i-1];
  //
  const __mmask16 rmask_01010101 = _mm512_int2mask(0xaa); // mask=10101010
  // replace unmasked data1 with swizzled data2
  // (hopefully compiler will see that this can be done in
  // a single vector instruction)
  data1 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data1),
    rmask_01010101, __m512d(data2),_MM_SWIZ_REG_CDAB));
  //cout << "swizzle for pattern 'cdab' with rmask=01010101\n"
  //  << data1 << endl;

  // data2: copy evens from odds of buff1,
  //   i.e.,
  // for(int i=0; i<7; i+=2) data2[i] = buff1[i+1];
  //
  const __mmask16 rmask_10101010 = _mm512_int2mask(0x55); // mask=01010101
  // replace unmasked data2 with swizzled buff1
  data2 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data2),
    rmask_10101010, __m512d(buff1),_MM_SWIZ_REG_CDAB));
  //cout << "swizzle for pattern 'cdab' with rmask=10101010\n"
  //  << data2 << endl;
}

// transpose 1x2 elements of each 2x4 block in the two rows of a 2x8 matrix
// used to transpose the 2x2 elements of the 4x4 blocks of the 8x8 matrix
void trans4x4(double in1[8], double in2[8])
{
  F64vec8& data1 = (F64vec8&)in1[0];
  F64vec8& data2 = (F64vec8&)in2[0];
  // copy the data that we will first overwrite
  const F64vec8 buff1 = data1;

  // data1: copy odd 1x2 elements from even 1x2 elements of data2,
  //   i.e.,
  // for(int i=2; i<8; i++) if((i/2)%2) data1[i] = data2[i-2];
  //
  // replace unmasked data1 with swizzled data2
  const __mmask16 rmask_00110011 = _mm512_int2mask(0xcc); // mask = 11001100
  data1 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data1),
    rmask_00110011, __m512d(data2),_MM_SWIZ_REG_BADC));
  //cout << "swizzle for pattern 'badc' with rmask=00110011\n"
  //  << outp1 << endl;

  // data2: copy even 1x2 elements from odd 1x2 from odds of buff1,
  //   i.e.,
  // for(int i=0; i<6; i++) if(!((i/2)%2)) data2[i] = buff1[i+2];
  //
  const __mmask16 rmask_11001100 = _mm512_int2mask(0x33); // mask=00110011
  data2 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data2),
    rmask_11001100, __m512d(buff1),_MM_SWIZ_REG_BADC));
  //cout << "swizzle for pattern 'badc' with rmask=11001100\n"
  //  << outp2 << endl;
}

// transpose 1x4 elements of each 2x4 block in the two rows of a 2x8 matrix
// used to transpose the 4x4 elements of the 8x8 matrix
void trans8x8(double in1[8], double in2[8])
{
  F64vec8& data1 = (F64vec8&)in1[0];
  F64vec8& data2 = (F64vec8&)in2[0];
  // copy the data that we will first overwrite
  const F64vec8 buff1 = data1;

  // for swizzle intel supports 256-bit lanes with 64-bit elements
  // (as well as 128-bit lanes with 32-bit elements),
  // but for shuffle intel only supports 128-bit lanes with 32-bit
  // elements, so we must cast to 32-bit data types, do the
  // shuffle, and cast back.
  //
  // The amount of (nested) casting needed here
  // in order to "convert" between 32-bit integer and
  // 64-bit double precision data seems kind of unbelievable...

  // data1: copy low (odd) 1x4 element from even (high) 1x4 element of data2,
  //   i.e.,
  // for(int i=0; i<4; i++) data1[i+4] = data2[i];
  //
  // replace unmasked data1 with shuffled data2
  const __mmask16 rmask_00001111 = _mm512_int2mask(0xff00); //mask=11110000
  data1 = F64vec8(_mm512_castps_pd(
    _mm512_mask_permute4f128_ps(_mm512_castpd_ps(__m512d(data1)),
      rmask_00001111, _mm512_castpd_ps(__m512d(data2)),_MM_PERM_BADC)));
  //cout << "inter lane shuffle for pattern 'badc' with rmask=11110000\n"
  //  << data << endl;

  // data2: copy high (even) 1x4 element from odd (low) 1x4 element of buff1,
  //   i.e.,
  // for(int i=0; i<4; i++) data2[i] = buff1[i+4];
  //
  // replace unmasked data2 with shuffled data1
  const __mmask16 rmask_11110000 = _mm512_int2mask(0x00ff); //mask=00001111
  data2 = F64vec8(_mm512_castps_pd(
    _mm512_mask_permute4f128_ps(_mm512_castpd_ps(__m512d(data2)),
      rmask_11110000, _mm512_castpd_ps(__m512d(buff1)),_MM_PERM_BADC)));
  //cout << "inter lane shuffle for pattern 'badc' with rmask=00001111\n"
  //  << data2 << endl;
}
void print_data(double data[8][8], int line)
{
  printf("=== at line %d: data is: ===\n", line);
  for(int i=0; i<8; i++)
  {
    printf("  ");
    for(int j=0; j<7; j++)
    {
      printf("%g, ", data[i][j]);
    }
    printf("%g)\n", data[i][7]);
  }
  printf("=== end of data ===\n");
}

// transpose with blocked algorithm in
// 24=8*3 = 8*log_2(8) vector instructions
// (ignoring 12=4*3 copies made to
// buffers to reduce number of registers needed)
//
void transpose_8x8_double(double data[8][8])
{
  ASSUME_ALIGNED(data);
  // 1. transpose each 2x2 block.
  for(int i=0; i<8; i+=2)
    trans2x2(data[i], data[i+1]);
  //// 2. transpose each 4x4 block of 2x2 elements
  trans4x4(data[0], data[2]);
  trans4x4(data[1], data[3]);
  trans4x4(data[4], data[6]);
  trans4x4(data[5], data[7]);
  //print_data(data,__LINE__);
  //// 3. swap lower left and upper right 4x4 elements
  for(int i=0; i<4; i+=1)
    trans8x8(data[i], data[i+4]);
}

void test_transpose_8x8_double()
{
  __declspec(align(64)) double data[8][8];
  // initialize with indexing data
  for(int i=0;i<8;i++)
  for(int j=0;j<8;j++)
    data[i][j]=(1+i)+.1*(j+1);

  print_data(data,__LINE__);
  transpose_8x8_double(data);
  print_data(data,__LINE__);
}

// copy applying mask
inline void masked_copy(__m512d& dst, __m512d src, __mmask8 mask)
{
    dst = _mm512_mask_mov_pd(dst,mask,src);
}
// copy low vector of src to low vector of dst
inline void copy012to012(F64vec8& dst, F64vec8 src)
{
  const __mmask8 rmask_11100000 = _mm512_int2mask(0x07); // mask=00000111
  // masked copy
  dst = _mm512_mask_mov_pd(dst,rmask_11100000,src);
}
// copy hgh vector of src to hgh vector of dst
inline void copy456to456(F64vec8& dst, F64vec8 src)
{
  const __mmask8 rmask_00001110 = _mm512_int2mask(0x70); // mask=01110000
  dst = _mm512_mask_mov_pd(dst,rmask_00001110,src);
}
// copy low vector of src to hgh vector of dst
inline void copy012to456(F64vec8& dst, F64vec8 src)
{
  // mask=0000000000111111 because 32-bit elements
  const __mmask16 rmask_1111110000000000 = _mm512_int2mask(0x03f);
  // write the first three elements of src beginning at dst[4]
  _mm512_mask_packstorelo_epi32(&dst[4],rmask_1111110000000000,*(__m512i*)&src);
}
// copy hgh vector of src to low vector of dst
inline void copy456to012(F64vec8& dst, F64vec8 src)
{
  // mask=0011111111111111 because 32-bit elements
  const __mmask16 rmask_1111111111111100 = _mm512_int2mask(0x3fff);
  // write all but the highest element to the left of dst[4].
  _mm512_mask_packstorehi_epi32(&dst[4],rmask_1111111111111100,*(__m512i*)&src);
}

int test_copy_methods()
{
  // const causes compiler error
  /*const*/ F64vec8 u(8,7,6,5,4,3,2,1);
  F32vec16 v = _mm512_cvtpd_pslo(u); 
  F32vec16 vinv = F32vec16(1.)/v;
  printexpr(vinv);
  F64vec8 vinvd = _mm512_cvtpslo_pd(vinv);
  printexpr(vinvd);
  printexpr(F64vec8(1.)/u);
  // store the part of the stream before the first 64-bit
  // byte-aligned address preceding 
  F64vec8 arr = F64vec8(0.);
  //double arr[8]ALLOC_ALIGNED;
  //for(int i=0;i<8;i++) arr[i] = 0.;
  print_double_arr(u);
  print_double_arr(arr);
  copy012to012(arr, u);
  print_double_arr(arr);
  copy456to456(arr, u);
  print_double_arr(arr);
  copy456to012(arr, u);
  print_double_arr(arr);
  copy012to456(arr, u);
  print_double_arr(arr);
  //_mm512_packstorehi_epi32(&arr[4],*(__m512i*)&u);
  //print_double_arr(arr);
  // _mm512_packstorelo_epi32(&arr[5],*(__m512i*)&u);

  const F64vec8 u1(8.1,7.1,6.1,5.1,4.1,3.1,2.1,1.1);
  const F64vec8 u2(8.2,7.2,6.2,5.2,4.2,3.2,2.2,1.2);
  const __mmask8 rmask_11101110 = _mm512_int2mask(0x77); //mask=01110111
  const F64vec8 out = _mm512_mask_mov_pd(u1,rmask_11101110,u2);
  print_double_arr(out);
}

// calculate the cross product of the corresponding
// physical vectors in the two halves of the arguments:
//
// for(int i=0;i<8;i+=4)
// {
//   w[i+0] = u[i+1]*v[i+2] - u[i+2]*v[i+1];
//   w[i+1] = u[i+2]*v[i+0] - u[i+0]*v[i+2];
//   w[i+2] = u[i+0]*v[i+1] - u[i+1]*v[i+0];
// }
inline F64vec8 cross_product(F64vec8 u, F64vec8 v)
{
  F64vec8 us = u.dacb();
  F64vec8 vs = v.dacb();
  return us*vs.dacb() - us.dacb()*vs;
}

void test_cross_product()
{
  // calculate u cross v
  F64vec8 u(1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1);
  F64vec8 v(2.8,2.7,2.6,2.5,2.4,2.3,2.2,2.1);
  F64vec8 w = cross_product(u,v);
  printexpr(w);
  for(int i=0;i<8;i+=4)
  {
    w[i+0] = u[i+1]*v[i+2] - u[i+2]*v[i+1];
    w[i+1] = u[i+2]*v[i+0] - u[i+0]*v[i+2];
    w[i+2] = u[i+0]*v[i+1] - u[i+1]*v[i+0];
  }
  printexpr(w);
}

int main()
{
  testI8();
  testF8();
  test_transpose_8x8_double();
  test_copy_methods();
  test_cross_product();
}
