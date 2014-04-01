#ifndef mic_basics_h
#define mic_basics_h

// documentation on intrinsics:
//
// * intel preprocessor macros:
//     http://goo.gl/CMS4WY
// * SIMD C++ classes:
//     http://software.intel.com/en-us/node/462716
// * AVX-512 intrinsics (for KNL):
//     http://software.intel.com/en-us/node/485150

#if defined(__MIC__)
  #include <micvec.h>
  #define printexpr(var) std::cout << "line " << __LINE__ << ": " \
    << #var << " = " << var << std::endl;
  #define MICVEC_DEFINE_OUTPUT_OPERATORS
  #include <iosfwd>
  // for some bizarre reason this, along with arithmetic and subscripting
  // operators, is not defined for this class, so I do it here.
  std::ostream& operator<<(std::ostream& os, const I32vec16 obj)
  {
    int* ref = (int*)&obj;
    os << "(";
    for(int i=0;i<15;i++)
      os << ref[i] << ",";
    os << ref[15] << ")";
    return os;
  }

  I32vec16 fused_multiply_add(I32vec16 v1, I32vec16 v2, I32vec16 v3)
  {
    return _mm512_fmadd_epi32(v1,v2,v3);
  }
  I32vec16 multiply(I32vec16 v1, I32vec16 v2)
  {
    return _mm512_mullo_epi32(v1,v2);
  }
  // unfortunately the built-in I32vec16 is missing
  // a lot of the functionality in e.g. F64vec8,
  // and unfortunately C++ does not provide a good
  // and simple mechanism to extend the interface of a class.
  class MyI32vec16: public I32vec16
  {
   public:
    //using I32vec16::operator+;
    int& operator[](int i)
    {
      int* memory = (int*)this;
      return memory[i];
    }
    // unfortunately in this case constructors aren't inherited
    MyI32vec16(int x15,int x14,int x13,int x12,int x11,int x10,int x9,int x8,
               int x7, int x6, int x5, int x4, int x3, int x2, int x1,int x0):
      I32vec16(x15,x14,x13,x12,x11,x10,x9,x8,x7,x6,x5,x4,x3,x2,x1,x0)
    {}
    MyI32vec16(): I32vec16()
    {}
    MyI32vec16(__m512i in): I32vec16(in)
    {}
    MyI32vec16(int in): I32vec16(in)
    {}
    MyI32vec16 operator*(MyI32vec16 in)
    {
      return _mm512_mulhi_epi32(*this,in);
    }
    MyI32vec16 operator+(MyI32vec16 in)
    {
      return _mm512_add_epi32(*this, in);
    }
    MyI32vec16 operator-(MyI32vec16 in)
    {
      return _mm512_sub_epi32(*this, in);
    }
  };

  inline F64vec8 maximum(F64vec8 u, F64vec8 v)
  {
    // which one?
    eprintf("implement me");
    return _mm512_max_pd(u,v);
    return _mm512_gmax_pd(u,v);
  }
  inline F64vec8 minimum(F64vec8 u, F64vec8 v)
  {
    eprintf("implement me");
    // which one?
    return _mm512_min_pd(u,v);
    return _mm512_gmin_pd(u,v);
  }

  inline F64vec8 cross_product(F64vec8 u, F64vec8 v)
  {
    F64vec8 us = u.dacb();
    F64vec8 vs = v.dacb();
    return us*vs.dacb() - us.dacb()*vs;
  }

  inline I32vec16 round_down(F64vec8 vec)
  { return _mm512_cvtfxpnt_roundpd_epi32lo(vec, _MM_ROUND_MODE_DOWN); }
  inline I32vec16 round_to_nearest(F64vec8 vec)
  { return _mm512_cvtfxpnt_roundpd_epi32lo(vec, _MM_ROUND_MODE_NEAREST); }

  inline F64vec8 make_F64vec8(double x, double y, double z, double t)
  { return F64vec8(t,z,y,x,t,z,y,x); }
  inline F64vec8 make_F64vec8(double x, double y, double z)
  { return F64vec8(0,z,y,x,0,z,y,x); }
  inline I32vec16 make_I32vec16(int x, int y, int z, int t)
  { return I32vec16(0,0,0,0,0,0,0,0,t,z,y,x,t,z,y,x); }
  inline I32vec16 make_I32vec16(int x, int y, int z)
  { return I32vec16(0,0,0,0,0,0,0,0,0,z,y,x,0,z,y,x); }

  // return concatenation of low half of lo and low half of hi
  inline __m512 cat_low_halves(__m512 lo, __m512 hi)
  {
    const __mmask16 rmask_00001111 = _mm512_int2mask(0xff00); //mask=11110000
    // low half uses lo, and high half uses low half of hi
    return _mm512_mask_permute4f128_ps(lo, rmask_00001111, hi,_MM_PERM_BADC);
  }
  
  // return concatenation of hgh half of lo and hgh half of hi
  inline __m512 cat_hgh_halves(__m512 lo, __m512 hi)
  {
    const __mmask16 rmask_11110000 = _mm512_int2mask(0x00ff); //mask=00001111
    // high half uses hi, and low half uses high half of lo
    return _mm512_mask_permute4f128_ps(hi, rmask_11110000, lo,_MM_PERM_BADC);
  }
  
  inline F64vec8 cat_low_halves(void* lo, void* hi)
  { return (F64vec8) cat_low_halves(*(__m512*)lo, *(__m512*)hi); }
  inline F64vec8 cat_hgh_halves(void* lo, void* hi)
  { return (F64vec8) cat_hgh_halves(*(__m512*)lo, *(__m512*)hi); }
  //
  inline F64vec8 cat_low_halves(F64vec8 lo, F64vec8 hi)
  { return (F64vec8) cat_low_halves(*(__m512*)(&lo), *(__m512*)(&hi)); }
  inline F64vec8 cat_hgh_halves(F64vec8 lo, F64vec8 hi)
  { return (F64vec8) cat_hgh_halves(*(__m512*)(&lo), *(__m512*)(&hi)); }
  
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
  
  // copy vectors of src to vectors of dst
  inline void copy012and456(F64vec8& dst, F64vec8 src)
  {
    const __mmask8 rmask_11101110 = _mm512_int2mask(0x77); //mask=01110111
    dst = _mm512_mask_mov_pd(dst,rmask_11101110,src);
  }
  
  // broadcast s0 into mask=0 slots, s1 into mask=1 slots
  inline F64vec8 broadcast_mask_blend(double s0, double s1, __mmask8 mask)
  {
    return _mm512_mask_blend_pd(mask, F64vec8(s0), F64vec8(s1));
    // broadcast s0 into vector
    //const F64vec8 t1 = _mm512_extload_pd(&s0, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
    // broadcast s1 into same vector with mask
    //return (F64vec8) _mm512_mask_extload_pd(t1, mask, &s1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
  }
  
  typedef I32vec16 Ivec;
  typedef F64vec8 Dvec;

//#elif defined(__AVX__)
//  // see http://software.intel.com/en-us/node/462716
//  #include <dvec.h>
//  typedef I32vec8 Ivec;
//  typedef F64vec4 Dvec;
//#elif defined(__SSE2__)
#else
  class Ivec8
  {
    int a[8];
  };
  class Dvec4
  {
    double a[4];
   public:
    Dvec4(double x, double y, double z, double t)
    { a[0]=x; a[1]=y; a[2]=z; a[3]=t; }
    Dvec4(double s)
    { Dvec4(s,s,s,s); }
  };
  inline Dvec4 initialize_vector(double x, double y, double z, double t)
  {
    return Dvec4(t,z,y,x);
  }
  inline Ivec8 round_down(Dvec4 dvec)
  {
    Ivec8 ret;
    for(int i=0;i<4;i++)
      ret[i] = floor(dvec[i]);
  }
  typedef Ivec8 Ivec;
  typedef Dvec4 Dvec;
#endif

#endif
