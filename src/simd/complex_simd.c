#include "../../header/simd/complex_simd.h"

complex* complex_arr_create_allign_16(int size)
{
    complex* complex_arr;
    double* re;
    double* im;
    
    complex_arr = malloc( sizeof( complex ) );

    re = _mm_malloc(size * sizeof( double ), 16);
    im = _mm_malloc(size * sizeof( double ), 16);

    for (int i = 0; i < size; i++)
    {
        re[i] = 0;
        im[i] = 0;
    }

    complex_arr->real = re;
    complex_arr->imag = im;

    return complex_arr;
};

reg_t complex_load_to_reg(double* re, double* im)
{
    reg_t a;

    a.re = _mm_load_pd(re);
    a.im = _mm_load_pd(im);

    return a;
}

void complex_store_reg(double* re, double* im, reg_t a)
{
    _mm_store_pd(re, a.re);
    _mm_store_pd(im, a.im);
}

void complex_arr_allign_16_destroy(complex* arr)
{
    _mm_free(arr->real);
    _mm_free(arr->imag);
    free(arr);
};

reg_t _mm_complexadd_pd(double* re1, double* im1, double* re2, double* im2)
{
    __m128d n1;
    __m128d n2;
    __m128d n3;
    __m128d n4;
    reg_t res;

    // grabbed two double values from re1 => a0 and a1
    n1 = _mm_loadu_pd(re1);
    n2 = _mm_loadu_pd(im1);
    n3 = _mm_loadu_pd(re2);
    n4 = _mm_loadu_pd(im2);
    
    // add n1 with n3
    res.re = _mm_add_pd(n1, n3);
    res.im = _mm_add_pd(n2, n4);

    return res;
}

reg_t _mm_complexsubs_pd(double* re1, double* im1, double* re2, double* im2)
{
    __m128d n1;
    __m128d n2;
    __m128d n3;
    __m128d n4;
    reg_t res;

    // grabbed two double values from re1 => a0 and a1
    n1 = _mm_loadu_pd(re1);
    n2 = _mm_loadu_pd(im1);
    n3 = _mm_loadu_pd(re2);
    n4 = _mm_loadu_pd(im2);
    
    // add n1 with n3
    res.re = _mm_sub_pd(n1, n3);
    res.im = _mm_sub_pd(n2, n4);

    return res;
}

reg_t _mm_complexmul_pd(double* re1, double* im1, double* re2, double* im2)
{
    __m128d n1;
    __m128d n2;
    __m128d n3;
    __m128d n4;

    reg_t res;

    // grabbed two double values from re1 => a0 and a1
    n1 = _mm_loadu_pd(re1);
    n2 = _mm_loadu_pd(im1);
    n3 = _mm_loadu_pd(re2);
    n4 = _mm_loadu_pd(im2);
    
    // a0c0 | a1c1
    // b0d0 | b1d1
    // a0c0 - b0d0 | a1c1 - b1d1
    res.re = _mm_sub_pd( _mm_mul_pd( n1, n3 ), _mm_mul_pd( n2, n4 ) ); 

    // a0d0 | a1d1
    // b0c0 | b1c1
    // a0d0 + b0c0 | a1d1 + b1c1
    res.im = _mm_add_pd( _mm_mul_pd( n1, n4 ), _mm_mul_pd( n2, n3 ) );
    return res;
}

reg_t _mm_complexdiv_pd(double* re1, double* im1, double* re2, double* im2)
{
    __m128d n1, n2, n3, n4, n3_c, n4_c, t;
    reg_t res;

    // grabbed two double values from re1 => a0 and a1
    n1 = _mm_loadu_pd(re1);
    n2 = _mm_loadu_pd(im1);

    n3 = _mm_loadu_pd(re2);
    n4 = _mm_loadu_pd(im2);

    n3_c = _mm_loadu_pd(re2);
    n4_c = _mm_loadu_pd(im2);

    t = _mm_add_pd( _mm_mul_pd( n3, n3_c ), _mm_mul_pd( n4, n4_c ) );

    // a0c0 | a1c1
    // b0d0 | b1d1
    // a0c0 + b0d0 | a1c1 + b1d1
    // aoco + b0d0 / co^2 + d0^2 | a1c1 + b1d1 / c1 ^ 2 + d1 ^ 2
    res.re = _mm_div_pd( _mm_add_pd( _mm_mul_pd( n1, n3 ), _mm_mul_pd( n2, n4 ) ), t); 

    // a0d0 | a1d1
    // b0c0 | b1c1
    // b0c0 - a0d0 | b1c1 - a1d1
    // b0c0 - a0c0 / co^2 + d0^2 | b1c1 - a1d1 / c1 ^ 2 + d1 ^ 2
    res.im = _mm_div_pd( _mm_sub_pd( _mm_mul_pd( n2, n3 ), _mm_mul_pd( n1, n4 ) ), t );
    
    return res;
}

reg_t _mm_complexadd_no_load_pd(reg_t a, reg_t b)
{
    reg_t res;

    // add n1 with n3
    res.re = _mm_add_pd(a.re, b.re);
    res.im = _mm_add_pd(a.im, b.im);

    return res;
}

reg_t _mm_complexsubs_no_load_pd(reg_t a, reg_t b)
{
    reg_t res;

    res.re = _mm_sub_pd(a.re, b.re);
    res.im = _mm_sub_pd(a.im, b.im);

    return res;
}

reg_t _mm_complexmul_no_load_pd(reg_t a, reg_t b)
{
    reg_t res;

    res.re = _mm_sub_pd( _mm_mul_pd( a.re, b.re ), _mm_mul_pd( a.im, b.im ) ); 
    res.im = _mm_add_pd( _mm_mul_pd( a.re, b.im ), _mm_mul_pd( a.im, b.re ) );

    return res;
}

reg_t _mm_complex_mulj_no_load_pd(reg_t a)
{
    reg_t res;

    res.re = _mm_add_pd(_mm_sub_pd(a.re, a.im), a.re);
    res.im = _mm_sub_pd(_mm_add_pd(a.re, a.im), a.im);

    return res;
}

////////////////////////////////// NORMAL COMPLEX OPERATION ///////////////////////////////////////////

complex_pair FFTLIBRARY_CALL complexadd(double re1, double re2, double im1, double im2)
{
    complex_pair r;
    r.real = re1 + re2;
    r.imag = im1 + im2;
    return r;
}

complex_pair FFTLIBRARY_CALL complexsub(double re1, double re2, double im1, double im2)
{
    complex_pair r;
    r.real = re1 - re2;
    r.imag = im1 - im2;
    return r;
}

complex_pair FFTLIBRARY_CALL complexmul(double re1, double re2, double im1, double im2)
{
    complex_pair r;

    r.real = re1 * re2 - im1 * im2;
    r.imag = re1 * im2 + re2 * im1;

    return r;
}

complex_pair FFTLIBRARY_CALL complexdiv(double re1, double re2, double im1, double im2)
{
    complex_pair r;
    double t;

    t = re2 * re2 + im2 * im2;
    
    r.real = (re1 * re2 + im1 * im2) / t;
    r.imag = (re2 * im1 - re1 * im2) / t;

    return r;
}

complex_pair FFTLIBRARY_CALL complexmul_j(double re, double im)
{
    complex_pair r;

    r.real = -im;
    r.imag = re;

    return r;
}

complex_pair FFTLIBRARY_CALL complexmul_minj(double re, double im)
{
    complex_pair r;

    r.real = im;
    r.imag = -re;

    return r;
}