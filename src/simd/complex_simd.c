#include "../../header/complex_simd.h"

complex* complex_arr_create_allign_16(int size)
{
    complex* complex_arr;
    double* re;
    double* im;
    
    complex_arr = malloc( sizeof( complex ) );

    re = _mm_malloc(size * sizeof( double ), 16);
    im = _mm_malloc(size * sizeof( double ), 16);

    complex_arr->real = re;
    complex_arr->imag = im;

    return complex_arr;
};

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
