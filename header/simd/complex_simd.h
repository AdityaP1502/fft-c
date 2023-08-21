#ifndef _COMPLEX_SIMD_H_
#define _COMPLEX_SIMD_H_

#include <immintrin.h>
#include <math.h>
#include "../dll_export_api.h"

struct complex_array {
    double* real;
    double* imag;
};

struct complex_pair {
    double real;
    double imag;
};

typedef struct _reg_t {
    __m128d re, im;
} reg_t;

typedef struct complex_array complex;
typedef struct complex_pair complex_pair;

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// create complex array usign two double arrays for each double and real
// each array is malloc using __mm_alloc with 16 byte allignment
FFTLIBRARY_API complex* FFTLIBRARY_CALL complex_arr_create_allign_16(int size);

// destroy complex array that are created with 16 byte alllignment
FFTLIBRARY_API void FFTLIBRARY_CALL complex_arr_allign_16_destroy(complex* arr);

// load real and imag to register
FFTLIBRARY_API reg_t FFTLIBRARY_CALL complex_load_to_reg(double* re, double* im);

FFTLIBRARY_API void FFTLIBRARY_CALL complex_store_reg(double* re, double* im, reg_t a);

// do two simultanious addition 
// (a0 + jb0) + (c0 + jd0) and (a1 + jb1) + (c1 + jd1)
FFTLIBRARY_API reg_t FFTLIBRARY_CALL _mm_complexadd_pd(double* re1, double* im1, double* re2, double* im2);

// do two simultanious substraction
// (a0 + jb0) - (c0 + jd0) and (a1 + jb1) - (c1 + jd1)
FFTLIBRARY_API reg_t FFTLIBRARY_CALL _mm_complexsubs_pd(double* re1, double* im1, double* re2, double* im2);

// do two simultanious multiplication
// (a0 + jb0) * (c0 + jd0) and (a1 + jb1) * (c1 + jd1)
FFTLIBRARY_API reg_t FFTLIBRARY_CALL _mm_complexmul_pd(double* re1, double* im1, double* re2, double* im2);

// do two simultanious division
// (a0 + jb0) / (c0 + jd0) and (a1 + jb1) / (c1 + jd1)
FFTLIBRARY_API reg_t FFTLIBRARY_CALL _mm_complexdiv_pd(double* re1, double* im1, double* re2, double* im2);

FFTLIBRARY_API reg_t FFTLIBRARY_CALL _mm_complexadd_no_load_pd(reg_t a, reg_t b);

FFTLIBRARY_API reg_t FFTLIBRARY_CALL _mm_complexsubs_no_load_pd(reg_t a, reg_t b);

FFTLIBRARY_API reg_t FFTLIBRARY_CALL _mm_complexmul_no_load_pd(reg_t a, reg_t b);

#endif