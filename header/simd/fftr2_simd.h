#ifndef _FFTR2_SIMD_H_
#define _FFTR2_SIMD_H_

#include "fft_simd.h"

complex* fftr_r2(double* r, complex* w, int arr_size, int fft_size);

complex* ifftr_r2(complex* c_arr, complex* w, int arr_size, int fft_size);

void ifftr_r2_inplace(complex* c_arr, complex* w, int fft_size);

complex** fftdr_r2(double* r1, double* r2, complex* w, int size_1, int size_2, int fft_size);

#endif