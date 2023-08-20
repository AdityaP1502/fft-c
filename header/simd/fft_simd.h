#ifndef _FFT_SIMD_H_
#define _FFT_SIMD_H_

#include "complex_simd.h"

FFTLIBRARY_API complex* FFTLIBRARY_CALL precompute_twiddle_factor(int length, int backward);

FFTLIBRARY_API complex* FFTLIBRARY_CALL precompute_twiddle_factor_radix_4(int length, int backward);

FFTLIBRARY_API void FFTLIBRARY_CALL fill_complex_arr_real(complex* c_arr, double* re, int size);

FFTLIBRARY_API void FFTLIBRARY_CALL fill_complex_arr_complex(complex* c_arr, double* re, double* im, int sizeR, int sizeIm);

FFTLIBRARY_API void FFTLIBRARY_CALL normalize_ifft(complex* c_arr, int size);

FFTLIBRARY_API void FFTLIBRARY_CALL seperate_combined_output(complex* src, complex* dest1, complex* dest2, int length);

#endif