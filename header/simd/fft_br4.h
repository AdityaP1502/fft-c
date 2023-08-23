#ifndef _FFT_BR4_H_
#define _FFT_BR4_H_

#include "fft_simd.h"

void fft_simd_br4_f(complex* c_arr, complex* W, int j0, int k, int ql, int hl);

void fft_simd_br4_b(complex* c_arr, complex* W, int j0, int k, int ql, int hl);

void fft_br4_f(complex* c_arr, complex* W, int j0, int k, int ql, int hl);

void fft_br4_b(complex* c_arr, complex* W, int j0, int k, int ql, int hl);

void fft_br4_0(complex* c_arr, int j0, int ql, int hl);

void fft_br4_f_N_over_8(complex* c_arr, int j0, int ql, int hl);

void fft_br4_b_N_over_16(complex* c_arr, int j0, int ql, int hl);

void fft_br4_f_N_over_16(complex* c_arr, complex* W, int j0, int ql, int hl, int k);

void fft_br4_b_N_over_8(complex* c_arr, complex* W, int j0, int ql, int hl, int k);

void fft_br4_f_3N_over_16(complex* c_arr, complex* W, int j0, int ql, int hl, int k);

void fft_br4_b_3N_over_16(complex* c_arr, complex* W, int j0, int ql, int hl, int k);

#endif
