#ifndef ITERATIVE_FFT_H
#define ITERATIVE_FFT_H

#include "fft.h"
#include "dll_export_api.h"

/* calculate fft of xn (xn and the result of the fft is ordered)
*order is the implementation type of the fft iterartive
*order can be "RN" which is reversed input and natural output
*order can be "NR" which is natural input and reversed output
*order can be "NN" which is natural iput and natural output 
*/
FFTLIBRARY_API fft_bins* FFTLIBRARY_CALL fft_iterative(double* xn, int length, char* order);

/*
*calculate ifft of xn (xn and the result of the fft is ordered)
*order is the implementation type of the fft iterartive
*order can be "RN" which is reversed input and natural output
*order can be "NR" which is natural input and reversed output
*order can be "NN" which is natural iput and natural output
*/
FFTLIBRARY_API fft_bins* FFTLIBRARY_CALL ifft_iterative(bins xk, int length, char* order);

/*
*calculate symmetric ifft of xn (xn and the result of the fft is ordered)
*order is the implementation type of the fft iterartive
*order can be "RN" which is reversed input and natural output
*order can be "NR" which is natural input and reversed output
*order can be "NN" which is natural iput and natural output 
*/
FFTLIBRARY_API ifft_symmetric_bins* FFTLIBRARY_CALL ifft_iterative_symmetric(bins Xk, int length, char* order);

/*
Calculate fft of two real signal at a price of one
*/
FFTLIBRARY_API fft_bins** FFTLIBRARY_CALL fft_double_real(double* xn_1, double* xn_2, int length_1, int length_2, char* order);


// Static N FFT Iterative
/*  Use this function when calculating FFT with a fixed size
    so no need for precomputing each fft calls. 
*/

FFTLIBRARY_API fft_bins* FFTLIBRARY_CALL fft_iterative_static_n(bins forward_twid_factors, int fft_size, double* xn, int length, char* order);

FFTLIBRARY_API fft_bins* FFTLIBRARY_CALL ifft_iterative_static_n(bins xk, bins backward_twid_factor, int fft_size, int length, char* order);

FFTLIBRARY_API ifft_symmetric_bins* FFTLIBRARY_CALL ifft_iterative_symmetric_static_n(bins Xk, bins backward_twid_factor, int fft_size, int length, char* order);

FFTLIBRARY_API fft_bins** FFTLIBRARY_CALL fft_double_real_static_n(bins forward_twiddle_factors, double* xn_1, double* xn_2, int length_1, int length_2, int fft_size, char* order);
#endif