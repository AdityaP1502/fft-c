#ifndef ITERATIVE_FFT_H
#define ITERATIVE_FFT_H

#include "fft.h"

/* calculate fft of xn (xn and the result of the fft is ordered)
*order is the implementation type of the fft iterartive
*order can be "RN" which is reversed input and natural output
*order can be "NR" which is natural input and reversed output
*order can be "NN" which is natural iput and natural output 
*/
fft_bins* fft_iterative(double* xn, int length, char* order);

/*
*calculate ifft of xn (xn and the result of the fft is ordered)
*order is the implementation type of the fft iterartive
*order can be "RN" which is reversed input and natural output
*order can be "NR" which is natural input and reversed output
*order can be "NN" which is natural iput and natural output
*/
fft_bins* ifft_iterative(bins xk, int length, char* order);

/*
*calculate symmetric ifft of xn (xn and the result of the fft is ordered)
*order is the implementation type of the fft iterartive
*order can be "RN" which is reversed input and natural output
*order can be "NR" which is natural input and reversed output
*order can be "NN" which is natural iput and natural output 
*/
ifft_symmetric_bins* ifft_iterative_symmetric(bins Xk, int length, char* order);

/*
Calculate fft of two real signal at a price of one
*/
fft_bins** fft_double_real(double* xn_1, double* xn_2, int length_1, int length_2, char* order);

#endif