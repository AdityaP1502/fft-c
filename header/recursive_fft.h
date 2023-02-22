#ifndef RECURSIVE_FFT_H
#define RECURSIVE_FFT_H

#include "fft.h"

// calculate fft of xn
fft_bins* fft_recursive(double* xn, int length);

// calculate inverse fft
fft_bins* ifft_recursive(bins Xk, int length);

// calculate inverse fft. Xk is symetric
ifft_symmetric_bins* ifft_recursive_symmetric(bins Xk, int length);

#endif