#include "fft.h"

// calculate fft of xn
fft_bin* fft_recursive(double* xn, int length);

// calculate inverse fft
fft_bin* ifft_recursive(bins Xk, int length);

// calculate inverse fft. Xk is symetric
ifft_symmetric_bin* ifft_recursive_symmetric(bins Xk, int length);

