#include "fft_simd.h"

// fft radix 4 
complex* fftr_r4(double* r, complex* W, int arr_size, int fft_size);