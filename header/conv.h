#ifndef CONV_H
#define CONV_H

#include "fft.h"

// convolve a and b
double* conv(double* a, double* b);

// convole a and b using fft
ifft_symmetric_bins* convfft(double* a, double* b, int length_a, int length_b);

ifft_symmetric_bins* convfft_overlap_save(double* a, double* b, int sample_length);

#endif
