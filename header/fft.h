#ifndef FFT_H
#define FFT_H

#include "complex.h"

typedef complex_number** bins;

typedef struct fft_result {
  bins fft_bins;
  int length;
} fft_bins;

typedef struct ifft_symmetric_result {
  double* bin;
  int length;
} ifft_symmetric_bins;

// free all data stored in bin
void destroy_bin(bins bin, int length);

// convert complex array to real array
// by omitting the imaginary part
double* convert_complex_to_real(bins complex_array, int length);

// convert real array to complex array
bins convert_real_to_complex(double* real_array, int length);

// find the nearest power of 2 of
// the given number
int nearest_power_of_2(int number);

// copy an array of complex of l = length to a new array
// if pad_length > 0, then the array will be padded with zero
bins deep_copy_bins_complex(bins xk, int length, int pad_length);

// copy an array of complex of l = length to a new array
// if pad_length > 0, then the array will be padded with zero
bins copy_bins_complex(bins xk, int length, int pad_length);

// copy an array of complex of l = length to a new array
// if pad_length > 0, then the array will be padded with zero
double* copy_bins_real(double* xn, int length, int pad_length);

// clear pad_length amount of zeros at the end of the array
void clear_pad(bins xk, int length, int pad_length);

// clear pad_length amount of zeros at the end of the array
double* clear_pad_real(double* xn, int length);

// precompute the twiddle factor (for length in power of 2)
// if forward is 1 then calculate the twiddle factor 
// for forward fft else for reverse fft (ifft)
// the output will always has element equal to half of the length passed
bins precompute_twiddle_factor(int length, int backward);

// calculate the magnitude spectrum
double* frequency_spectrum_magnitude(bins Xk, int length);

// calculate the phase spectrum
double* frequency_spectrum_phase(bins Xk, int length);

#endif