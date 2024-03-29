#ifndef FFT_H
#define FFT_H

#include "complex.h"
#include "dll_export_api.h"

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
FFTLIBRARY_API void FFTLIBRARY_CALL destroy_bin(bins bin, int length);

// convert complex array to real array
// by omitting the imaginary part
FFTLIBRARY_API double* FFTLIBRARY_CALL convert_complex_to_real(bins complex_array, int length);

// convert real array to complex array
FFTLIBRARY_API bins FFTLIBRARY_CALL convert_real_to_complex(double* real_array, int length);

// find the nearest power of 2 of
// a given number
FFTLIBRARY_API int FFTLIBRARY_CALL nearest_power_of_2(int number);

// find the nearest power of 4 of
// a given number
FFTLIBRARY_API int FFTLIBRARY_CALL nearest_power_of_4(int number);

// copy an array of complex of l = length to a new array
// if pad_length > 0, then the array will be padded with zero
FFTLIBRARY_API bins FFTLIBRARY_CALL deep_copy_bins_complex(bins xk, int length, int pad_length);

// copy an array of complex of l = length to a new array
// if pad_length > 0, then the array will be padded with zero
FFTLIBRARY_API bins FFTLIBRARY_CALL copy_bins_complex(bins xk, int length, int pad_length);

// copy an array of complex of l = length to a new array
// if pad_length > 0, then the array will be padded with zero
FFTLIBRARY_API double* FFTLIBRARY_CALL copy_bins_real(double* xn, int length, int pad_length);

// clear pad_length amount of zeros at the end of the array
FFTLIBRARY_API void FFTLIBRARY_CALL clear_pad(bins xk, int length, int pad_length);

// clear pad_length amount of zeros at the end of the array
FFTLIBRARY_API double* FFTLIBRARY_CALL clear_pad_real(double* xn, int length);

// precompute the twiddle factor (for length in power of 2)
// if backward is 0 then calculate the twiddle factor 
// for forward fft else for reverse fft (ifft)
// the output will always has element equal to half of the length passed
FFTLIBRARY_API bins FFTLIBRARY_CALL precompute_twiddle_factor(int length, int backward);

// precompute the twiddle factor (for length in power of 2)
// if backward is 0 then calculate the twiddle factor 
// for forward fft else for reverse fft (ifft)
// the output will always has element equal to half of the length passed
FFTLIBRARY_API bins FFTLIBRARY_CALL precompute_twiddle_factor_radix_4(int length, int backward);

// combined two real input into one complex output
FFTLIBRARY_API bins FFTLIBRARY_CALL combined_two_real_input(double *xn_1, double *xn_2, int max_length, int min_length);

// Seperated combined output from double FFt 
// The output will be placed in dest 
FFTLIBRARY_API void FFTLIBRARY_CALL seperate_combined_output(bins combined_output, int length, fft_bins **dest);

// calculate the magnitude spectrum
FFTLIBRARY_API double* FFTLIBRARY_CALL frequency_spectrum_magnitude(bins Xk, int length);

// calculate the phase spectrum
FFTLIBRARY_API double* FFTLIBRARY_CALL frequency_spectrum_phase(bins Xk, int length);

FFTLIBRARY_API void FFTLIBRARY_CALL normalize_ifft(bins xk, int length);

#endif