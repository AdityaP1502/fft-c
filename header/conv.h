#ifndef CONV_H
#define CONV_H

#include "fft.h"
#include "dll_export_api.h"

// convole a and b using fft
FFTLIBRARY_API ifft_symmetric_bins* FFTLIBRARY_CALL convfft(double* a, double* b, int length_a, int length_b);

FFTLIBRARY_API ifft_symmetric_bins* FFTLIBRARY_CALL convfft_overlap_save(double* a, double* b, int sample_length);

FFTLIBRARY_API ifft_symmetric_bins* FFTLIBRARY_CALL convfft_static_n(double* a, double* b, int length_a, int length_b, bins forward_twid_factor, bins backward_twid_factor, int fft_size);

FFTLIBRARY_API ifft_symmetric_bins* FFTLIBRARY_CALL convfft_overlap_save_static_n(double* a, double* b, bins forward_twid_factor, bins backward_twid_factor, int fft_size);

FFTLIBRARY_API ifft_symmetric_bins* FFTLIBRARY_CALL convfft4_static_n(double* a, double* b, int length_a, int length_b, bins forward_twid_factor, bins backward_twid_factor, int fft_size);

FFTLIBRARY_API ifft_symmetric_bins* FFTLIBRARY_CALL convfft4_overlap_save_static_n(double* a, double* b, bins forward_twid_factor, bins backward_twid_factor, int fft_size);
#endif
