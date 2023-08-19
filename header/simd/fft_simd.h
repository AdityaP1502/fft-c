#include "complex_simd.h"

complex* FFTLIBRARY_CALL precompute_twiddle_factor(int length, int backward);

complex* FFTLIBRARY_CALL precompute_twiddle_factor_radix_4(int length, int backward);