#include "fft.h"
#include "dll_export_api.h"

/*  Calculate FFT of xn using radix 4 fft. Xn is ordered and the result is ordered
    The order of calculation will be using the NR mode.  
*/
FFTLIBRARY_API fft_bins* FFTLIBRARY_CALL fft_radix_4_iterative_static_n(bins forward_twid_factors, int fft_size, double* xn, int length);

/*  Calculate IFFT of Xk using radix 4 fft. Xk is ordered and the result is ordered
    The order of calculation will be using the NR mode.  
*/
FFTLIBRARY_API fft_bins* FFTLIBRARY_CALL ifft_radix_4_iterative_static_n(bins backward_twid_factor, int fft_size, bins xk, int length);

/*  Calculate IFFT of xn using radix 4 fft. Xk is ordered and the result is ordered.
    The order of calculation will be using the NR mode.  
    The result use the assumption that Xk is symmetric, therefore the imaginary part of the answer is discarded
*/
FFTLIBRARY_API ifft_symmetric_bins* FFTLIBRARY_CALL ifft_radix_4_iterative_symmetric_static_n(bins backward_twid_factor, int fft_size, bins xk, int length);

/*  Calculate FFT of xn using radix 4 fft. Xn is ordered and the result is ordered
    The order of calculation will be using the NR mode.  
    Calculate two fft at the price of one. 
*/
FFTLIBRARY_API fft_bins** FFTLIBRARY_CALL fft_radix_4_double_real_static_n(bins forward_twiddle_factors, double* xn_1, double* xn_2, int length_1, int length_2, int fft_size);