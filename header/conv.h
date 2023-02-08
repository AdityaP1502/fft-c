#include <stdlib.h>
#include <math.h>

// convolve a and b
double* conv(double* a, double* b);

// convole a and b using fft
ifft_symmetric_bin* convfft(double* a, double* b, int length_a, int length_b);
