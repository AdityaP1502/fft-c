#include <stdio.h>
#include "../header/simd/fftr2_simd.h"

void fft_r2_test()
{
    int size = 8;
    double arr[] = { 1, 3, 1, 1, 1, 1, 1, 1 };
    complex* w = precompute_twiddle_factor(size, 0);
    complex* s = fftr_r2(arr, w, size, size);

    for (int i = 0; i < size; i++)
    {
        printf("%f %f\n", s->real[i], s->imag[i]);
    }

    complex_arr_allign_16_destroy(w);
    complex_arr_allign_16_destroy(s);
}

void ifft_r2_test()
{
    int size = 8;
    double arr[] = { 1, 3, 1, 1, 1, 1, 1, 1 };

    complex* w = precompute_twiddle_factor(size, 0);
    complex* wb = precompute_twiddle_factor(size, 1);
    complex* s = fftr_r2(arr, wb, size, size);

    ifftr_r2_inplace(s, w, size);

    for (int i = 0; i < size; i++)
    {
        printf("%f %f\n", s->real[i], s->imag[i]);
    }
    
    complex_arr_allign_16_destroy(w);
    complex_arr_allign_16_destroy(wb);
    complex_arr_allign_16_destroy(s);
}
int main()
{
    fft_r2_test();
    ifft_r2_test();
    return 0;
}