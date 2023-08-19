#include <stdio.h>
#include "../header/simd/fft_simd.h"

void test_for_twiddle_forward()
{
    int size = 8;
    complex* twid = precompute_twiddle_factor(size, 0);

    // expected
    // 1.00000 0.00000
    // 0.70711 -0.70711
    // 0.00000 -1.00000
    // -0.70711 -0.70711

    for (int i = 0; i < size >> 1; i++)
    {
        printf("%f + j%f\n", twid->real[i], twid->imag[i]);
    }
    
    // received 
    // 1.000000 + j0.000000
    // 0.707107 + j-0.707107
    // 0.000000 + j-1.000000
    // -0.707107 + j-0.707107
    complex_arr_allign_16_destroy(twid);
}

void test_for_twiddle_backward()
{
    int size = 8;
    complex* twid = precompute_twiddle_factor(size, 1);

    // expected
    // 1.00000 0.00000
    // 0.70711 -0.70711
    // 0.00000 -1.00000
    // -0.70711 -0.70711

    for (int i = 0; i < size >> 1; i++)
    {
        printf("%f + j%f\n", twid->real[i], twid->imag[i]);
    }
    
    // received 
    // 1.000000 + j0.000000
    // 0.707107 + j-0.707107
    // 0.000000 + j-1.000000
    // -0.707107 + j-0.707107
    complex_arr_allign_16_destroy(twid);
}

void test_for_twiddle_forward_r4()
{
    int size = 8;
    complex* twid = precompute_twiddle_factor_radix_4(size, 0);

    // expected
    // 1.00000 0.00000
    // 0.70711 -0.70711
    // 0.00000 -1.00000
    // -0.70711 -0.70711

    for (int i = 0; i < (3 * size >> 2); i++)
    {
        printf("%f + j%f\n", twid->real[i], twid->imag[i]);
    }
    
    // received 
    // 1.000000 + j0.000000
    // 0.707107 + j-0.707107
    // 0.000000 + j-1.000000
    // -0.707107 + j-0.707107
    complex_arr_allign_16_destroy(twid);
}

void test_for_twiddle_backward_r4()
{
    int size = 8;
    complex* twid = precompute_twiddle_factor_radix_4(size, 1);

    // expected
    // 1.00000 0.00000
    // 0.70711 -0.70711
    // 0.00000 -1.00000
    // -0.70711 -0.70711

    for (int i = 0; i < (3 * size >> 2); i++)
    {
        printf("%f + j%f\n", twid->real[i], twid->imag[i]);
    }
    
    // received 
    // 1.000000 + j0.000000
    // 0.707107 + j-0.707107
    // 0.000000 + j-1.000000
    // -0.707107 + j-0.707107
    complex_arr_allign_16_destroy(twid);
}
int main()
{
    printf("Radix 2 F\n");
    test_for_twiddle_forward();
    printf("Radix 2 B\n");
    test_for_twiddle_backward();
    printf("Radix 4 F\n");
    test_for_twiddle_forward_r4();
    printf("Radix 4 B\n");
    test_for_twiddle_backward_r4();
    return 0;
}