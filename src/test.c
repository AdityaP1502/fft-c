/* Testing modules*/

#include "../header/complex.h"
#include "../header/iterative_fft.h"
#include "../header/conv.h"

int test_for_ifft(double* xk, int length)
{
    char* a;
    fft_bins* bin;
    fft_bins* res;

    bin = fft_iterative(xk, length, "NR");

    printf("Iterative FFT in using NR\n");

    for (int i = 0; i < bin->length; i++)
    {
        a = complex_number_to_string(bin->fft_bins[i]);
        printf("%s\n", a);
        free(a);
    }

    printf("ifft result\n");

    res = ifft_iterative(bin->fft_bins, bin->length, "NR");

    for (int i = 0; i < bin->length; i++)
    {
        a = complex_number_to_string(res->fft_bins[i]);
        printf("%s\n", a);
        free(a);
    }

    destroy_bin(bin->fft_bins, bin->length);
    free(bin);

    destroy_bin(res->fft_bins, res->length);
    free(res);

    return 0;
}

int test_for_symmetric_ifft(double* xk, int length)
{
    char* a;
    fft_bins* bin;
    ifft_symmetric_bins* res;

    bin = fft_iterative(xk, length, "NR");

    printf("Iterative FFT in using NR\n");

    for (int i = 0; i < bin->length; i++)
    {
        a = complex_number_to_string(bin->fft_bins[i]);
        printf("%s\n", a);
        free(a);
    }

    printf("ifft result\n");

    res = ifft_iterative_symmetric(bin->fft_bins, bin->length, "NR");

    for (int i = 0; i < res->length; i++)
    {
        printf("%f\n", res->bin[i]);
    }

    destroy_bin(bin->fft_bins, bin->length);
    free(bin);

    free(res->bin);
    free(res);

    return 0;
}

int test_for_double_real_fft(double* a, double* b, int length_a, int length_b)
{
    char* string;
    fft_bins* a_bin;
    fft_bins* b_bin;
    fft_bins** res;

    res = fft_double_real(a, b, length_a, length_b, "NR");

    a_bin = res[0];
    b_bin = res[1];

    for (int i = 0; i < a_bin->length; i++)
    {
        string = complex_number_to_string(a_bin->fft_bins[i]);
        printf("%s\n", string);
        free(string);
    }

    destroy_bin(a_bin->fft_bins, a_bin->length);
    free(a_bin);

    for (int i = 0; i < b_bin->length; i++)
    {
        string = complex_number_to_string(b_bin->fft_bins[i]);
        printf("%s\n", string);
        free(string);
    }

    destroy_bin(b_bin->fft_bins, b_bin->length);
    free(b_bin);

    free(res);

    return 0;
}

int test_for_conv(double* a, double* b, int length_a, int length_b)
{
    ifft_symmetric_bins* res;

    res = convfft(a, b, length_a, length_b);
    
    for (int i = 0; i < res->length; i++)
    {
        printf("%f\n", res->bin[i]);
    }

    free(res->bin);
    free(res);

    return 0;
}

int test_for_overlap_save(double* a, double* b, int sample_length)
{
	ifft_symmetric_bins* res;
	res = convfft_overlap_save(a, b, sample_length);
	for (int i = 0; i < res->length; i++) 
	{
		printf("%f\n", res->bin[i]);
	}

	free(res->bin);
	free(res);
	
	return 0;
}
int main() 
{
    int length_2, length_1;
    // fft_bin* ifft_bin;

    // length = 8;
    // backward = 1;

    // bins twiddle_factors = precompute_twiddle_factor(length, backward);
    // for (int i = 0; i < length >> 1; i++)
    // {
    //     complex_number* a = twiddle_factors[i];
    //     char* b = complex_number_to_string(a);
    //     printf("%s\n", b);
    //     free(b);
    // }

    // destroy_bin(twiddle_factors, length >> 1);

    // int length = 32;
    // for (int i = 0; i < length; i++)
    // {
    //     printf("%d,%d\n", i, reverse_bit(i, length));
    // }

    double xk_1[] = { 2, 2, 2, 2 };
    length_1 = sizeof(xk_1) / sizeof(xk_1[0]);

    double xk_2[] = { 2, 3, 3, 3 };
    length_2 =  sizeof(xk_2) / sizeof(xk_2[0]);

    //test_for_ifft(xk_2, length_2);
    // test_for_symmetric_ifft(xk_1, length_1);
    // test_for_double_real_fft(xk_1, xk_2, length_1, length_2);
    // test_for_conv(xk_1, xk_2, length_1, length_2);
    test_for_overlap_save(xk_1, xk_2, length_1);

    return 0;
}
