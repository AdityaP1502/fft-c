/* Testing modules*/

#include "../header/complex.h"
#include "../header/iterative_fft.h"
#include "../header/fft_radix_4.h"
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


int test_for_ifft_static_n(double* xk, int length, bins forward_twid_factor, int fft_size) {
    char* a;
    fft_bins* bin;
    fft_bins* res;


    bin = fft_iterative_static_n(forward_twid_factor, fft_size, xk, length, "NR");
    
    printf("Iterative FFT with fft size %d\n", fft_size);

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

int test_for_symmetric_ifft_static_n(double* xk, int length)
{
    char* a;
    fft_bins* bin;
    ifft_symmetric_bins* res;
    int fft_size;
    bins backward_twid_factor;

    bin = fft_iterative(xk, length, "NR");

    printf("Iterative FFT in using NR\n");

    for (int i = 0; i < bin->length; i++)
    {
        a = complex_number_to_string(bin->fft_bins[i]);
        printf("%s\n", a);
        free(a);
    }

    fft_size = bin->length;
    backward_twid_factor = precompute_twiddle_factor(fft_size, 0);

    printf("ifft result\n");

    res = ifft_iterative_symmetric_static_n(bin->fft_bins, backward_twid_factor, fft_size, bin->length, "NR");

    for (int i = 0; i < res->length; i++)
    {
        printf("%f\n", res->bin[i]);
    }

    destroy_bin(bin->fft_bins, bin->length);
    free(bin);

    free(res->bin);
    free(res);

    destroy_bin(backward_twid_factor, fft_size >> 1);


    return 0;
}

int test_for_double_real_fft_static_n(double* a, double* b, int length_a, int length_b, int fft_size, bins forward_twid_factor)
{
    char* string;
    fft_bins* a_bin;
    fft_bins* b_bin;
    fft_bins** res;

    res = fft_double_real_static_n(forward_twid_factor, a, b, length_a, length_b, fft_size, "NR");

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

int test_for_conv_static_n(double* a, double* b, int length_a, int length_b, bins forward_twid_factor, bins backward_twid_factor, int fft_size)
{
    ifft_symmetric_bins* res;

    res = convfft_static_n(a, b, length_a, length_b, forward_twid_factor, backward_twid_factor, fft_size);
    
    for (int i = 0; i < res->length; i++)
    {
        printf("%f\n", res->bin[i]);
    }

    free(res->bin);
    free(res);

    return 0;
}

int test_for_overlap_save_static_n(double* a, double* b, int sample_length, bins forward_twid_factor, bins backward_twid_factor)
{
	ifft_symmetric_bins* res;
	res = convfft_overlap_save_static_n(a, b, forward_twid_factor, backward_twid_factor, sample_length);

	for (int i = 0; i < res->length; i++) 
	{
		printf("%f\n", res->bin[i]);
	}

	free(res->bin);
	free(res);
	
	return 0;
}

int test_for_it_radix_4(double* a, int length, int sample_length, bins forward_twid_factor)
{
    fft_bins* res;
    char* string;

    res = fft_radix_4_iterative_static_n(forward_twid_factor, sample_length, a, length);

    for (int i = 0; i < res->length; i++)
    {
        string = complex_number_to_string(res->fft_bins[i]);
        printf("%s\n", string);
        free(string);
    }

    destroy_bin(res->fft_bins, res->length);
    free(res);

    return 0;
}
int test_for_ifft_radix_4_static_n(double* xk, int length, bins forward_twid_factor, bins backward_twid_factor, int fft_size) {
    char* a;
    fft_bins* bin;
    fft_bins* res;


    bin = fft_radix_4_iterative_static_n(forward_twid_factor, fft_size, xk, length);
    
    printf("Iterative FFT with fft size %d\n", fft_size);

    for (int i = 0; i < bin->length; i++)
    {
        a = complex_number_to_string(bin->fft_bins[i]);
        printf("%s\n", a);
        free(a);
    }

    printf("ifft result\n");

    res = ifft_radix_4_iterative_static_n(backward_twid_factor, fft_size, bin->fft_bins, bin->length);

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

int test_for_symmetric_ifft_radix_4_static_n(double* xk, int length, bins forward_twid_factor, bins backward_twid_factor, int fft_size)
{
    char* a;
    fft_bins* bin;
    ifft_symmetric_bins* res;

    bin = fft_radix_4_iterative_static_n(forward_twid_factor, fft_size, xk, length);

    printf("Iterative FFT in using NR\n");

    for (int i = 0; i < bin->length; i++)
    {
        a = complex_number_to_string(bin->fft_bins[i]);
        printf("%s\n", a);
        free(a);
    }

    fft_size = bin->length;
    printf("ifft result\n");

    res = ifft_radix_4_iterative_symmetric_static_n(backward_twid_factor, fft_size, bin->fft_bins, bin->length);

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

int test_for_double_real_fft_radix_4_static_n(double* a, double* b, int length_a, int length_b, int fft_size, bins forward_twid_factor)
{
    char* string;
    fft_bins* a_bin;
    fft_bins* b_bin;
    fft_bins** res;

    res = fft_radix_4_double_real_static_n(forward_twid_factor, a, b, length_a, length_b, fft_size);

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

// int test_for_conv_static_n(double* a, double* b, int length_a, int length_b, bins forward_twid_factor, bins backward_twid_factor, int fft_size)
// {
//     ifft_symmetric_bins* res;

//     res = convfft_static_n(a, b, length_a, length_b, forward_twid_factor, backward_twid_factor, fft_size);
    
//     for (int i = 0; i < res->length; i++)
//     {
//         printf("%f\n", res->bin[i]);
//     }

//     free(res->bin);
//     free(res);

//     return 0;
// }

// int test_for_overlap_save_static_n(double* a, double* b, int sample_length, bins forward_twid_factor, bins backward_twid_factor)
// {
// 	ifft_symmetric_bins* res;
// 	res = convfft_overlap_save_static_n(a, b, forward_twid_factor, backward_twid_factor, sample_length);

// 	for (int i = 0; i < res->length; i++) 
// 	{
// 		printf("%f\n", res->bin[i]);
// 	}

// 	free(res->bin);
// 	free(res);
	
// 	return 0;
// }

int main() 
{
    int length_2, length_1;

    bins twiddle_factors_forward;
    bins twiddle_factors_backward;
    char* string;

    // fft_bin* ifft_bin;

    // length = 8;
    // backward = 1;

    // bins twiddle_factors = precompute_twiddle_factor(length, backward);factors_backwad, fft_size);
    // test_for_overlap_save_static_n(xk_1, xk_2, length_1, twiddle_factors_forward, twiddle_factors_backwad);

    // for (int i = 0; i < length >> 1; i++)
    // {bins twid_factor,
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

    double xk_1[] = { 1, 3, 2, 1, 4, 4, 4, 5, 6, 7, 1, 2, 3, 7, 5, 11};
    length_1 = sizeof(xk_1) / sizeof(xk_1[0]);

    double xk_2[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    length_2 =  sizeof(xk_2) / sizeof(xk_2[0]);

    // test_for_ifft(xk_2, length_2);
    // test_for_symmetric_ifft(xk_1, length_1);
    // test_for_double_real_fft(xk_1, xk_2, length_1, length_2);
    // test_for_conv(xk_1, xk_2, length_1, length_2);
    // test_for_overlap_save(xk_1, xk_2, length_1);

    /* Test for static fft size */
    int fft_size = 16;

    twiddle_factors_forward = precompute_twiddle_factor_radix_4(fft_size, 0);
    twiddle_factors_backward = precompute_twiddle_factor_radix_4(fft_size, 1);

    // for (int i = 0; i < 12; i++)
    // {
    //     string = complex_number_to_string(twiddle_factors_backward[i]);
    //     printf("%s\n", string);
    //     free(string);
    // }

    // test_for_conv_static_n(xk_1, xk_2, length_1, length_2, twiddle_factors_forward, twiddle_factors_backwad, fft_size);
    // test_for_overlap_save_static_n(xk_1, xk_2, length_1, twiddle_factors_forward, twiddle_factors_backwad);

    // test_for_it_radix_4(xk_1, length_1, fft_size, twiddle_factors_forward);
    // test_for_ifft_radix_4_static_n(xk_1, length_1, twiddle_factors_forward, twiddle_factors_backward, fft_size);
    // test_for_symmetric_ifft_radix_4_static_n(xk_1, length_1, twiddle_factors_forward, twiddle_factors_backward, fft_size);
    // test_for_double_real_fft_radix_4_static_n(xk_1, xk_2, length_1, length_2, fft_size, twiddle_factors_forward);

    // test_for_conv(xk_1, xk_2, length_1, length_2);
    // test_for_overlap_save(xk_1, xk_2, length_1);

    // for radix 4
    destroy_bin(twiddle_factors_forward, 3 * (fft_size >> 2));
    destroy_bin(twiddle_factors_backward, 3 * (fft_size >> 2));
    
    return 0;
}
