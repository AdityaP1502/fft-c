#include "../header/iterative_fft.h"
#include "../header/fft_radix_4.h"
#include "../header/conv.h"

/*
    This is a module that implement
    convolutional 1D using Iterative FFT
*/

// each input already has the linear length N = (L + M - 1)
// and sample length is power of 2
ifft_symmetric_bins *FFTLIBRARY_CALL convfft_overlap_save(double *a, double *b, int sample_length)
{


    int out_length;

    fft_bins *placeholder;
    fft_bins **double_res;

    ifft_symmetric_bins *conv_out;

    double *a_padded;
    double *b_padded;
    double *temp;

    a_padded = copy_bins_real(a, sample_length, 0);
    b_padded = copy_bins_real(b, sample_length, 0);

    double_res = fft_double_real(a_padded, b_padded, sample_length, sample_length, "NR");
    out_length = double_res[0]->length;

    placeholder = double_res[0];

    for (int i = 0; i < out_length; i++)
    {
        complex_multiply(NULL, placeholder->fft_bins[i], double_res[1]->fft_bins[i]);
    }

    destroy_bin(double_res[1]->fft_bins, out_length);
    free(double_res[1]);
    free(double_res);

    conv_out = ifft_iterative_symmetric(placeholder->fft_bins, out_length, "NR");
    temp = conv_out->bin;
    conv_out->bin = clear_pad_real(temp, sample_length);
    conv_out->length = sample_length;

    destroy_bin(placeholder->fft_bins, placeholder->length);
    free(placeholder);

    free(a_padded);
    free(b_padded);
    free(temp);

    return conv_out;
}

ifft_symmetric_bins *FFTLIBRARY_CALL convfft(double *a, double *b, int length_a, int length_b)
{
    // conv using fft == circular convolution
    // linear convolution == circular convolution when N = N1 + N2 - 1

    int conv_length, pad_length_a, pad_length_b, out_length;

    fft_bins *placeholder;
    fft_bins **double_res;

    ifft_symmetric_bins *conv_out;

    double *a_padded;
    double *b_padded;
    double *temp;

    conv_length = length_a + length_b - 1;

    pad_length_a = conv_length - length_a;
    pad_length_b = conv_length - length_b;

    a_padded = copy_bins_real(a, length_a, pad_length_a);
    b_padded = copy_bins_real(b, length_b, pad_length_b);

    double_res = fft_double_real(a_padded, b_padded, conv_length, conv_length, "NR");
    out_length = double_res[0]->length;

    placeholder = double_res[0];

    for (int i = 0; i < out_length; i++)
    {
        complex_multiply(NULL, placeholder->fft_bins[i], double_res[1]->fft_bins[i]);
    }

    destroy_bin(double_res[1]->fft_bins, out_length);
    free(double_res[1]);
    free(double_res);

    // printf("%d\n", out_length);

    conv_out = ifft_iterative_symmetric(placeholder->fft_bins, out_length, "NR");
    temp = conv_out->bin;
    conv_out->bin = clear_pad_real(temp, conv_length);
    conv_out->length = conv_length;

    destroy_bin(placeholder->fft_bins, placeholder->length);
    free(placeholder);

    free(a_padded);
    free(b_padded);
    free(temp);

    return conv_out;
}

ifft_symmetric_bins *FFTLIBRARY_CALL convfft_overlap_save_static_n(double *a, double *b, bins forward_twid_factor, bins backward_twid_factor, int fft_size)
{
    // each input already has the linear length N = (L + M - 1)
    // and sample length is power of 2

    int out_length;

    fft_bins *placeholder;
    fft_bins **double_res;

    ifft_symmetric_bins *conv_out;

    double *a_padded;
    double *b_padded;
    double *temp;

    a_padded = copy_bins_real(a, fft_size, 0);
    b_padded = copy_bins_real(b, fft_size, 0);

    double_res = fft_double_real_static_n(forward_twid_factor, a_padded, b_padded, fft_size, fft_size, fft_size, "NR");
    out_length = double_res[0]->length;

    placeholder = double_res[0];

    for (int i = 0; i < out_length; i++)
    {
        complex_multiply(NULL, placeholder->fft_bins[i], double_res[1]->fft_bins[i]);
    }

    destroy_bin(double_res[1]->fft_bins, out_length);
    free(double_res[1]);
    free(double_res);

    conv_out = ifft_iterative_symmetric_static_n(placeholder->fft_bins, backward_twid_factor, fft_size, out_length, "NR");
    temp = conv_out->bin;
    conv_out->bin = clear_pad_real(temp, fft_size);
    conv_out->length = fft_size;

    destroy_bin(placeholder->fft_bins, placeholder->length);
    free(placeholder);

    free(a_padded);
    free(b_padded);
    free(temp);

    return conv_out;
}

ifft_symmetric_bins *FFTLIBRARY_CALL convfft_static_n(double *a, double *b, int length_a, int length_b, bins forward_twid_factor, bins backward_twid_factor, int fft_size)
{
    // conv using fft == circular convolution
    // linear convolution == circular convolution when N = N1 + N2 - 1

    int conv_length, pad_length_a, pad_length_b, out_length;

    fft_bins *placeholder;
    fft_bins **double_res;

    ifft_symmetric_bins *conv_out;

    double *a_padded;
    double *b_padded;
    double *temp;

    conv_length = length_a + length_b - 1;

    pad_length_a = conv_length - length_a;
    pad_length_b = conv_length - length_b;

    a_padded = copy_bins_real(a, length_a, pad_length_a);
    b_padded = copy_bins_real(b, length_b, pad_length_b);

    double_res = fft_double_real_static_n(forward_twid_factor, a_padded, b_padded, conv_length, conv_length, fft_size, "NR");
    out_length = double_res[0]->length;

    placeholder = double_res[0];

    for (int i = 0; i < out_length; i++)
    {
        complex_multiply(NULL, placeholder->fft_bins[i], double_res[1]->fft_bins[i]);
    }

    destroy_bin(double_res[1]->fft_bins, out_length);
    free(double_res[1]);
    free(double_res);

    conv_out = ifft_iterative_symmetric_static_n(placeholder->fft_bins, backward_twid_factor, fft_size, out_length, "NR");
    temp = conv_out->bin;
    conv_out->bin = clear_pad_real(temp, conv_length);
    conv_out->length = conv_length;

    destroy_bin(placeholder->fft_bins, placeholder->length);
    free(placeholder);

    free(a_padded);
    free(b_padded);
    free(temp);

    return conv_out;
}

ifft_symmetric_bins *FFTLIBRARY_CALL convfft4_overlap_save_static_n(double *a, double *b, bins forward_twid_factor, bins backward_twid_factor, int fft_size)
{
    // each input already has the linear length N = (L + M - 1)
    // and sample length is power of 2

    int out_length;

    fft_bins *placeholder;
    fft_bins **double_res;

    ifft_symmetric_bins *conv_out;

    double *a_padded;
    double *b_padded;
    double *temp;

    a_padded = copy_bins_real(a, fft_size, 0);
    b_padded = copy_bins_real(b, fft_size, 0);

    double_res = fft_radix_4_double_real_static_n(forward_twid_factor, a_padded, b_padded, fft_size, fft_size, fft_size);
    out_length = double_res[0]->length;

    placeholder = double_res[0];

    for (int i = 0; i < out_length; i++)
    {
        complex_multiply(NULL, placeholder->fft_bins[i], double_res[1]->fft_bins[i]);
    }

    destroy_bin(double_res[1]->fft_bins, out_length);
    free(double_res[1]);
    free(double_res);

    conv_out = ifft_radix_4_iterative_symmetric_static_n(backward_twid_factor, fft_size, placeholder->fft_bins, out_length);
    temp = conv_out->bin;
    conv_out->bin = clear_pad_real(temp, fft_size);
    conv_out->length = fft_size;

    destroy_bin(placeholder->fft_bins, placeholder->length);
    free(placeholder);

    free(a_padded);
    free(b_padded);
    free(temp);

    return conv_out;
}

ifft_symmetric_bins *FFTLIBRARY_CALL convfft4_static_n(double *a, double *b, int length_a, int length_b, bins forward_twid_factor, bins backward_twid_factor, int fft_size)
{
    // conv using fft == circular convolution
    // linear convolution == circular convolution when N = N1 + N2 - 1

    int conv_length, pad_length_a, pad_length_b, out_length;

    fft_bins *placeholder;
    fft_bins **double_res;

    ifft_symmetric_bins *conv_out;

    double *a_padded;
    double *b_padded;
    double *temp;

    conv_length = length_a + length_b - 1;

    pad_length_a = conv_length - length_a;
    pad_length_b = conv_length - length_b;

    a_padded = copy_bins_real(a, length_a, pad_length_a);
    b_padded = copy_bins_real(b, length_b, pad_length_b);

    double_res = fft_radix_4_double_real_static_n(forward_twid_factor, a_padded, b_padded, conv_length, conv_length, fft_size);
    out_length = double_res[0]->length;

    placeholder = double_res[0];

    for (int i = 0; i < out_length; i++)
    {
        complex_multiply(NULL, placeholder->fft_bins[i], double_res[1]->fft_bins[i]);
    }

    destroy_bin(double_res[1]->fft_bins, out_length);
    free(double_res[1]);
    free(double_res);

    conv_out = ifft_radix_4_iterative_symmetric_static_n(backward_twid_factor, fft_size, placeholder->fft_bins, out_length);
    temp = conv_out->bin;
    conv_out->bin = clear_pad_real(temp, conv_length);
    conv_out->length = conv_length;

    destroy_bin(placeholder->fft_bins, placeholder->length);
    free(placeholder);

    free(a_padded);
    free(b_padded);
    free(temp);

    return conv_out;
}

