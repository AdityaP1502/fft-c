#include <string.h>
#include "../header/iterative_fft.h"

/*
    This module provide implemnetation
    on the iterative dif fft algorithm
    for the RN order, NR order, and the NN order
*/

static char RN[] = "RN";
static char NR[] = "NR";
static char NN[] = "NN";

// static int FFTLIBRARY_CALL max(int a, int b)
// {
//     return a > b ? a : b;
// }

static int FFTLIBRARY_CALL reverse_bit(int number, int length)
{
    // reverse the bit of the number
    // length is a power of 2

    int reversed_num;

    reversed_num = 0;
    length >>= 1;

    while (length)
    {
        reversed_num <<= 1;
        reversed_num |= number & 1;
        number >>= 1;
        length >>= 1;
    }

    return reversed_num & (length - 1);
}

static void FFTLIBRARY_CALL bit_reverse_array(bins dest, bins input, int length)
{
    // bit reverse the array
    // natural order -> bit reversed order
    // bit reversed order -> natural order

    int j;

    if (dest)
    {
        for (int i = 0; i < length; i++)
        {
            j = reverse_bit(i, length);
            dest[i] = input[j];
        }

        return;
    }

    complex_number *temp;

    for (int i = 0; i < length; i++)
    {
        j = reverse_bit(i, length);

        if (j < i)
        {
            // already been visited
            continue;
        }

        temp = input[i];
        input[i] = input[j];
        input[j] = temp;
    }
}

static void FFTLIBRARY_CALL gentleman_sandle_butterfly(bins xk, int j1, int j2, complex_number *twid_factor)
{
    complex_number *temp;
    complex_number *placeholder;

    // double real, imag;

    temp = xk[j1];
    placeholder = create_complex_number(0, 0);

    // xk[j] = xk[j] + xk[j + N/2]
    complex_add(placeholder, temp, xk[j2]);
    xk[j1] = placeholder;

    // xk[j + N / 2] = (xk[j] - xk[j + N / 2]) * W;
    complex_substract(NULL, temp, xk[j2]);
    complex_multiply(NULL, temp, twid_factor);

    placeholder = xk[j2];

    xk[j2] = temp;

    free(placeholder);
}

static void FFTLIBRARY_CALL nr_subproblems(bins xk, bins twiddle_factors, int half_length, int offset, int N_problems)
{
    complex_number *twid_factor;

    for (int j = 0; j < half_length; j++)
    {
        twid_factor = twiddle_factors[N_problems * j];
        gentleman_sandle_butterfly(xk, j + offset, j + offset + half_length, twid_factor);
    }
}

static void FFTLIBRARY_CALL do_in_place_fft_nr(bins xk, int length, bins twiddle_factors)
{
    // do in-place fft
    int N_problems, N_element, half_length, offset;

    N_element = length;
    N_problems = 1;

    while (N_element > 1)
    {
        offset = 0;
        half_length = N_element >> 1;

        for (int i = 0; i < N_problems; i++)
        {
            nr_subproblems(xk, twiddle_factors, half_length, offset, N_problems);
            offset += N_element;
        }

        // each stage the elements are halfed
        // whereas the problems are doubled
        N_element >>= 1;
        N_problems <<= 1;
    }
}

static void FFTLIBRARY_CALL rn_subproblems(bins xk, bins twiddle_factors, int start, int distance, int length)
{
    int j, m;
    // complex_number *temp;
    complex_number *twid_factor;

    j = start;
    m = 0;

    while (j < length - 1)
    {
        twid_factor = twiddle_factors[m];
        gentleman_sandle_butterfly(xk, j, j + distance, twid_factor);
        m++;
        j += 2 * distance;
    }
}

static void FFTLIBRARY_CALL do_in_place_fft_rn(bins xk, int length, bins twiddle_factors)
{
    // do in place fft using rn order
    int N_element, N_problems;

    N_element = length;
    N_problems = 1;

    while (N_element > 1)
    {
        for (int i = 0; i < N_problems; i++)
        {
            rn_subproblems(xk, twiddle_factors, i, N_problems, length);
        }

        N_element >>= 1;
        N_problems <<= 1;
    }
}

static bins FFTLIBRARY_CALL do_fft_iterative(bins input_array, int length, int backward, char *order)
{
    bins output_array;
    bins twiddle_factors;

    output_array = NULL;
    twiddle_factors = precompute_twiddle_factor(length, backward);

    if (!strncmp(order, RN, 2))
    {
        // order == RN
        // reversed the order of the input array (in-place)
        bit_reverse_array(NULL, input_array, length);
        bit_reverse_array(NULL, twiddle_factors, length >> 1);
        do_in_place_fft_rn(input_array, length, twiddle_factors);
        output_array = input_array;
    }
    else if (!strncmp(order, NR, 2))
    {
        // order == NR
        do_in_place_fft_nr(input_array, length, twiddle_factors);
        bit_reverse_array(NULL, input_array, length);
        output_array = input_array;
    }
    else if (!strncmp(order, NN, 2))
    {
        // order == NN
    }

    destroy_bin(twiddle_factors, length >> 1);

    return output_array;
}

fft_bins *FFTLIBRARY_CALL fft_iterative(double *xn, int length, char *order)
{
    int pad_length, input_length;
    bins input_array, output_array, temp;

    fft_bins *res = malloc(sizeof(fft_bins));

    // convert double to complex
    input_array = convert_real_to_complex(xn, length);

    // pad the input if not power of 2
    pad_length = nearest_power_of_2(length) - length;

    temp = input_array;
    input_array = copy_bins_complex(input_array, length, pad_length);
    free(temp);

    input_length = length + pad_length;

    output_array = do_fft_iterative(input_array, input_length, 0, order);

    res->fft_bins = output_array;
    res->length = input_length;

    return res->fft_bins ? res : NULL;
}

fft_bins *FFTLIBRARY_CALL ifft_iterative(bins xk, int length, char *order)
{
    int pad_length, input_length;
    bins input_array, output_array;

    // printf("%d\n", length);
    // for (int i = 0; i < length; i++)
    // {
    //     str = complex_number_to_string(xk[i]);
    //     printf("%s\n", str);
    //     free(str);
    // }

    output_array = NULL;
    input_array = xk;

    fft_bins *res = malloc(sizeof(fft_bins));

    // pad the input if not power of 2
    pad_length = nearest_power_of_2(length) - length;
    input_array = deep_copy_bins_complex(input_array, length, pad_length);
    input_length = length + pad_length;

    output_array = do_fft_iterative(input_array, input_length, 1, order);

    if (output_array)
    {
        normalize_ifft(output_array, length);
    }

    res->fft_bins = output_array;
    res->length = input_length;

    return res->fft_bins ? res : NULL;
}

ifft_symmetric_bins *FFTLIBRARY_CALL ifft_iterative_symmetric(bins xk, int length, char *order)
{
    fft_bins *output_ifft;
    ifft_symmetric_bins *res;
    double *real_ifft;

    res = malloc(sizeof(ifft_symmetric_bins));

    output_ifft = ifft_iterative(xk, length, order);
    real_ifft = convert_complex_to_real(output_ifft->fft_bins, output_ifft->length);

    res->bin = real_ifft;
    res->length = output_ifft->length;

    destroy_bin(output_ifft->fft_bins, output_ifft->length);
    free(output_ifft);

    return res;
}

fft_bins **FFTLIBRARY_CALL fft_double_real(double *xn_1, double *xn_2, int length_1, int length_2, char *order)
{
    int input_length, pad_length;

    bins combined_bin;
    bins combined_result_bin;
    bins temp;

    fft_bins **res;
    fft_bins *placeholder;

    if (length_2 > length_1)
    {
        res = fft_double_real(xn_2, xn_1, length_2, length_1, order);
        placeholder = res[0];
        res[0] = res[1];
        res[1] = placeholder;
        return res;
    }

    res = malloc(2 * sizeof(fft_bins *));

    combined_bin = combined_two_real_input(xn_1, xn_2, length_1, length_2);
    input_length = length_1;

    // pad the input if not power of 2
    pad_length = nearest_power_of_2(input_length) - input_length;

    temp = combined_bin;
    combined_bin = copy_bins_complex(combined_bin, input_length, pad_length);
    free(temp);

    input_length += pad_length;

    combined_result_bin = do_fft_iterative(combined_bin, input_length, 0, order);

    seperate_combined_output(combined_result_bin, input_length, res);

    destroy_bin(combined_result_bin, input_length);

    return res;
}

// LENGTH MUST EQUAL WITH INPUT LENGTH. 
 // CHECK MUST BE DONE BEFORE CALLING THIS FUNCTION TO ENSURE PROPER EXECUTION
static bins FFTLIBRARY_CALL do_fft_iterative_static_n(bins input_array, bins twiddle_factors, int length, int backward, char *order)
{
  

    bins output_array;

    output_array = NULL;

    if (!strncmp(order, RN, 2))
    {
        // order == RN
        // reversed the order of the input array (in-place)
        bit_reverse_array(NULL, input_array, length);
        bit_reverse_array(NULL, twiddle_factors, length >> 1);
        do_in_place_fft_rn(input_array, length, twiddle_factors);
        output_array = input_array;
    }
    else if (!strncmp(order, NR, 2))
    {
        // order == NR
        do_in_place_fft_nr(input_array, length, twiddle_factors);
        bit_reverse_array(NULL, input_array, length);
        output_array = input_array;
    }
    else if (!strncmp(order, NN, 2))
    {
        // order == NN
    }

    return output_array;
}


// Sample length is the fixed size determined by the user
// twid factor length is determined by the sample length
fft_bins *FFTLIBRARY_CALL fft_iterative_static_n(bins forward_twid_factor, int fft_size, double *xn, int length, char *order)
{

    int pad_length, input_length;
    bins input_array, output_array, temp;

    fft_bins *res = malloc(sizeof(fft_bins));

    // convert double to complex
    input_array = convert_real_to_complex(xn, length);

    // pad the input if not power of 2
    pad_length = nearest_power_of_2(length) - length;

    temp = input_array;
    input_array = copy_bins_complex(input_array, length, pad_length);
    free(temp);

    input_length = length + pad_length;

    if (fft_size != input_length) {
        // TODO: Terminate when the sample length isn't the same with the input length
        // call error function
    } 

    output_array = do_fft_iterative_static_n(input_array, forward_twid_factor, input_length, 0, order);

    res->fft_bins = output_array;
    res->length = input_length;

    return res->fft_bins ? res : NULL;
}

// Sample length is the fixed size determined by the user
// twid factor length is determined by the sample length
fft_bins *FFTLIBRARY_CALL ifft_iterative_static_n(bins xk, bins backward_twid_factor, int fft_size, int length, char *order)
{
    int pad_length, input_length;
    bins input_array, output_array;

    output_array = NULL;
    input_array = xk;

    fft_bins *res = malloc(sizeof(fft_bins));

    // pad the input if not power of 2
    pad_length = nearest_power_of_2(length) - length;
    input_array = deep_copy_bins_complex(input_array, length, pad_length);
    input_length = length + pad_length;
    
    if (fft_size != input_length) {
        // TODO: Terminate when the sample length isn't the same with the input length
        // call error function
    } 

    output_array = do_fft_iterative_static_n(input_array, backward_twid_factor, input_length, 1, order);

    if (output_array)
    {
        normalize_ifft(output_array, length);
    }

    res->fft_bins = output_array;
    res->length = input_length;

    return res->fft_bins ? res : NULL;
}

ifft_symmetric_bins *FFTLIBRARY_CALL ifft_iterative_symmetric_static_n(bins xk, bins backward_twid_factor, int fft_size, int length, char *order)
{
    fft_bins *output_ifft;
    ifft_symmetric_bins *res;
    double *real_ifft;

    res = malloc(sizeof(ifft_symmetric_bins));

    output_ifft = ifft_iterative_static_n(xk, backward_twid_factor, fft_size, length, order);
    real_ifft = convert_complex_to_real(output_ifft->fft_bins, output_ifft->length);

    res->bin = real_ifft;
    res->length = output_ifft->length;

    destroy_bin(output_ifft->fft_bins, output_ifft->length);
    free(output_ifft);

    return res;
}

fft_bins **FFTLIBRARY_CALL fft_double_real_static_n(bins forward_twid_factor, double *xn_1, double *xn_2, int length_1, int length_2, int fft_size, char *order)
{
    int input_length, pad_length;

    bins combined_bin;
    bins combined_result_bin;
    bins temp;

    fft_bins **res;
    fft_bins *placeholder;

    temp = NULL;

    // length_1 must be bigger than length 2
    if (length_2 > length_1)
    {
        res = fft_double_real_static_n(forward_twid_factor, xn_2, xn_1, length_2, length_1, fft_size, order);
        placeholder = res[0];
        res[0] = res[1];
        res[1] = placeholder;
        return res;
    }

    res = malloc(2 * sizeof(fft_bins *));

    combined_bin = combined_two_real_input(xn_1, xn_2, length_1, length_2);
    input_length = length_1;

    // pad the input if not power of 2
    pad_length = nearest_power_of_2(input_length) - input_length;

    temp = combined_bin;
    combined_bin = copy_bins_complex(combined_bin, input_length, pad_length);
    free(temp);

    input_length += pad_length; // sample length must be the same as the length of twid factor

    if (fft_size != input_length) {
        // TODO: Terminate when the sample length isn't the same with the input length
        // call error function
        printf("%d %d\n", fft_size, input_length);
        printf("Mismatch length\n");
    } 

    combined_result_bin = do_fft_iterative_static_n(combined_bin, forward_twid_factor, input_length, 0, order);

    seperate_combined_output(combined_result_bin, input_length, res);

    destroy_bin(combined_result_bin, input_length);

    return res;
}

