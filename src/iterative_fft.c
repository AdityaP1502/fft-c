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

static int max(int a, int b)
{
    return a > b ? a : b;
}

static int reverse_bit(int number, int length) 
{
    // reverse the bit of the number
    // where length_of_bit == log(length)
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

    return reversed_num;
}

static void bit_reverse_array(bins dest, bins input, int length)
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

    complex_number* temp;

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

static void gentleman_sandle_butterfly(bins xk, int j1, int j2, complex_number* twid_factor)
{
    complex_number* temp;
    complex_number* placeholder;

    double real, imag;

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

static void nr_subproblems(bins xk, bins twiddle_factors, int half_length, int offset, int N_problems)
{
    complex_number* twid_factor;

    for (int j = 0; j < half_length; j++) 
    { 
        twid_factor = twiddle_factors[N_problems * j];
        gentleman_sandle_butterfly(xk, j + offset, j + offset + half_length, twid_factor);
    }
}

static void do_in_place_fft_nr(bins xk, int length, bins twiddle_factors)
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

static void rn_subproblems(bins xk, bins twiddle_factors, int start, int distance, int length) 
{
    int j, m;
    complex_number* temp;
    complex_number* twid_factor;

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

static void do_in_place_fft_rn(bins xk, int length, bins twiddle_factors)
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

static bins do_fft_iterative(bins input_array, int length, int backward, char* order)
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

fft_bins* fft_iterative(double* xn, int length, char* order)
{
    int pad_length, input_length;
    bins input_array, output_array, temp, twiddle_factors;

    fft_bins* res = malloc(sizeof(fft_bins));

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

static void normalize_ifft(bins xk, int length)
{
    for (int i = 0; i < length; i++)
    {
        xk[i]->real = xk[i]->real / length;
        xk[i]->imag = xk[i]->imag / length;
    }
}

fft_bins* ifft_iterative(bins xk, int length, char* order)
{
    int pad_length, input_length;
    char* str;
    bins input_array, output_array, temp, twiddle_factors;

    // printf("%d\n", length);
    // for (int i = 0; i < length; i++)
    // {
    //     str = complex_number_to_string(xk[i]);
    //     printf("%s\n", str);
    //     free(str);
    // }

    temp = NULL;
    output_array = NULL;
    input_array = xk;

    fft_bins* res = malloc(sizeof(fft_bins));

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

ifft_symmetric_bins* ifft_iterative_symmetric(bins xk, int length, char* order)
{
    fft_bins* output_ifft;
    ifft_symmetric_bins* res;
    double* real_ifft;

    res = malloc(sizeof(ifft_symmetric_bins));

    output_ifft = ifft_iterative(xk, length, order);
    real_ifft = convert_complex_to_real(output_ifft->fft_bins, output_ifft->length);

    res->bin = real_ifft;
    res->length = output_ifft->length;

    destroy_bin(output_ifft->fft_bins, output_ifft->length);
    free(output_ifft);

    return res;
}

static bins combined_two_real_input(double* xn_1, double* xn_2, int max_length, int min_length)
{
    bins combined_bin = malloc(max_length * sizeof(complex_number*));

    // make X[k] = xn_1[k] + jxn_2[k]
    for (int i = 0; i < min_length; i++)
    {
        combined_bin[i] = create_complex_number(xn_1[i], xn_2[i]);
    }

    for (int i = min_length; i < max_length; i++)
    {
        combined_bin[i] = create_complex_number(xn_1[i], 0);
    }

    return combined_bin;
}

static void seperate_combined_output(bins combined_output, int length, fft_bins** dest)
{
    fft_bins* result_1;
    fft_bins* result_2;

    bins fr;
    bins gr;

    double real, imag;

    fr = malloc(length * sizeof(complex_number*));
    gr = malloc(length * sizeof(complex_number*));

    result_1 = malloc(sizeof(fft_bins));
    result_2 = malloc(sizeof(fft_bins));

    fr[0] = create_complex_number(combined_output[0]->real, 0);
    gr[0] = create_complex_number(combined_output[0]->imag, 0);

    for (int i = 1; i < length; i++)
    {
        // fr = 0.5 (xr + x*(N - r))
        real = 0.5 * (combined_output[i]->real + combined_output[length - i]->real);
        imag = 0.5 * (combined_output[i]->imag - combined_output[length - i]->imag);
        
        fr[i] = create_complex_number(real, imag);

        // gr = 0.5j (x*(N - r) - x(r) -> 0.5(jA* - jB)
        // gr = 0.5 (jA.real + A.imag - (-B.imag + jB.real)) -> 1/2(A.imag + B.imag) + 1/2j(A.real - B.real)

        real = 0.5 * (combined_output[length - i]->imag + combined_output[i]->imag);
        imag = 0.5 * (combined_output[length - i]->real - combined_output[i]->real);

        gr[i] = create_complex_number(real, imag);
    }   

    result_1->fft_bins = fr;
    result_1->length = length;

    result_2->fft_bins = gr;
    result_2->length = length;

    dest[0] = result_1;
    dest[1] = result_2;
}

fft_bins** fft_double_real(double* xn_1, double* xn_2, int length_1, int length_2, char* order)
{
    int input_length, pad_length;

    bins combined_bin;
    bins combined_result_bin;
    bins temp;

    fft_bins** res;
    fft_bins* placeholder;


    temp = NULL;
    
    // length_1 must be bigger than length 2
    if (length_2 > length_1) 
    {
        res = fft_double_real(xn_2, xn_1, length_2, length_1, order);
        placeholder = res[0];
        res[0] = res[1];
        res[1] = placeholder;
        return res;
    }

    res = malloc(2 * sizeof(fft_bins*));

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


