#include "../header/fft_radix_4.h"

/*
    This module provide implemnetation
    on the iterative dif fft algorithm in radix 4
    The implementation only provide NR mode
*/

static int FFTLIBRARY_CALL quaternary_reversed(int num, int length)
{
    // reversed the quaternary digit
    int q_reversed, count;

    count = length;
    q_reversed = num;

    // remove the first digit from num
    num >>= 2;
    count >>= 2;

    while (count)
    {
        q_reversed <<= 2;
        q_reversed |= num & 3;
        num >>= 2;
        count >>= 2;
    }

    return q_reversed & (length);
}

static void FFTLIBRARY_CALL quaternary_reversed_array(bins xk, int length)
{
    // in place reversed array
    complex_number* temp;
    int reversed_i;

    for (int i = 0; i < length; i++)
    {
        reversed_i = quaternary_reversed(i, length);

        if (reversed_i < i) continue;

        temp = xk[reversed_i];

        xk[reversed_i] = xk[i];
        xk[i] = temp;
    }
}

static void FFTLIBRARY_CALL radix_4_butterfly(bins xk, int j1, int j2, int j3, int j4, bins wl, bins temp)
{
    // j1 = l, j2 = l + N/4, j3 = l + N/2, j4 = l + 3N/4
    complex_number *a1 = temp[0];
    complex_number *a2 = temp[1];
    complex_number *a3 = temp[2];
    complex_number *a4 = temp[3];

    complex_add(a1, xk[j1], xk[j3]);       // xl + x(l + N/2)
    complex_add(a2, xk[j2], xk[j4]);       // x(l + N/4) + x[l + 3N/4]
    complex_substract(a3, xk[j1], xk[j3]); // xl - x(l + N/2)
    complex_substract(a4, xk[j2], xk[j4]); // x(l + N/4) - x[l + 3N/4]

    // Yl calculation
    complex_add(xk[j1], a1, a2); // yl = a1 + a2

    // Gl calculation
    complex_substract(xk[j3], a1, a2);
    complex_multiply(NULL, xk[j3], wl[1]); // gl = (a1 - a2) * w2l

    // Zl calculation
    complex_multiply_with_j(NULL, a4, 1); // j * ( x(l + N/4) - x(l + 3N/4) )
    complex_substract(xk[j2], a3, a4);
    complex_multiply(NULL, xk[j2], wl[0]); // zl = (a3 - ja4) * wl

    // Hl calculation
    complex_add(xk[j4], a3, a4);
    complex_multiply(NULL, xk[j4], wl[2]); // hl = (a3 + ja4) * w3l
}

static void FFTLIBRARY_CALL do_in_place_fft_radix_4_nr(bins xk, int length, bins twiddle_factors)
{
    int n_problems, n_element, Jstart, half_length, quarter_length, three_quarter_length;
    int n1, n2, n3; // twiddle index
    int j1, j2, j3, j4; // xk index
    complex_number* t; // temp (to store twid factor that surpass boundary)

    bins wl;
    bins temp;

    wl = malloc(sizeof(complex_number *) * 3);
    temp = malloc(sizeof(complex_number *) * 4);

    temp[0] = create_complex_number(0, 0);
    temp[1] = create_complex_number(0, 0);
    temp[2] = create_complex_number(0, 0);
    temp[3] = create_complex_number(0, 0);

    n_problems = 1;
    n_element = length;

    while (n_element > 1)
    {
        half_length = n_element / 2;
        quarter_length = half_length / 2;
        three_quarter_length = 3 * quarter_length;

        j1 = 0;
        
        for (int i = 0; i < n_problems; i++)
        {
            j2 = j1 + quarter_length;
            j3 = j1 + half_length;
            j4 = j1 + three_quarter_length;

            for (int j = 0; j < quarter_length; j++)
            {
                n1 = n_problems * j;
                n2 = n1 * 2;
                n3 = n2 + n1;

                wl[0] = twiddle_factors[n1];
                wl[1] = twiddle_factors[n2];
                wl[2] = twiddle_factors[n3];

                radix_4_butterfly(xk, j1, j2, j3, j4, wl, temp);
                j1++;
                j2++;
                j3++;
                j4++;
            }

            j1 = j4;
        }

        n_problems <<= 2;
        n_element >>= 2;
    }
}

fft_bins* fft_radix_4_iterative_static_n(bins forward_twid_factor, int fft_size, double* xn, int length)
{
    int pad_length, input_length;
    bins input_array, output_array, temp, twiddle_factors;

    fft_bins *res = malloc(sizeof(fft_bins));

    // convert double to complex
    input_array = convert_real_to_complex(xn, length);

    // pad the input if not power of 2
    pad_length = nearest_power_of_4(length) - length;

    temp = input_array;
    input_array = copy_bins_complex(input_array, length, pad_length);
    free(temp);

    input_length = length + pad_length;

    if (fft_size != input_length) {
        // TODO: Terminate when the sample length isn't the same with the input length
        // call error function
    } 

    do_in_place_fft_radix_4_nr(input_array, input_length, forward_twid_factor);

    quaternary_reversed_array(input_array, input_length - 1);

    res->fft_bins = input_array;
    res->length = input_length;

    return res->fft_bins ? res : NULL;
}
