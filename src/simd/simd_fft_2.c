/* 
    FFT Radix 4 implementation using SIMD. 
    SIMD will be used to boost complex operation. 
*/

#include "../../header/simd/fft_simd.h"

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

static void FFTLIBRARY_CALL bit_reverse_array(complex* c_arr, int length)
{
    // bit reverse the array
    // natural order -> bit reversed order
    // bit reversed order -> natural order

    int j;
    double tempR, tempIm;

    for (int i = 0; i < length; i++)
    {
        j = reverse_bit(i, length);

        if (j < i)
        {
            // already been visited
            continue;
        }

        tempR = c_arr->real[i];
        tempIm = c_arr->imag[i];

        c_arr->real[i] = c_arr->real[j];
        c_arr->imag[i] = c_arr->imag[j];

        c_arr->real[j] = tempR;
        c_arr->imag[j] = tempIm;
    }
}

static void fft_br2(complex* c_arr, complex* W, int i, int j, int k)
{
    // j = i + N' / 2
    reg_t x1, x2, w, t1;

    /// load
    x1 = complex_load_to_reg(c_arr->real + i, c_arr->imag + i);
    x2 = complex_load_to_reg(c_arr->real + j, c_arr->imag + j);
    w = complex_load_to_reg(W->real + k, W->imag + k);

    // xi = xi + xj
    t1 = _mm_complexadd_no_load_pd(x1, x2);
    complex_store_reg(c_arr->real + i, c_arr->imag + i, t1);

    // xj = (x1 - xj) * w1
    t1 = _mm_complexsubs_no_load_pd(x1, x2);
    t1 = _mm_complexmul_no_load_pd(t1, w);
    complex_store_reg(c_arr->real + j, c_arr->imag + j, t1);
}

static void solve_fft_r2(complex* c_arr, complex* w, int fft_size)
{
    int half_length, n_element, n_problems, j0, tR, tIm;

    n_element = fft_size;
    n_problems = 1;
    half_length = n_element / 2;


    while (n_element > 2)
    {
        j0 = 0;

        for (int i = 0; i < n_problems; i++)
        {
            for (int j = 0; j < half_length; j += 2)
            {
                fft_br2(c_arr, w, j0, j0 + half_length, j);
                j0 += 2;
            }

            j0 += half_length;
        }

        // reorder twid factor
        for (int j = 1; j < half_length / 2; j++)
        {
            w->real[j] = w->real[2 * j];
            w->imag[j] = w->imag[2 * j];
        }

        n_element >>= 1;
        n_problems <<= 1;
        half_length >>= 1;
    }

    // solve for n_element == 2
    for (int i = 0; i < fft_size; i += 2)
    {
        tR = c_arr->real[i];
        tIm = c_arr->imag[i];

        c_arr->real[i] = tR + c_arr->real[i + 1];
        c_arr->imag[i] = tIm + c_arr->imag[i + 1];

        c_arr->real[i + 1] = tR - c_arr->real[i + 1];
        c_arr->imag[i + 1] = tIm - c_arr->imag[i + 1];
    }
}

complex* fftr_r2(double* r, complex* w, int arr_size, int fft_size)
{
    complex* c_arr;

    c_arr = complex_arr_create_allign_16(fft_size);
    fill_complex_arr_real(c_arr, r, arr_size);
    solve_fft_r2(c_arr, w, fft_size);
    bit_reverse_array(c_arr, fft_size);
    return c_arr;
}

complex* ifftr_r2(complex* c_arr, complex* w, int arr_size, int fft_size)
{
    complex* c_arr_copy;

    c_arr_copy = complex_arr_create_allign_16(fft_size);
    fill_complex_arr_complex(c_arr_copy, c_arr->real, c_arr->imag, arr_size, arr_size);
    solve_fft_r2(c_arr_copy, w, fft_size);
    bit_reverse_array(c_arr_copy, fft_size);
    normalize_ifft(c_arr, fft_size);
    return c_arr;
}

void ifftr_r2_inplace(complex* c_arr, complex* w, int fft_size)
{
    solve_fft_r2(c_arr, w, fft_size);
    bit_reverse_array(c_arr, fft_size);
    normalize_ifft(c_arr, fft_size);
}

complex** fftdr_r2(double* r1, double* r2, complex* w, int size_1, int size_2, int fft_size)
{
    complex* c_arr1;
    complex* c_arr2;
    complex** res;

    res = malloc( 2 * sizeof(complex*) );
    c_arr1 = complex_arr_create_allign_16(fft_size);
    c_arr2 = complex_arr_create_allign_16(fft_size);

    res[0] = c_arr1;
    res[1] = c_arr2;

    fill_complex_arr_complex(c_arr1, r1, r2, size_1, size_2);
    solve_fft_r2(c_arr1, w, fft_size);
    bit_reverse_array(c_arr1, fft_size);
    seperate_combined_output(c_arr1, c_arr1, c_arr2, fft_size);

    return res;
}   