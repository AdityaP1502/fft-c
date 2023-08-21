/* 
    FFT Radix 4 implementation using SIMD. 
    SIMD will be used to boost complex operation. 
*/
#include "../../header/simd/fft_simd.h"
#include "../../header/simd/fft_br4.h"

static void subproblem_routine_optimized_f(complex* c_arr, complex* w, int ql, int hl, int j0)
{  
    int eighthl = ql / 2;
    int sixtenthl = eighthl / 2;

    int iter = eighthl - 1;
    int start = j0;

    start++; // start at j1

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < iter; j++)
        {
            // do butterfly with double complex
            fft_simd_br4_f(c_arr, w, start, j, ql, hl);
            start += 2;
        }

        // normal butterfly
        fft_br4_f(c_arr, w, start, 2 * iter, ql, hl);
        start += 2;
    }

    // deal with special cases
    fft_br4_0(c_arr, j0, ql, hl);
    fft_br4_f_N_over_16(c_arr, j0 + sixtenthl, ql, hl);
    fft_br4_f_N_over_8(c_arr, w, j0 + eighthl, ql, hl, 2);
    fft_br4_f_3N_over_16(c_arr, w, j0 + (3 * sixtenthl), ql, hl, 3);
}

static void solve_N_4_fft_r4(complex* c_arr, complex* w, int n_problems)
{
    int j0 = 0;
    int ql = 1;
    int hl = 2;

    for (int i = 0; i < n_problems; i++)
    {
        fft_br4_0(c_arr, j0, ql, hl);
    }
}

static void solve_N_16_fft_r4(complex* c_arr, complex* w, int n_problems)
{
    int j0 = 0;
    int ql = 4;
    int hl = 8;
    
    for (int i = 0; i < n_problems; i++)
    {
        fft_br4_0(c_arr, j0, ql, hl);
        fft_br4_f_N_over_16(c_arr, j0, ql, hl);
        fft_br4_f_N_over_8(c_arr, w, j0, ql, hl, 2);
        fft_br4_f_3N_over_16(c_arr, w, j0, ql, hl, 3);
        j0 += 16;
    }
}

static void solve_fft_r4(complex* c_arr, complex* w, int fft_size)
{
    // W memory structure
    // w[0....N/4 - 1] = wl
    // w[N/4....N/2-1] = w2l
    // w[N/2...3N/4-1] = w3l
    
    int l, hl, ql, j0;
    int n_elements, n_problems;

    l = fft_size;
    hl = l / 2;
    ql = l / 4;

    n_elements = l;
    n_problems = 1;

    while (n_elements > 16)
    {
        j0 = 0;
        for (int i = 0; i < n_problems; i++)
        {
            subproblem_routine_optimized_f(c_arr, w, ql, hl, j0);
            j0 += n_elements;
        }
    }

    if (n_elements == 16)
    {
        solve_N_16_fft_r4(c_arr, w, n_problems);
        n_problems *= 4;
    }

    solve_N_4_fft_r4(c_arr, w, n_problems);
}