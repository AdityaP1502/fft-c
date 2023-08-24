/* SIMD Radix 4 Butterfly */

#include "../../header/simd/fft_br4.h"

// 0.707106781
static double _W_8 = 0.707106781;

typedef complex_pair (*cmul_f) (double, double, double, double); 

static void fft_simd_br4_p(complex* c_arr, reg_t n[6], int j0, int j1, int j2, int j3)
{
    // first stage (will be the same for b or f)
    n[0] = complex_load_to_reg(c_arr->real + j1, c_arr->imag + j1);
    n[1] = complex_load_to_reg(c_arr->real + j3, c_arr->imag + j3);
    n[4] = complex_load_to_reg(c_arr->real + j0, c_arr->imag + j0);
    n[5] = complex_load_to_reg(c_arr->real + j2, c_arr->imag + j2);

    n[2] = _mm_complexadd_no_load_pd(n[4], n[5]); //  j0 + jN /2 
    n[3] = _mm_complexsubs_no_load_pd(n[4], n[5]); // j0 - jN / 2
    n[4] = _mm_complexadd_no_load_pd(n[0], n[1]); // jN / 4 + j3N / 4
    n[5] = _mm_complexsubs_no_load_pd(n[0], n[1]); // jN/4 - j3N/4

    // second stage that are the same for b and f
    n[0] = _mm_complexadd_no_load_pd(n[2], n[4]); 
    n[1] = _mm_complexsubs_no_load_pd(n[2], n[4]);
    n[2] = n[3];
    n[3] =  _mm_complex_mulj_no_load_pd(n[5]); // j (jN/4 - j3N/4)
}


void fft_simd_br4_f(complex* c_arr, complex* W, int j0, int k, int ql, int hl)
{
    // W[k] = W_k
    // W[ql + k] = W_2k
    // W[hl + k] = W_3k

    int j1, j2, j3;

    j1 = j0 + ql;
    j2 = j1 + ql;
    j3 = j2 + ql;

    reg_t n1, n2, n3, n4, w;
    reg_t n[6];

    fft_simd_br4_p(c_arr, n, j0, j1, j2, j3);

    n1 = n[0]; //0 + n/2 + n/4 + 3n/4
    n2 = n[1]; // 0 + n/2 - (n/4 + 3n/4)
    n3 = n[2]; // 0 - n/2
    n4 = n[3]; // j (n/4 - 3n/4)
    
    w = complex_load_to_reg(W->real + (ql + k), W->imag + (ql + k));
    n2 = _mm_complexmul_no_load_pd(n2, w); // result 2

    complex_store_reg(c_arr->real + j0, c_arr->imag + j0, n1);
    complex_store_reg(c_arr->real + j2, c_arr->imag + j2, n2);

    n1 = _mm_complexsubs_no_load_pd(n3, n4);
    n2 = _mm_complexadd_no_load_pd(n3, n4);

    w = complex_load_to_reg(W->real + k, W->imag + k);
    n1 = _mm_complexmul_no_load_pd(n1, w);

    w = complex_load_to_reg(W->real + (hl + k), W->imag + (hl + k));
    n2 = _mm_complexmul_no_load_pd(n2, w);

    complex_store_reg(c_arr->real + j1, c_arr->imag + j1, n1);
    complex_store_reg(c_arr->real + j3, c_arr->imag + j3, n2);
}

static void fft_br4_p(complex* c_arr, complex_pair t[4], int j0, int j1, int j2, int j3)
{
    complex_pair t0, t1, t2, t3, t4;

    t0 = complexadd(c_arr->real[j0], c_arr->real[j2], c_arr->imag[j0], c_arr->imag[j2]);
    t1 = complexsub(c_arr->real[j0], c_arr->real[j2], c_arr->imag[j0], c_arr->imag[j2]);
    t2 = complexadd(c_arr->real[j1], c_arr->real[j3], c_arr->imag[j1], c_arr->imag[j3]);
    t3 = complexsub(c_arr->real[j1], c_arr->real[j3], c_arr->imag[j1], c_arr->imag[j3]);

    t4 = t0;
    t0 = complexadd(t0.real, t2.real, t0.imag, t2.imag); 
    t2 = complexsub(t4.real, t2.real, t4.imag, t2.imag);
    t3 = complexmul_j(t3.real, t3.imag);

    t[0] = t0;
    t[1] = t1;
    t[2] = t2;
    t[3] = t3;
}

static complex_pair complexmul_w_w4_f(double re, double dummy, double im, double dummy1)
{
    return complexmul_minj(re, im);
}

static complex_pair complexmul_w_w8_f(double re, double dummy, double im, double dummy1)
{
    complex_pair r;

    r.real = (re + im) * _W_8;
    r.imag = (im - re) * _W_8;

    return r;
}

static complex_pair complexmul_w_3w8_f(double re, double dummy, double im, double dummy1)
{
    complex_pair r;

    r.real = -1 * (re - im) * _W_8;
    r.imag = -1 * (re + im) * _W_8;

    return r;
}

void fft_br4_f_0(complex* c_arr, int j0, int ql, int hl)
{
    complex_pair t0, t1, t2, t3, t4;
    complex_pair a[4];
    int j1, j2, j3;

    j1 = j0 + ql;
    j2 = j1 + ql;
    j3 = j2 + ql;

    // first stage
    fft_br4_p(c_arr, a, j0, j1, j2, j3);
    
    t0 = a[0];
    t1 = a[1];
    t2 = a[2];
    t3 = a[3];

    t4 = t1;
    t1 = complexsub(t1.real, t3.real, t1.imag, t3.imag);
    t3 = complexadd(t4.real, t3.real, t4.imag, t3.imag);

    // save result
    c_arr->real[j0] = t0.real;
    c_arr->real[j1] = t1.real;
    c_arr->real[j2] = t2.real;
    c_arr->real[j3] = t3.real;

    c_arr->imag[j0] = t0.imag;
    c_arr->imag[j1] = t1.imag;
    c_arr->imag[j2] = t2.imag;
    c_arr->imag[j3] = t3.imag;
}

static void fft_br4_f_t(complex* c_arr, complex* W, int j0, int ql, int hl, int k, cmul_f f_wl, cmul_f f_w2l, cmul_f f_w3l)
{
    complex_pair t0, t1, t2, t3, t4;
    complex_pair a[4];
    int j1, j2, j3;

    j1 = j0 + ql;
    j2 = j1 + ql;
    j3 = j2 + ql;

    // first stage
    fft_br4_p(c_arr, a, j0, j1, j2, j3);
    
    t0 = a[0];
    t1 = a[1];
    t2 = a[2];
    t3 = a[3];

    t2 = f_w2l(t2.real, W->real[ql + k], t2.imag, W->imag[ql + k]);

    t4 = t1;
    t1 = complexsub(t1.real, t3.real, t1.imag, t3.imag);
    t1 = f_wl(t1.real, W->real[k], t1.imag, W->imag[k]);

    t3 = complexadd(t4.real, t3.real, t4.imag, t3.imag);
    t3 = f_w3l(t3.real, W->real[hl + k], t3.imag, W->imag[hl + k]);

    // save result
    c_arr->real[j0] = t0.real;
    c_arr->real[j1] = t1.real;
    c_arr->real[j2] = t2.real;
    c_arr->real[j3] = t3.real;

    c_arr->imag[j0] = t0.imag;
    c_arr->imag[j1] = t1.imag;
    c_arr->imag[j2] = t2.imag;
    c_arr->imag[j3] = t3.imag;
}

void fft_br4_f(complex* c_arr, complex* W, int j0, int k, int ql, int hl)
{
    fft_br4_f_t(c_arr, W, j0, ql, hl, k, complexmul, complexmul, complexmul);
}

void fft_br4_f_N_over_8(complex* c_arr, complex* W, int j0, int ql, int hl, int k)
{
    fft_br4_f_t(c_arr, W, j0, ql, hl, k, complexmul_w_w8_f, complexmul_w_w4_f, complexmul_w_3w8_f);
}

void fft_br4_f_N_over_16(complex* c_arr, complex* W, int j0, int ql, int hl, int k)
{
    fft_br4_f_t(c_arr, W, j0, ql, hl, k, complexmul, complexmul_w_w8_f, complexmul);
}

void fft_br4_f_3N_over_16(complex* c_arr, complex* W, int j0, int ql, int hl, int k)
{
    fft_br4_f_t(c_arr, W, j0, ql, hl, k, complexmul, complexmul_w_3w8_f, complexmul);
}

///////////////////////////////// IFFT Radix 4 Butterfly /////////////////////////////////////////////////

static complex_pair complexmul_w_w4_b(double re, double dummy, double im, double dummy1)
{
    return complexmul_j(re, im);
}

static complex_pair complexmul_w_w8_b(double re, double dummy, double im, double dummy1)
{
    complex_pair r;

    r.real = (re - im) * _W_8;
    r.imag = (im + re) * _W_8;

    return r;
}

static complex_pair complexmul_w_3w8_b(double re, double dummy, double im, double dummy1)
{
    complex_pair r;

    r.real = -1 * (re + im) * _W_8;
    r.imag = (re - im) * _W_8;

    return r;
}

void fft_simd_br4_b(complex* c_arr, complex* W, int j0, int k, int ql, int hl)
{
    // W[k] = W_k
    // W[ql + k] = W_2k
    // W[hl + k] = W_3k

    int j1, j2, j3;

    j1 = j0 + ql;
    j2 = j1 + ql;
    j3 = j2 + ql;

    reg_t n1, n2, n3, n4, w;
    reg_t n[6];

    fft_simd_br4_p(c_arr, n, j0, j1, j2, j3);

    n1 = n[0]; //0 + n/2 + n/4 + 3n/4
    n2 = n[1]; // 0 + n/2 - (n/4 + 3n/4)
    n3 = n[2]; // 0 - n/2
    n4 = n[3]; // j (n/4 - 3n/4)
    
    w = complex_load_to_reg(W->real + (ql + k), W->imag + (ql + k));
    n2 = _mm_complexmul_no_load_pd(n2, w); // result 2

    complex_store_reg(c_arr->real + j0, c_arr->imag + j0, n1);
    complex_store_reg(c_arr->real + j2, c_arr->imag + j2, n2);

    n1 = _mm_complexadd_no_load_pd(n3, n4);
    n2 = _mm_complexsubs_no_load_pd(n3, n4);

    w = complex_load_to_reg(W->real + k, W->imag + k);
    n1 = _mm_complexmul_no_load_pd(n1, w);

    w = complex_load_to_reg(W->real + (hl + k), W->imag + (hl + k));
    n2 = _mm_complexmul_no_load_pd(n2, w);

    complex_store_reg(c_arr->real + j1, c_arr->imag + j1, n1);
    complex_store_reg(c_arr->real + j3, c_arr->imag + j3, n2);
}

static void fft_br4_b_t(complex* c_arr, complex* W, int j0, int ql, int hl, int k, cmul_f f_wl, cmul_f f_w2l, cmul_f f_w3l)
{
    complex_pair t0, t1, t2, t3, t4;
    complex_pair a[4];
    int j1, j2, j3;

    j1 = j0 + ql;
    j2 = j1 + ql;
    j3 = j2 + ql;

    // first stage
    fft_br4_p(c_arr, a, j0, j1, j2, j3);
    
    t0 = a[0];
    t1 = a[1];
    t2 = a[2];
    t3 = a[3];

    t2 = f_w2l(t2.real, W->real[ql + k], t2.imag, W->imag[ql + k]);

    t4 = t1;
    t1 = complexadd(t1.real, t3.real, t1.imag, t3.imag);
    t1 = f_wl(t1.real, W->real[k], t1.imag, W->imag[k]);

    t3 = complexsub(t4.real, t3.real, t4.imag, t3.imag);
    t3 = f_w3l(t3.real, W->real[hl + k], t3.imag, W->imag[hl + k]);

    // save result
    c_arr->real[j0] = t0.real;
    c_arr->real[j1] = t1.real;
    c_arr->real[j2] = t2.real;
    c_arr->real[j3] = t3.real;

    c_arr->imag[j0] = t0.imag;
    c_arr->imag[j1] = t1.imag;
    c_arr->imag[j2] = t2.imag;
    c_arr->imag[j3] = t3.imag;
}

void fft_br4_b_0(complex* c_arr, int j0, int ql, int hl)
{
    complex_pair t0, t1, t2, t3, t4;
    complex_pair a[4];
    int j1, j2, j3;

    j1 = j0 + ql;
    j2 = j1 + ql;
    j3 = j2 + ql;

    // first stage
    fft_br4_p(c_arr, a, j0, j1, j2, j3);
    
    t0 = a[0];
    t1 = a[1];
    t2 = a[2];
    t3 = a[3];

    t4 = t1;
    t1 = complexadd(t1.real, t3.real, t1.imag, t3.imag);
    t3 = complexsub(t4.real, t3.real, t4.imag, t3.imag);

    // save result
    c_arr->real[j0] = t0.real;
    c_arr->real[j1] = t1.real;
    c_arr->real[j2] = t2.real;
    c_arr->real[j3] = t3.real;

    c_arr->imag[j0] = t0.imag;
    c_arr->imag[j1] = t1.imag;
    c_arr->imag[j2] = t2.imag;
    c_arr->imag[j3] = t3.imag;
}

void fft_br4_b(complex* c_arr, complex* W, int j0, int k, int ql, int hl)
{
    fft_br4_b_t(c_arr, W, j0, ql, hl, k, complexmul, complexmul, complexmul);
}

void fft_br4_b_N_over_8(complex* c_arr, complex* W, int j0, int ql, int hl, int k)
{
    fft_br4_b_t(c_arr, W, j0, ql, hl, k, complexmul_w_w8_b, complexmul_w_w4_b, complexmul_w_3w8_b);
}

void fft_br4_b_N_over_16(complex* c_arr, complex* W, int j0, int ql, int hl, int k)
{
    fft_br4_b_t(c_arr, W, j0, ql, hl, k, complexmul, complexmul_w_w8_b, complexmul);
}

void fft_br4_b_3N_over_16(complex* c_arr, complex* W, int j0, int ql, int hl, int k)
{
    fft_br4_b_t(c_arr, W, j0, ql, hl, k, complexmul, complexmul_w_3w8_b, complexmul);
}