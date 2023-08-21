#include <stdio.h>
#include "../header/simd/complex_simd.h"

void test_complex_add()
{
    double x_re[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double x_im[] = {-1, 2.13, 3.14, 2, 2, 1, 7, 1};

    // (1 -j + 3 + j3.14) and (2 + j2.13) + (4 + j2)
    reg_t res = _mm_complexadd_pd(x_re, x_im, x_re + 2, x_im + 2);

    _mm_storeu_pd(x_re, res.re);
    _mm_storeu_pd(x_im, res.im);

    // expected (4 + j2.14) and (6 + j4.13)
    printf("%f+j%f\n", x_re[0], x_im[0]);
    printf("%f+j%f\n", x_re[1], x_im[1]);

    // received 4 + j2.140000
    // received 6 + j4.130000

    // PASSED
}

void test_complex_subs()
{
    double x_re[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double x_im[] = {-1, 2.13, 3.14, 2, 2, 1, 7, 1};

    // (1 -j - 3 + j3.14) and (2 + j2.13) - (4 + j2)
    reg_t res = _mm_complexsubs_pd(x_re, x_im, x_re + 2, x_im + 2);

    _mm_storeu_pd(x_re, res.re);
    _mm_storeu_pd(x_im, res.im);

    // expected (-2 - j4.14) and (-2 + j0.13)
    printf("%f+j%f\n", x_re[0], x_im[0]);
    printf("%f+j%f\n", x_re[1], x_im[1]);

    // received -2 - j4.140000
    // received -2 + j0.130000

    // PASSED
}

void test_complex_mul()
{
    double x_re[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double x_im[] = {-1, 2.13, 3.14, 2, 2, 1, 7, 1};

    // (1 -j) * (3 + j3.14) and (2 + j2.13) * (4 + j2)
    reg_t res = _mm_complexmul_pd(x_re, x_im, x_re + 2, x_im + 2);

    _mm_storeu_pd(x_re, res.re);
    _mm_storeu_pd(x_im, res.im);

    // expected (6.14 + j0.14) and (4 + j12.52)
    printf("%f+j%f\n", x_re[0], x_im[0]);
    printf("%f+j%f\n", x_re[1], x_im[1]);

    // received 6.14 + j0.140000
    // received 3.74 + j12.52000

    // PASSED
}

void test_complex_div()
{
    double x_re[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double x_im[] = {-1, 2.13, 3.14, 2, 2, 1, 7, 1};

    // (1 -j) / (3 + j3.14) and (2 + j2.13)  / (4 + j2)
    reg_t res = _mm_complexdiv_pd(x_re, x_im, x_re + 2, x_im + 2);

    _mm_storeu_pd(x_re, res.re);
    _mm_storeu_pd(x_im, res.im);

    // expected (0.0074 - j0.32) and (0.613 + j0.226)
    printf("%f+j%f\n", x_re[0], x_im[0]);
    printf("%f+j%f\n", x_re[1], x_im[1]);

    // received -0.007423 -j0.325564
    // received 0.613000 + j0.226000

    // PASSED
}

int main()
{
    test_complex_add();
    test_complex_subs();
    test_complex_mul();
    test_complex_div();
    return 0;
}