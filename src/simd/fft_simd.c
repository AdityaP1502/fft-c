#include "../../header/simd/complex_simd.h"


complex* FFTLIBRARY_CALL precompute_twiddle_factor(int length, int backward)
{
  /*
    Twiddle factor is symmetric
    therefore only need to compute the first quadrant
    because after theta=45, the value of sin and cos if flip
    only need to calculate the first half of the first quadrant

    due to the nature of fft, we only need N/2 point of the twiddle factor
  */

  double theta;

  double magic_factor_1;
  double magic_factor_2;

  int half_length;
  int N_over_eight;
  int N_over_four;

  complex* twiddle_factors;

  theta = backward ? (M_PI * 2) / length : -(M_PI * 2) / length;
  magic_factor_1 = 1 - 2 * pow(sin(theta / 2), 2); // C
  magic_factor_2 = sin(theta);                     // S

  half_length = length >> 1;
  N_over_four = half_length >> 1;
  N_over_eight = N_over_four >> 1;

  twiddle_factors = complex_arr_create_allign_16(half_length);

  twiddle_factors->real[0] = 1;
  twiddle_factors->imag[0] = 0;

  switch (length)
  {
  case 2:
    break;

  case 4:
    // twiddle_factors[1] = backward ? create_complex_number(0, 1) : create_complex_number(0, -1);
    twiddle_factors->real[1] = 0;
    twiddle_factors->imag[1] = backward ? 1 : -1;
    break;

  default:
    for (int i = 1; i <= N_over_eight; i++)
    {
      // using singleteon
      // cos ((r + 1)theta) = C * cos(rtheta) - S * sin(rtheta)
      // sin((r + 1)theta) = S * cos(rtheta) + C * sin(rtheta)

      twiddle_factors->real[i] = magic_factor_1 * twiddle_factors->real[i - 1] - magic_factor_2 * twiddle_factors->imag[i - 1];
      twiddle_factors->imag[i] = magic_factor_2 * twiddle_factors->real[i - 1] + magic_factor_1 * twiddle_factors->imag[i - 1];
    }

    for (int i = 1; i < N_over_eight; i++)
    {
      // cos((N/8 + k) * theta) = sin((N/8 - k) * theta)
      // sin((N/8 + k) * theta) = cos((N/8 - k) * theta)

      twiddle_factors->real[N_over_eight + i] = backward ? twiddle_factors->imag[N_over_eight - i] : -twiddle_factors->imag[N_over_eight - i];
      twiddle_factors->imag[N_over_eight + i] = backward ? twiddle_factors->real[N_over_eight - i] : -twiddle_factors->real[N_over_eight - i];
    }

    // calculate the other N / 4 points in second quadrand
    for (int i = 0; i < N_over_four; i++)
    {
      // cos((N/4 + k) * theta) = -sin(k * theta)
      // sin((N/4 + k) * theta) = cos(k * theta)

      twiddle_factors->real[N_over_four + i] = backward ? -twiddle_factors->imag[i] : twiddle_factors->imag[i];
      twiddle_factors->imag[N_over_four + i] = backward ? twiddle_factors->real[i] : -twiddle_factors->real[i];
    }
  }
  return twiddle_factors;
}

complex* FFTLIBRARY_CALL precompute_twiddle_factor_radix_4(int length, int backward)
{
  /*
    Twiddle factor is symmetric
    therefore only need to compute the first quadrant
    because after theta=45, the value of sin and cos if flip
    only need to calculate the first half of the first quadrant

    due to the nature of fft, we only need N/2 point of the twiddle factor
  */

  double theta;

  double magic_factor_1;
  double magic_factor_2;

  int half_length;
  int N_over_eight;
  int N_over_four;

  complex* twiddle_factors;

  theta = backward ? (M_PI * 2) / length : -(M_PI * 2) / length;
  magic_factor_1 = 1 - 2 * pow(sin(theta / 2), 2); // C
  magic_factor_2 = sin(theta);                     // S

  half_length = length >> 1;
  N_over_four = half_length >> 1;
  N_over_eight = N_over_four >> 1;

  twiddle_factors = complex_arr_create_allign_16(3 * N_over_four);

  twiddle_factors->real[0] = 1;
  twiddle_factors->imag[0] = 0;

  switch (length)
  {
  case 2:
    break;

  case 4:
    // twiddle_factors[1] = backward ? create_complex_number(0, 1) : create_complex_number(0, -1);
    twiddle_factors->real[1] = 0;
    twiddle_factors->imag[1] = backward ? 1 : -1;
    break;

  default:
    for (int i = 1; i <= N_over_eight; i++)
    {
      // using singleteon
      // cos ((r + 1)theta) = C * cos(rtheta) - S * sin(rtheta)
      // sin((r + 1)theta) = S * cos(rtheta) + C * sin(rtheta)

      twiddle_factors->real[i] = magic_factor_1 * twiddle_factors->real[i - 1] - magic_factor_2 * twiddle_factors->imag[i - 1];
      twiddle_factors->imag[i] = magic_factor_2 * twiddle_factors->real[i - 1] + magic_factor_1 * twiddle_factors->imag[i - 1];
    }

    for (int i = 1; i < N_over_eight; i++)
    {
      // cos((N/8 + k) * theta) = sin((N/8 - k) * theta)
      // sin((N/8 + k) * theta) = cos((N/8 - k) * theta)

      twiddle_factors->real[N_over_eight + i] = backward ? twiddle_factors->imag[N_over_eight - i] : -twiddle_factors->imag[N_over_eight - i];
      twiddle_factors->imag[N_over_eight + i] = backward ? twiddle_factors->real[N_over_eight - i] : -twiddle_factors->real[N_over_eight - i];
    }

    // calculate the other N / 4 points in second quadrand
    for (int i = 0; i <= N_over_four; i++)
    {
      // cos((N/4 + k) * theta) = -sin(k * theta)
      // sin((N/4 + k) * theta) = cos(k * theta)

      twiddle_factors->real[N_over_four + i] = backward ? -twiddle_factors->imag[i] : twiddle_factors->imag[i];
      twiddle_factors->imag[N_over_four + i] = backward ? twiddle_factors->real[i] : -twiddle_factors->real[i];
    }

    // Extra N/4 
    // calculate the other N / 4 points in third quadrand
    for (int i = 1; i < N_over_four; i++)
    {
      // cos((N/2 + k) * theta) = -cos((N/2 - k) * theta)
      // sin((N/2 + k) * theta) = -sin((N/2 - k) * theta)

      twiddle_factors->real[half_length + i] = twiddle_factors->real[half_length - i];
      twiddle_factors->imag[half_length + i] = -1 * twiddle_factors->imag[half_length - i];
    }
  }
  return twiddle_factors;
}