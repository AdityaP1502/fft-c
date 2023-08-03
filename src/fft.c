/*
    This module provide basic operation
    and definition that will be used
    to implement fft in recursive and iterative
*/

#include "../header/fft.h"

void FFTLIBRARY_CALL destroy_bin(bins bin, int length)
{
  for (int i = 0; i < length; i++)
  {
    free(bin[i]);
  }

  free(bin);
}

bins FFTLIBRARY_CALL deep_copy_bins_complex(bins xn, int length, int pad_length)
{
  // padded xn with zero
  bins padded_xn = malloc((length + pad_length) * sizeof(complex_number *));

  for (int i = 0; i < length; i++)
  {
    padded_xn[i] = create_complex_number(xn[i]->real, xn[i]->imag);
  }

  // create a 0 complex number
  for (int i = length; i < length + pad_length; i++)
  {
    padded_xn[i] = create_complex_number(0, 0);
  }

  return padded_xn;
}

bins FFTLIBRARY_CALL copy_bins_complex(bins xn, int length, int pad_length)
{
  // padded xn with zero
  bins padded_xn = malloc((length + pad_length) * sizeof(complex_number *));

  for (int i = 0; i < length; i++)
  {
    padded_xn[i] = xn[i];
  }

  // create a 0 complex number
  for (int i = length; i < length + pad_length; i++)
  {
    padded_xn[i] = create_complex_number(0, 0);
  }

  return padded_xn;
}

double *FFTLIBRARY_CALL copy_bins_real(double *xn, int length, int pad_length)
{
  double *padded_xn = calloc(length + pad_length, sizeof(double));

  for (int i = 0; i < length; i++)
  {
    padded_xn[i] = xn[i];
  }

  return padded_xn;
}

void FFTLIBRARY_CALL clear_pad(bins bin, int length, int pad_length)
{
  for (int i = length; i < length + pad_length; i++)
  {
    free(bin[i]);
  }
}

double *FFTLIBRARY_CALL clear_pad_real(double *xn, int length)
{
  double *new_xn = malloc(length * sizeof(double));

  for (int i = 0; i < length; i++)
  {
    new_xn[i] = xn[i];
  }

  return new_xn;
}

int FFTLIBRARY_CALL nearest_power_of_2(int length)
{
  // return the difference between length and the closees biggest power of 2
  // 0 if length is a power of 2

  // init -> 2 ^ 0
  int power_of_2 = 1;

  while (length > power_of_2)
  {
    power_of_2 = power_of_2 << 1;
  }

  return power_of_2;
}

int FFTLIBRARY_CALL nearest_power_of_4(int length)
{
  // return the difference between length and the closeest biggest power of 4
  // 0 if length is a power of 4

  int power_of_4 = 1;

  while (length > power_of_4)
  {
    power_of_4 <<= 2;
  }

  return power_of_4;
}

bins FFTLIBRARY_CALL convert_real_to_complex(double *xn, int length)
{
  bins complexs = malloc(length * sizeof(complex_number *));

  for (int i = 0; i < length; i++)
  {
    complexs[i] = create_complex_number(xn[i], 0);
  }

  return complexs;
}

double *FFTLIBRARY_CALL convert_complex_to_real(bins xn, int length)
{
  double *real = malloc(length * sizeof(double));

  // Assumed that bin doesn't contain any imaginary part
  for (int i = 0; i < length; i++)
  {
    real[i] = xn[i]->real;
  }

  return real;
}

bins FFTLIBRARY_CALL precompute_twiddle_factor(int length, int backward)
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

  double real_part;
  double imag_part;

  int half_length;
  int N_over_eight;
  int N_over_four;

  bins twiddle_factors;

  theta = backward ? (M_PI * 2) / length : -(M_PI * 2) / length;
  magic_factor_1 = 1 - 2 * pow(sin(theta / 2), 2); // C
  magic_factor_2 = sin(theta);                     // S

  half_length = length >> 1;
  N_over_four = half_length >> 1;
  N_over_eight = N_over_four >> 1;

  twiddle_factors = malloc(half_length * sizeof(complex_number *));

  // calculate the first N / 8 points

  // init the first factor -> r * theta = 0
  twiddle_factors[0] = create_complex_number(1, 0);
  switch (length)
  {
  case 2:
    break;

  case 4:
    twiddle_factors[1] = backward ? create_complex_number(0, 1) : create_complex_number(0, -1);
    break;

  default:
    for (int i = 1; i <= N_over_eight; i++)
    {
      // using singleteon
      // cos ((r + 1)theta) = C * cos(rtheta) - S * sin(rtheta)
      // sin((r + 1)theta) = S * cos(rtheta) + C * sin(rtheta)

      real_part = magic_factor_1 * twiddle_factors[i - 1]->real - magic_factor_2 * twiddle_factors[i - 1]->imag;
      imag_part = magic_factor_2 * twiddle_factors[i - 1]->real + magic_factor_1 * twiddle_factors[i - 1]->imag;

      twiddle_factors[i] = create_complex_number(real_part, imag_part);
    }

    for (int i = 1; i < N_over_eight; i++)
    {
      // cos((N/8 + k) * theta) = sin((N/8 - k) * theta)
      // sin((N/8 + k) * theta) = cos((N/8 - k) * theta)

      real_part = backward ? twiddle_factors[N_over_eight - i]->imag : -twiddle_factors[N_over_eight - i]->imag;
      imag_part = backward ? twiddle_factors[N_over_eight - i]->real : -twiddle_factors[N_over_eight - i]->real;

      twiddle_factors[N_over_eight + i] = create_complex_number(real_part, imag_part);
    }

    // calculate the other N / 4 points in second quadrand
    for (int i = 0; i < N_over_four; i++)
    {
      // cos((N/4 + k) * theta) = -sin(k * theta)
      // sin((N/4 + k) * theta) = cos(k * theta)

      real_part = backward ? -twiddle_factors[i]->imag : twiddle_factors[i]->imag;
      imag_part = backward ? twiddle_factors[i]->real : -twiddle_factors[i]->real;

      twiddle_factors[N_over_four + i] = create_complex_number(real_part, imag_part);
    }
  }

  return twiddle_factors;
}

bins FFTLIBRARY_CALL precompute_twiddle_factor_radix_4(int length, int backward)
{
  /*
    Twiddle factor is symmetric
    therefore only need to compute the first quadrant
    because after theta=45, the value of sin and cos if flip
    only need to calculate the first half of the first quadrant

    Because of high multiplier in radix 4 fft, we need 3N/4 point
  */

  double theta;

  double magic_factor_1;
  double magic_factor_2;

  double real_part;
  double imag_part;

  int half_length;
  int N_over_eight;
  int N_over_four;

  bins twiddle_factors;

  theta = backward ? (M_PI * 2) / length : -(M_PI * 2) / length;
  magic_factor_1 = 1 - 2 * pow(sin(theta / 2), 2); // C
  magic_factor_2 = sin(theta);                     // S

  half_length = length >> 1;
  N_over_four = half_length >> 1;
  N_over_eight = N_over_four >> 1;

  twiddle_factors = malloc((3 * N_over_four) * sizeof(complex_number *));

  // calculate the first N / 8 points

  // init the first factor -> r * theta = 0
  twiddle_factors[0] = create_complex_number(1, 0);
  switch (length)
  {
  case 2:
    break;

  case 4:
    twiddle_factors[1] = backward ? create_complex_number(0, 1) : create_complex_number(0, -1);
    break;

  default:
    for (int i = 1; i <= N_over_eight; i++)
    {
      // using singleteon
      // cos ((r + 1)theta) = C * cos(rtheta) - S * sin(rtheta)
      // sin((r + 1)theta) = S * cos(rtheta) + C * sin(rtheta)

      real_part = magic_factor_1 * twiddle_factors[i - 1]->real - magic_factor_2 * twiddle_factors[i - 1]->imag;
      imag_part = magic_factor_2 * twiddle_factors[i - 1]->real + magic_factor_1 * twiddle_factors[i - 1]->imag;

      twiddle_factors[i] = create_complex_number(real_part, imag_part);
    }

    for (int i = 1; i < N_over_eight; i++)
    {
      // cos((N/8 + k) * theta) = sin((N/8 - k) * theta)
      // sin((N/8 + k) * theta) = cos((N/8 - k) * theta)

      real_part = backward ? twiddle_factors[N_over_eight - i]->imag : -1 * twiddle_factors[N_over_eight - i]->imag;
      imag_part = backward ? twiddle_factors[N_over_eight - i]->real : -1 * twiddle_factors[N_over_eight - i]->real;

      twiddle_factors[N_over_eight + i] = create_complex_number(real_part, imag_part);
    }

    // calculate the other N / 4 points in second quadrand
    for (int i = 0; i <= N_over_four; i++)
    {
      // cos((N/4 + k) * theta) = -sin(k * theta)
      // sin((N/4 + k) * theta) = cos(k * theta)

      real_part = backward ? -twiddle_factors[i]->imag : twiddle_factors[i]->imag;
      imag_part = backward ? twiddle_factors[i]->real : -1 * twiddle_factors[i]->real;

      twiddle_factors[N_over_four + i] = create_complex_number(real_part, imag_part);
    }

    // Extra N/4 
    // calculate the other N / 4 points in third quadrand
    for (int i = 1; i < N_over_four; i++)
    {
      // cos((N/2 + k) * theta) = -cos((N/2 - k) * theta)
      // sin((N/2 + k) * theta) = -sin((N/2 - k) * theta)

      real_part = twiddle_factors[half_length - i]->real;
      imag_part = -1 * twiddle_factors[half_length - i]->imag;

      twiddle_factors[half_length + i] = create_complex_number(real_part, imag_part);
    }
  }

  return twiddle_factors;
}

bins FFTLIBRARY_CALL combined_two_real_input(double *xn_1, double *xn_2, int max_length, int min_length)
{
    bins combined_bin = malloc(max_length * sizeof(complex_number *));

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

void FFTLIBRARY_CALL seperate_combined_output(bins combined_output, int length, fft_bins **dest)
{
    fft_bins *result_1;
    fft_bins *result_2;

    bins fr;
    bins gr;

    double real, imag;

    fr = malloc(length * sizeof(complex_number *));
    gr = malloc(length * sizeof(complex_number *));

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

void FFTLIBRARY_CALL normalize_ifft(bins xk, int length)
{
    for (int i = 0; i < length; i++)
    {
        xk[i]->real = xk[i]->real / length;
        xk[i]->imag = xk[i]->imag / length;
    }
}