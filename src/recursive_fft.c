#include "../header/recursive_fft.h"

/*
  This module provide the implmenentation of
  recursive dit fft algorithm 
*/

static bins do_fft_dit_recursive(bins xn, int length, int forward)
{
  // xn is zero paddded
  // length is guarenteed a power of 2
  // fft bin hold references to complex_number struct pointer

  int c, half_length;
  double angle;
  double dt;

  complex_number *a;
  complex_number *b;
  complex_number *w;

  bins ye;
  bins yo;
  bins even_xn;
  bins odd_xn;
  bins res;

  res = malloc(sizeof(complex_number *) * length);

  // init res
  for (int i = 0; i < length; i++)
  {
    res[i] = malloc(sizeof(complex_number));
  }

  // base case
  if (length == 1)
  {
    res[0]->real = xn[0]->real;
    res[0]->imag = xn[0]->imag;

    return res;
  }

  half_length = length >> 1;
  angle = 0;
  dt = forward ? -(2 * M_PI) / (length) : (2 * M_PI) / (length);

  // decomposition mem alloc
  even_xn = malloc(half_length * sizeof(complex_number *));
  odd_xn = malloc(half_length * sizeof(complex_number *));

  // fill even and odd coeeficient
  c = 0;
  for (int i = 0; i < half_length; i++)
  {
    even_xn[i] = xn[c];
    odd_xn[i] = xn[c + 1];
    c += 2;
  }

  // do fft for even and odd term
  ye = do_fft_dit_recursive(even_xn, half_length, forward);
  yo = do_fft_dit_recursive(odd_xn, half_length, forward);

  // combine
  for (int i = 0; i < half_length; i++)
  {
    a = ye[i];
    b = yo[i];
    w = complex_exponent(angle);  // calculate w^i
    complex_multiply(NULL, w, b); // calculate w^i * yo[i]

    // res[i] = ye[i] + w^i * yo[i]
    // res[i + half_length] = ye[i] - w^i * yo[i]
    complex_add(res[i], a, w);
    complex_substract(res[i + half_length], a, w);
    free(w);     // complex_exponent malloc a memory
    angle += dt; // increment angle
  }

  // free intermediate value content
  destroy_bin(ye, half_length);
  destroy_bin(yo, half_length);

  // just free the holder, don't need to destroy the content
  free(even_xn);
  free(odd_xn);

  return res;
}

fft_bins* fft_recursive(double *xn, int length)
{
  fft_bins* result;
  bins xn_complex;
  bins padded_xn;
  bins bin;
  int pad_length;
  double dt; 
  
  result = malloc(sizeof(fft_bins));

  // convert xn to complex
  xn_complex = convert_real_to_complex(xn, length);

  // check if length is power of 2
  pad_length = nearest_power_of_2(length) - length;
  

  if (pad_length)
  {
    // paddded xn if not power of 2
    padded_xn = deep_copy_bins_complex(xn_complex, length, pad_length);
    bin = do_fft_dit_recursive(padded_xn, length + pad_length, 1);
    destroy_bin(padded_xn, length + pad_length);
  }
  else 
  {
    bin = do_fft_dit_recursive(xn_complex, length, 1);
    destroy_bin(xn_complex, length);
  }
  
  result->fft_bins = bin;
  result->length = length + pad_length;

  return result;
}

fft_bins* ifft_recursive(bins Xk, int length)
{
  fft_bins* result;
  bins padded_Xk;
  bins bin;
  int pad_length;

  result = malloc(sizeof(fft_bins));

  // check if length is power of 2
  pad_length = nearest_power_of_2(length) - length;

  if (pad_length)
  {
    // paddded xn if not power of 2
    padded_Xk = deep_copy_bins_complex(Xk, length, pad_length);
    bin = do_fft_dit_recursive(padded_Xk, length + pad_length, 0);
    clear_pad(padded_Xk, length, pad_length);
  }
  else 
  {
    bin = do_fft_dit_recursive(Xk, length, 0);
  }
  
  result->fft_bins = bin;
  result->length = length + pad_length;

  // normalize result
  for (int i = 0; i < result->length; i++)
  {
    result->fft_bins[i]->real /= result->length;
    result->fft_bins[i]->imag /= result->length;
  }

  return result;
}

ifft_symmetric_bins* ifft_recursive_symmetric(bins Xk, int length) 
{
  ifft_symmetric_bins* real_bin = malloc(sizeof(ifft_symmetric_bins));

  fft_bins* ifft_result = ifft_recursive(Xk, length);
  double* real_result = convert_complex_to_real(ifft_result->fft_bins, ifft_result->length);

  real_bin->bin = real_result;
  real_bin->length = ifft_result->length;
  
  destroy_bin(ifft_result->fft_bins, ifft_result->length);
  free(ifft_result);

  return real_bin;
}