#include "../header/complex.h"

static char COMPLEX_NUMBER_FORMAT_PRECISION[] = "%.5f";

complex_number *create_complex_number(double real, double imag)
{
  complex_number *complex = malloc(sizeof(complex_number));
  complex->real = real;
  complex->imag = imag;

  return complex;
}

complex_number *create_complex_number_from_angle(double angle)
{
  complex_number *complex = malloc(sizeof(complex_number));
  complex->real = cos(angle);
  complex->imag = sin(angle);

  return complex;
}

static void remove_trailing_zeros_in_decimal(char *number)
{
  if (!number)
    return;

  char *is_decimal;
  int length, zero_is_trail_zero, i;

  length = strlen(number);

  if (length)
  {
    zero_is_trail_zero = number[length - 1] == '0';
    i = length - 1;

    while (zero_is_trail_zero)
    {
      // will guarenteed to break before i == 0
      // because will break after number[i] == '.'
      number[i] = '\0';
      i--;
      zero_is_trail_zero = number[i] == '0';
    }

    if (number[i] == '.')
    {
      number[i] = '\0';
    }
  }
}

static void format_complex_number(char *buffer, char *real_number, char *imag_number, int is_imag_negative)
{
  int j;

  remove_trailing_zeros_in_decimal(real_number);
  remove_trailing_zeros_in_decimal(imag_number);

  if (real_number && imag_number)
  {
    if (is_imag_negative)
    {
      j = snprintf(buffer, MAXIMUM_TO_STRING_BUFFER, "%s-j%s", real_number, imag_number);
    }
    else
    {
      j = snprintf(buffer, MAXIMUM_TO_STRING_BUFFER, "%s+j%s", real_number, imag_number);
    }
  }
  else if (real_number && !imag_number)
  {
    j = snprintf(buffer, MAXIMUM_TO_STRING_BUFFER, "%s", real_number);
  }
  else if (!real_number && imag_number)
  {
    if (is_imag_negative)
    {
      j = snprintf(buffer, MAXIMUM_TO_STRING_BUFFER, "-j%s", imag_number);
    }
    else
    {
      j = snprintf(buffer, MAXIMUM_TO_STRING_BUFFER, "j%s", imag_number);
    }
  }
  else
  {
    j = snprintf(buffer, MAXIMUM_TO_STRING_BUFFER, "%s", "0");
  }
}

static void complex_number_negative_formatting(complex_number *complex_number, char *real, char *imag)
{
  int j;

  // negative formatting for real number
  if (complex_number->real > -0.00001 && complex_number->real <= 0)
  {
    // would be conveted to 0
    // omit the sign
    j = snprintf(real, MAXIMUM_DIGIT_CHAR, COMPLEX_NUMBER_FORMAT_PRECISION, fabs(complex_number->real));
  }
  else
  {
    j = snprintf(real, MAXIMUM_DIGIT_CHAR, COMPLEX_NUMBER_FORMAT_PRECISION, complex_number->real);
  }

  // imaginary number placed before the j, so the number sign must be omitted
  j = snprintf(imag, MAXIMUM_DIGIT_CHAR, COMPLEX_NUMBER_FORMAT_PRECISION, fabs(complex_number->imag));
}

char *complex_number_to_string(complex_number *complex_number)
{
  int is_imag_negative;
  char *buffer;
  char *real_number;
  char *imag_number;

  buffer = malloc(MAXIMUM_TO_STRING_BUFFER * sizeof(char));
  real_number = malloc(MAXIMUM_DIGIT_CHAR * sizeof(char));
  imag_number = malloc(MAXIMUM_DIGIT_CHAR * sizeof(char));

  // based on snprintf design
  // all string will be null terminated

  complex_number_negative_formatting(complex_number, real_number, imag_number);

  is_imag_negative = complex_number->imag <= -0.00001;

  if (!strncmp(real_number, "0.0000", 6))
  {
    free(real_number);
    real_number = NULL;
  }

  if (!strncmp(imag_number, "0.0000", 6))
  {
    free(imag_number);
    imag_number = NULL;
  }

  format_complex_number(buffer, real_number, imag_number, is_imag_negative);

  // free non null resources
  if (real_number)
  {
    free(real_number);
  }

  if (imag_number)
  {
    free(imag_number);
  }

  return buffer;
}

double complex_magnitude(complex_number *a)
{
  return sqrt(a->real * a->real + a->imag * a->imag);
}

double complex_angle(complex_number *a)
{
  if (a->imag == 0)
  {
    return 0;
  }

  if (a->real == 0)
  {
    return M_PI / 2;
  }

  return atan(a->imag / a->real);
}

void complex_add(complex_number *dst, complex_number *a, complex_number *b)
{
  if (dst)
  {
    dst->real = a->real + b->real;
    dst->imag = a->imag + b->imag;
  }
  else
  {
    a->real = a->real + b->real;
    a->imag = a->imag + b->imag;
  }
}

void complex_substract(complex_number *dst, complex_number *a, complex_number *b)
{
  if (dst)
  {
    dst->real = a->real - b->real;
    dst->imag = a->imag - b->imag;
  }
  else
  {
    a->real = a->real - b->real;
    a->imag = a->imag - b->imag;
  }
}

void complex_multiply(complex_number *dst, complex_number *a, complex_number *b)
{
  // (a + bi) * (c + di) = (ac - bd) + (ad + bd)i
  if (dst)
  {
    dst->real = a->real * b->real - a->imag * b->imag;
    dst->imag = a->real * b->imag + b->real * a->imag;
  }
  else
  {
    double t;
    t = a->real;
    a->real = a->real * b->real - a->imag * b->imag;
    a->imag = t * b->imag + b->real * a->imag;
  }
}

void complex_divide(complex_number *dst, complex_number *a, complex_number *b)
{
  // (a + bi) / (c + di) = (a + bi) *  (c - di) / (c^2 + d^2)
  // (a + bi) *  (c - di) / (c^2 + d^2) = (ac + bd) + (bc - ad)i / mag_squared(b)
  double mag_squared = b->real * b->real + b->imag * b->imag;
  if (dst)
  {
    dst->real = (a->real * b->real + a->imag * b->imag) / mag_squared;
    dst->imag = (-a->real * b->imag + b->real * a->imag) / mag_squared;
  }
  else
  {
    double t;
    t = a->real;
    a->real = (a->real * b->real + a->imag * b->imag) / mag_squared;
    a->imag = (-t * b->imag + b->real * a->imag) / mag_squared;
  }
}

complex_number *complex_exponent(double angle)
{
  // calculate e^(jx)
  // using euler identity
  // e^(jx) = cos x + j * sin(x)

  complex_number *res = malloc(sizeof(complex_number));
  res->real = cos(angle);
  res->imag = sin(angle);

  return res;
}
