#ifndef COMPLEX_H
#define COMPLEX_H

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "dll_export_api.h"

typedef struct complex_ {
  double real;
  double imag;
} complex_number;

// use  5 digit after decimal seperator

// maximum digit/char given to a number when convert to string
#define MAXIMUM_DIGIT_CHAR 25

// maximum character in a complex string
#define MAXIMUM_TO_STRING_BUFFER 55 

// create a complex number
FFTLIBRARY_API complex_number* FFTLIBRARY_CALL create_complex_number(double real, double imag);

// create a complex number from angle(radians)
FFTLIBRARY_API complex_number* FFTLIBRARY_CALL create_complex_number_from_angle(double angle);

// to string
FFTLIBRARY_API char* FFTLIBRARY_CALL complex_number_to_string(complex_number* complex_number);

// calculate complex exponential
FFTLIBRARY_API complex_number* FFTLIBRARY_CALL complex_exponent(double angle);

// find the magnitude of complex number
FFTLIBRARY_API double FFTLIBRARY_CALL complex_magnitude(complex_number* a);

// find the angle of complex number in radian
FFTLIBRARY_API double FFTLIBRARY_CALL complex_angle(complex_number* a);

/*  Add two complex number, if dst is not null, will return a new number, 
*   else would overwrite a
*/
FFTLIBRARY_API void FFTLIBRARY_CALL complex_add(complex_number* dst, complex_number* a, complex_number* b);

/*  multiply two complex number, if dst is not null, will return a new number, 
*   else would overwrite a
*/
FFTLIBRARY_API void FFTLIBRARY_CALL complex_multiply(complex_number* dst, complex_number* a, complex_number* b);

/*  substract two complex number, if dst is not null, will return a new number, 
*   else would overwrite a
*/
FFTLIBRARY_API void FFTLIBRARY_CALL complex_substract(complex_number* dst, complex_number* a, complex_number* b);

/*  divide two complex number, if dst is not null, will return a new number, 
*   else would overwrite a
*/
FFTLIBRARY_API void FFTLIBRARY_CALL complex_divide(complex_number* dst, complex_number* a, complex_number* b);

#endif