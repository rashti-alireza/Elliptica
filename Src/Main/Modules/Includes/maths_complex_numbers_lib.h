#ifndef maths_complex_numbers_LIB_H
#define maths_complex_numbers_LIB_H


#include <complex.h>
#include <tgmath.h>//for the type generate macros.

double complex Ylm(const unsigned l, int m, const double theta, const double phi);
double complex dYlm_dphi(const unsigned l, const int m, const double theta, const double phi);
double complex dYlm_dtheta(const unsigned l, const int m, const double theta, const double phi);
void *alloc_double_complex(const unsigned N);
#endif




