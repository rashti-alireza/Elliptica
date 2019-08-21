#include "core_lib.h"
#include "maths_general_lib.h"
#include <complex.h>

void rft_1d_ChebyshevExtrema_coeffs(double *const values ,double *const coeffs,const unsigned n);
void rft_1d_ChebyshevNodes_coeffs(double *const values ,double *const coeffs,const unsigned n);
void r2cft_1d_EquiSpaced_coeffs(const double *const value,double complex *const coeffs,const unsigned n,const unsigned l);
