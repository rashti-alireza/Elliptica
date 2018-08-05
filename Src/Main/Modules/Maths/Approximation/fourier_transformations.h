#include "core_lib.h"
#include "maths_general_lib.h"
#include <fftw3.h>

void fftw_1d_ChebyshevExtrema_coeffs(double *const values,double *const coeffs,const unsigned n);
void fftw_1d_ChebyshevExtrema_values(double *const values,double *const coeffs,const unsigned n);
void fftw_3d_ChebyshevExtrema_coeffs(double *const values,double *const coeffs,const unsigned *const n);
void fftw_3d_ChebyshevExtrema_values(double *const values,double *const coeffs,const unsigned *const n);
void rft_1d_ChebyshevExtrema_coeffs(double *const values ,double *const coeffs,const unsigned n);
