/*
// Alireza Rashti
// July 2018
*/

#include "fftw.h"

/* finding coefficients of a function with Chebyshev basis expansion
// such that grid's node collocated at Chebyshev extrema.
//
// input explanation: 
// =================
//
// o. *f is the value of function at each point
// o. coeffs are the coefficients
// o. *n is n[3] i.e. number of points in each direction
// 
// some definitions:
// =================
//
// o. REDFT00 is for :
// Y_k = X_0 + (-1)^k*X_(n-1)+2*\sum_{j=1}^{n-2}(X_j*cos(Pi*j*k/(n-1)))
// where Y_k s are coeffes and X_k s are field value. 
// moreover the inverse transformation is the same, i.e. REDFT00.
// o. for each transformation the result is unnormalized and to normalize
// it one must divide the result by N = 2*(n-1).
//
// ->return value: EXIT_SUCCESS
*/
int fftw_3d_ChebyshevExtrema_coeffs(double *const f,double *const coeffs,const int *const n)
{
  fftw_plan p;
  
  p = fftw_plan_r2r_3d(n[0],n[1],n[2],f,coeffs,
            FFTW_REDFT00,FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE);
  fftw_execute(p);
        
  /*Freeing*/
  fftw_destroy_plan(p);
  fftw_cleanup();

  return EXIT_SUCCESS;
}

/* it is same as above, but this one calculates the values of field
// based on given coeffs.
// ->return value: EXIT_SUCCESS.
*/
int fftw_3d_ChebyshevExtrema_values(double *const f,double *const coeffs,const int *const n)
{
  fftw_plan p;
  
  p = fftw_plan_r2r_3d(n[0],n[1],n[2],coeffs,f,
            FFTW_REDFT00,FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE);
  fftw_execute(p);
        
  /*Freeing*/
  fftw_destroy_plan(p);
  fftw_cleanup();

  return EXIT_SUCCESS;
}
