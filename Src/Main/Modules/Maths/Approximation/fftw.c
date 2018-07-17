/*
// Alireza Rashti
// July 2018
*/

#include "fftw.h"

/* finding NORMALIZED coefficients of a function with Chebyshev basis of 
// first kind expansion
// such that grid's node collocated at Chebyshev extrema.
//
// input explanation: 
// =================
//
// o. *values refers to the value of function at each point
// o. coeffs are the coefficients
// o. n is number of points in each direction
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
*/
void fftw_1d_ChebyshevExtrema_coeffs(double *const values,double *const coeffs,const unsigned n)
{
  fftw_plan p;
  const double Nrm = 2*(n-1);
  const unsigned B = n;
  unsigned i;
  
  p = fftw_plan_r2r_1d((int)n,values,coeffs,FFTW_REDFT00,FFTW_ESTIMATE);
  fftw_execute(p);
        
  /*Freeing*/
  fftw_destroy_plan(p);
  fftw_cleanup();
  
  for (i = 0; i < B; ++i)
    coeffs[i] /= Nrm;
}

/* finding NORMALIZED coefficients of a function with Chebyshev basis of 
// first kind expansion with Chebyshev extrema collocation, 3-D one.
*/
void fftw_3d_ChebyshevExtrema_coeffs(double *const values,double *const coeffs,const unsigned *const n)
{
  fftw_plan p;
  const double Nrm = 8*(n[0]-1)*(n[1]-1)*(n[2]-1);
  const unsigned B = n[0]*n[1]*n[2];
  unsigned ijk;
  
  p = fftw_plan_r2r_3d((int)n[0],(int)n[1],(int)n[2],values,coeffs,
            FFTW_REDFT00,FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE);
  fftw_execute(p);
        
  /*Freeing*/
  fftw_destroy_plan(p);
  fftw_cleanup();
  
  /* normalizing coeffs */
  for (ijk = 0; ijk < B; ++ijk)
    coeffs[ijk] /= Nrm;
}

/* calculates the NORMALIZED values of field based on given coeffs. */
void fftw_3d_ChebyshevExtrema_values(double *const values,double *const coeffs,const unsigned *const n)
{
  fftw_plan p;
  const double Nrm = 8*(n[0]-1)*(n[1]-1)*(n[2]-1);
  const unsigned B = n[0]*n[1]*n[2];
  unsigned ijk;
  
  p = fftw_plan_r2r_3d((int)n[0],(int)n[1],(int)n[2],coeffs,values,
            FFTW_REDFT00,FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE);
  fftw_execute(p);
        
  /*Freeing*/
  fftw_destroy_plan(p);
  fftw_cleanup();
  
  /* normalizing values */
  for (ijk = 0; ijk < B; ++ijk)
    values[ijk] /= Nrm;
}

/* calculates the NORMALIZED values of field based on given coeffs. */
void fftw_1d_ChebyshevExtrema_values(double *const values,double *const coeffs,const unsigned n)
{
  fftw_plan p;
  const double Nrm = 2*(n-1);
  const unsigned B = n;
  unsigned i;
  
  p = fftw_plan_r2r_1d((int)n,coeffs,values,FFTW_REDFT00,FFTW_ESTIMATE);
  fftw_execute(p);
        
  /*Freeing*/
  fftw_destroy_plan(p);
  fftw_cleanup();
  
  for (i = 0; i < B; ++i)
    values[i] /= Nrm;
}
