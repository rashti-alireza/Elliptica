/*
// Alireza Rashti
// July 2018
*/

#include "fourier_transformation.h"


/* 1-D regular Fourier transformation -as opposed to 
// fast Fourier transformation.
// getting values and number of values,
// it finds NORMALIZED coefficients of a function with Chebyshev basis of 
// first kind expansion with Chebyshev Extrema collocation.
//
// what it calculates are c's which stand for coeffs as follows :
// v(i) = c(0) + (-1)^i*c(n-1) + 2 \sum_{j=1}^{n-2} c(j)*cos(i*j*PI/(n-1))
// the algorithm is simple : one needs to write the above c's in 
// terms of v's as the following:
// c(i)*scale = v(0) + (-1)^i*v(n-1) + 2 \sum_{j=1}^{n-2} v(j)*cos(i*j*PI/(n-1))
// in wich scale is 2*(n-1).
// note: the coeffs are NORMALIZED.
*/
void rft_1d_ChebyshevExtrema_coeffs(double *const values ,double *const coeffs,const unsigned n)
{
  unsigned i,j;
  const double SIGN[2] = {1.0,-1.0};
  const double th0 = M_PI/(n-1);
  const double scale = 2*(n-1);
  
  for (i = 0; i < n; ++i)
  {
    double th = i*th0;
    double c = 0;/* coeffs */
    for (j = 1; j < n-1; ++j)
      c += values[j]*cos(j*th);
    
    c *= 2;
    c += values[0]+SIGN[i%2]*values[n-1];
    coeffs[i] = c/scale;
  }
}


/* 1-D regular Fourier transformation -as opposed to 
// fast Fourier transformation.
// getting values and number of values,
// it finds NORMALIZED coefficients of a function with Chebyshev basis of 
// first kind expansion with Chebyshev Nodes collocation
//
// what it calculates are c's which stand for coeffs as follows :
// v(i) = c(0)+2*\sum_{j=1}^{n-1} c(j)*cos(j*PI(i+1/2)/n)
// the algorithm is simple : one needs to write the above c's in 
// terms of v's as the following:
// c(i)*scale = 2*\sum_{j=0}^{n-1} v(j)*cos(i*PI(j+1/2)/n)
// in wich scale is 2n.
// note: the coeffs are NORMALIZED.
*/
void rft_1d_ChebyshevNodes_coeffs(double *const values ,double *const coeffs,const unsigned n)
{
  unsigned i,j;
  const double th0 = M_PI/n;
  const double scale = n;
  
  for (i = 0; i < n; ++i)
  {
    double th = i*th0;
    double c = 0;/* coeffs */
    for (j = 0; j < n; ++j)
      c += values[j]*cos(th*(j+0.5));
      
    coeffs[i] = c/scale;
  }
}
