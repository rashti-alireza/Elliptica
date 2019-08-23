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

/* fourier transformation from real value to complex coeffs:
// f(x) = \sum_{m=-l+1}^{l-1} c(m)*exp(I*m*x), where x is in [0,2*pi],
// and l = n/2+1 (if n is odd, it is rounded down) thus: we have:
// c(m) = 1/(2*pi)*\integral_{0}^{2*pi} f(x)*exp(-I*m*x) dx.
//
// some notes:
// ============
// note: l = n/2 +1 if n is odd it is rounded down.
// note: since f(x) is real we have: c(-m) = c*(m) 
//       and we only calculate c(j), j = 0,...,l-1
// note: argument n is the number of grid points in x direction.
// note: of course, works best when f function is periodic in [0,2*pi]
// note: it allocates memory so at some point you need to free it.
// note: we use trapezoidal rule to carry out the intergral. 
// ->return value: c(m) */
void *r2cft_1d_EquiSpaced_coeffs(const double *const value,const unsigned n)
{
  if (n == 0)
    abortEr("Fourier Transformation: No points!\n");
  
  const unsigned l = n/2+1;/* note: if n is not even, it is rounded down */
  double complex *const coeffs = alloc_double_complex(l);
  const double complex x0 = -2.*I*M_PI/n;/* - included */
  unsigned m;
  
  for (m = 0; m < l; ++m)/* note: l is excluded, otherwise one have aliasing and then error */
  {
    unsigned i;
    
    coeffs[m] = 0;
    for (i = 0; i < n; ++i)
      coeffs[m] += value[i]*cexp(m*i*x0);
    
    coeffs[m] /= n;
  }
  
  return coeffs;
}

/* fourier transformation from complex coeffs to real values:
// f(x) = \sum_{m=-l+1}^{l-1} c(m)*exp(I*m*x), where x is in [0,2*pi],
// and l = n/2+1 (if n is odd, it is rounded down) thus: we have:
// c(m) = 1/(2*pi)*\integral_{0}^{2*pi} f(x)*exp(-I*m*x) dx. 
// ->return value : f(x) */
double *c2rft_1d_EquiSpaced_values(void *const coeffs,const unsigned N)
{
  if (N == 0)
    abortEr("Fourier Transformation: No points!\n");
  
  const double complex *const c = coeffs;
  double *f = alloc_double(N);
  const unsigned l = N/2+1;
  const double complex x0 = 2.*I*M_PI/N;
  unsigned i,j;
  
  for (i = 0; i < N; ++i)
  {
    double complex x = i*x0;
    f[i] = creal(c[0]);
    for (j = 1; j < l; ++j)
      f[i] += 2.*creal(c[j]*cexp((double)j*x));
  }

  return f;
}
