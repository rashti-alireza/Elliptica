/*
// Alireza Rashti
// June 2019
*/

#include "integration.h"

/* synopsis:
// =========
//
// ** filling the integration struct **
// Integration_T *I = init_integration();
//
// *** for example if you want to use Composite Simpson's Rule 1D: ***
// I->type = "Composite Simpson's Rule 1D"
// I->Composite_Simpson_1D->a = 1.5;// e.g
// I->Composite_Simpson_1D->b = M_PI;// e.g
// I->Composite_Simpson_1D->n = 11;// e.g
// I->Composite_Simpson_1D->f = array;
//
// *** for example if you want to use Gaussian Quadrature Chebyshev Extrema: ***
// I->type = "Gaussian Quadrature Chebyshev Extrema"
// I->GQ_ChebyshevExtrema->f = array;
// I->GQ_ChebyshevExtrema->n = 10;// the dimension of array
//
// *** for example if you want to use Gaussian Quadrature Lobatto method: ***
// I->type = "Gaussian Quadrature Lobatto"
// I->GQ_Lobatto->f = array;
// I->GQ_Lobatto->n = 10;// the dimension of array
// 
// ** planning the appropriate function for integration **
// plan_integration(I);
//
// ** evaluating integration **
// double integral = execute_integration(I);
//
// ** freeing **
// free_integration(I);
*/


/* initializing an Integration_T struct with calloc.
// ->return value: a pristine struct */
Integration_T *init_integration(void)
{
  Integration_T *I = calloc(1,sizeof(*I));
  pointerEr(I);
  
  return I;
}

/* ->return value: integral */
double execute_integration(Integration_T *const I)
{
  return I->integration_func(I);
}

/* given the information it decides how to perform the integral */
void plan_integration(Integration_T *const I)
{
  if (strcmp_i(I->type,"Composite Simpson's Rule 1D"))
  {
    I->integration_func = Composite_Simpson_1D;
  }
  else if (strcmp_i(I->type,"Gaussian Quadrature Chebyshev Extrema"))
  {
    I->integration_func = GaussQuadrature_ChebyshevExtrema;
  }
  else if (strcmp_i(I->type,"Gaussian Quadrature Lobatto"))
  {
    I->integration_func = GaussQuadrature_Lobatto;
  }
  else
    abortEr(NO_OPTION);
}

/* free the integral struct */
void free_integration(Integration_T *I)
{
  if (!I)
    return;
    
  free(I);
}

/* Composite Simpson's Rule in 1D
// ->return value: the result of inegral */
static double Composite_Simpson_1D(Integration_T *const I)
{
  if (I->Composite_Simpson_1D->n % 2 != 1)
    abortEr("Composite Simpson's Rule requires odd number of points.\n");
    
  const double h = (I->Composite_Simpson_1D->b-I->Composite_Simpson_1D->a)/(I->Composite_Simpson_1D->n-1);
  const double *const f = I->Composite_Simpson_1D->f;
  double i0 = 0,i1 = 0,i2 = 0;
  unsigned j;
  
  I->err = fabs(I->Composite_Simpson_1D->b-I->Composite_Simpson_1D->a)/180*pow(h,4)*100;
  
  i0 = f[0]+f[I->Composite_Simpson_1D->n-1];
  for (j = 1; j < I->Composite_Simpson_1D->n-1; ++j)
  {
    if (j%2 == 0)
      i2 += f[j];
    else
      i1 += f[j];
  }
  
  return h*(i0+2*i2+4*i1)/3;
}

/* performing the integral \integral_{-1}^{1} dx f(x)/sqrt(1-x^2) 
// using Guassian Quadrature method.
// note: to get a good result for this method, you need many points
//       unless f(x) is a polynomial which works best.
// note: the collocation points for f(x) are Chebyshev Extrema
// note: the integration is from -1 to 1 then the expected order of f(x)
//       is f(-1) = f[0] and f(1) = f[n-1].
// ->return value: \integral_{-1}^{1} dx f(x)/sqrt(1-x^2) */
static double GaussQuadrature_ChebyshevExtrema(Integration_T *const I)
{
  double i0 = 0;
  const double *const f = I->GQ_ChebyshevExtrema->f;
  const unsigned n      = I->GQ_ChebyshevExtrema->n;
  const double   w      = M_PI/(n-1);
  double err = M_PI;
  unsigned i;
  
  err /= Factorial(2*(int)n);
  err *= L_inf(n,f);/* approximately */
  err /= pow(2,2*n-1);
  I-> err = err;
  
  for (i = 1; i <= n-2; ++i)
    i0 += f[i];
  i0 += (f[0]+f[n-1])/2;
  i0 *= w;
  
  return i0;
}

/* performing the integral \integral_{-1}^{1} dx f(x)
// using Guassian Quadrature Lobatto's method.
// note: the collocation points for f(x) are zeros of d(Legendre(x,n-1))/dx
// note: the integration is from -1 to 1 then the expected order of f(x)
//       is f(-1) = f[0] and f(1) = f[n-1].
// ->return value: \integral_{-1}^{1} f(x)dx */
static double GaussQuadrature_Lobatto(Integration_T *const I)
{
  double i0 = 0;
  const double *const f = I->GQ_Lobatto->f;
  const unsigned n      = I->GQ_Lobatto->n;
  const int ni          = (int)n;
  double (*w)(const double x, const unsigned n) = Lobbatto_weight_function;
  double err;
  unsigned i;
  
  err = n*pow(n-1,3)/(2*n-1);
  err *= pow(2,2*n-1)/Factorial(2*ni-2);
  err *= Factorial(ni-2);
  err *= Factorial(ni-2);
  err /= Factorial(2*ni-2);
  err *= Factorial(ni-2);
  err /= Factorial(2*ni-2);
  err *= Factorial(ni-2);
  err *= L_inf(n,f);/* approximately */
  I-> err = err;

  for (i = 1; i <= n-2; ++i)
    i0 += w(Lobbatto_root_function(i-1,n-1),n)*f[i];
  i0 += 2*(f[0]+f[n-1])/(n*(n-1));
  
  return i0;
}
