/*
// Alireza Rashti
// June 2019
*/

#include "calculus_tests.h"

/* test integral */
int integration_tests(Grid_T *const grid)
{
  int status;
  
  if (DO)
  {
    printf("Integration test: Composite Simpson's Rule 1D => \n");
    status = csr_1d(grid);
    check_test_result(status);
  }
  if (DO)
  {
    printf("Integration test: Gaussian Quadrature Chebyshev Extrema => \n");
    status = GQ_ChebExtrema(grid);
    check_test_result(status);
  }
  if (DO)
  {
    printf("Integration test: Gaussian Quadrature Lobatto method => \n");
    status = GQ_Lobatto(grid);
    check_test_result(status);
  }
  
  return EXIT_SUCCESS;
}


/* testing Gaussian Quadrature Chebyshev Extrema integration.
// ->return value: if successful => TEST_SUCCESSFUL */
static int GQ_ChebExtrema(Grid_T *const grid)
{
  unsigned N = 0;
  Integration_T *I = init_integration();
  const char *const par = GetParameterS_E("Test_Integration");
  double *f;
  double sf,an;/* resultant */
  double t0,x;
  unsigned i;
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    N = (unsigned)atoi(s);
    _free(s);
  }
  else
    N = 20000;
 
  f = alloc_double(N);/* integrant */
  t0 = M_PI/(N-1);
  I->type = "Gaussian Quadrature Chebyshev Extrema";
  plan_integration(I);

  /* [-1,1] */
  for (i = 0; i < N; ++i)
  {
    x    = -cos(i*t0);
    f[i] = sqrt(1-SQR(x))*(sin(x)+pow(x,2));/* integral sin(x)+x^2 dx */
  }
    
  an = 2./3.;/* analytic answer from -1 to 1 */

  I->GQ_ChebyshevExtrema->n = N;
  I->GQ_ChebyshevExtrema->f = f;
  sf = execute_integration(I);
  free(f);
  free_integration(I);
  
  printf("Numeric = %0.15f, Analytic = %0.15f, N = %u\n=> ",sf,an,N);
  
  UNUSED(grid);
  return TEST_SUCCESSFUL;
}

/* testing Gaussian Quadrature Lobatto method integration.
// ->return value: if successful => TEST_SUCCESSFUL */
static int GQ_Lobatto(Grid_T *const grid)
{
  unsigned N = 0;
  Integration_T *I = init_integration();
  const char *const par = GetParameterS_E("Test_Integration");
  double *f;
  double sf,an;/* resultant */
  double x;
  unsigned i;
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    N = (unsigned)atoi(s);
    _free(s);
  }
  else
    N = 14;
 
  f = alloc_double(N);/* integrant */
  I->type = "Gaussian Quadrature Lobatto";
  plan_integration(I);

  /* [-1,1] */
  x = -1.; f[0]   = sqrt(1-SQR(x))*(sin(x)+pow(x,2));
  x = 1. ; f[N-1] = sqrt(1-SQR(x))*(sin(x)+pow(x,2));
  for (i = 1; i <= N-2; ++i)
  {
    x    = Lobbatto_root_function(i-1,N-1);
    f[i] = sqrt(1-SQR(x))*(sin(x)+pow(x,2));/* integral sqrt(1-x^2)*(sin(x)+x^2) dx */
  }
    
  an = 0.39269908169872414;/* analytic answer from -1 to 1 */

  I->GQ_Lobatto->n = N;
  I->GQ_Lobatto->f = f;
  sf = execute_integration(I);
  free(f);
  free_integration(I);
  
  printf("Numeric = %0.15f, Analytic = %0.15f, N = %u\n=> ",sf,an,N);
  
  UNUSED(grid);
  return TEST_SUCCESSFUL;
}

/* testing composite simpson's Rule 1D.
// ->return value: if successful => TEST_SUCCESSFUL, otherwise TEST_UNSUCCESSFUL */
static int csr_1d(Grid_T *const grid)
{
  unsigned N = 0;
  const double a = M_PI;
  const double b = 3./2.*M_PI;
  Integration_T *I = init_integration();
  const char *const par = GetParameterS_E("Test_Integration");
  double *f;/* integrant */
  double sf,an;/* resultant */
  double dx;
  double err;/* expected error from theory */
  unsigned i;
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    N = (unsigned)atoi(s);
    _free(s);
  }
  else
    N = 2019;
 
  f = alloc_double(N);
  UNUSED(grid);
  /* testing composite simpson's Rule 1D */
  dx = (b-a)/(N-1);
  for (i = 0; i < N; ++i)
    f[i] = sin(a+i*dx)+pow(a+i*dx,2);/* sin(x)+x^2 */
    
  an = 23.546635705237353;/* from Pi to 3/2*Pi */
  I->type = "Composite Simpson's Rule 1D";
  I->Composite_Simpson_1D->b = b;
  I->Composite_Simpson_1D->a = a;
  I->Composite_Simpson_1D->n = N;
  I->Composite_Simpson_1D->f = f;
  plan_integration(I);
  sf = execute_integration(I);
  err = fabs((b-a)/180*pow(dx,4));
  free(f);
  free_integration(I);
  
  printf("Numeric = %0.15f, Analytic = %0.15f, Expected error = %g\n=> ",sf,an,err);
  
  if (!EQL(fabs(sf-an),err))
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}
