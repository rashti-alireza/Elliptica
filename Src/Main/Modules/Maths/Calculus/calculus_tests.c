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
  
  return EXIT_SUCCESS;
}

/* testing composite simpson's Rule 1D.
// ->return value: if successful => TEST_SUCCESSFUL, otherwise TEST_UNSUCCESSFUL */
static int csr_1d(Grid_T *const grid)
{
  unsigned N = 2019;
  const double a = M_PI;
  const double b = 3./2.*M_PI;
  Integration_T *I = init_integration();
  double *f = alloc_double(N);/* integrant */
  double sf,an;/* resultant */
  double dx;
  double err;/* expected error from theory */
  unsigned i;
  
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
