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
// I->Composite_Simpson_1D->a = 10;
// I->Composite_Simpson_1D->b = 20;
// I->Composite_Simpson_1D->n = 11;
// I->Composite_Simpson_1D->f = array;
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

