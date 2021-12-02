/*
// Alireza Rashti
// August 2019
*/

#include "misc_functions.h"

/* calculate factorial for integer number/
// ->return value: n*(n-1)*...*1. */
int Factorial(const int n)
{
  if (n < 0)
    Error0("Factorial argument is negative!\n");
  if (n == 0)
    return 1;
    
  return n*Factorial(n-1);
}


/* second derivative of Cheb_Tn. 
// ->return value: second derivative of Tn
*/
double d2Cheb_Tn_dx2(const int n, const double x)
{
  double d = DBL_MAX;
  
  if (n == 0 || n == 1)
    d = 0;
  else if (n == 2)
    d = 4;
  else if (EQL(x,1))
    d = n*n*(n*n-1)/3.0;
  else if (EQL(x,-1))
  {
    if (n%2)
      d = -n*n*(n*n-1)/3.0;
    else
      d = n*n*(n*n-1)/3.0;
  }
  else
  {
    d = n*((n+1)*Cheb_Tn(n,x)-Cheb_Un(n,x))/(x*x-1);
  }
  
  return d;
}

/* ->return value: d(Cheb_Tn(x))/dx */
double dCheb_Tn_dx(const int n,const double x)
{
  return n*Cheb_Un(n-1,x);
}

