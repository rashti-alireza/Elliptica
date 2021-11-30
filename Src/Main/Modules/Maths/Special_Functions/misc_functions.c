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

/* Chebyshev polynomial of second kind Un(x). x MUST be normalized value.
// ->return value: Un(x)
*/
double Cheb_Un(const int n, const double x)
{
  double u = DBL_MAX;
  
  if (n == 0) 
    u = 1;
  else if (EQL(x,1))
    u = n+1;
  else if (EQL(x,-1)) 
  {
    if (n%2)
      u = -n-1;
    else
      u = n+1;
  }  
  else
  {
    double th = acos(x);
    u = sin((n+1)*th)/sin(th);
  }
  
  return u;
}

/* Chebyshev polynomial of first kind Tn(x). x MUST be normalized value.
// ->return value: Tn(x)
*/
double Cheb_Tn(const int n, const double x)
{
  double t = DBL_MAX;
  
  if (n == 0)
    t = 1;
  else if (EQL(x,1))
    t = 1;
  else if (EQL(x,-1))
  {
    if (n%2)
      t = -1;
    else
      t = 1;
  }
  else
  {
    double th = acos(x);
    t = cos(n*th);
  }
  
  return t;
}

/* second derivative of Cheb_Tn. 
// ->return value: second derivative of Tn
*/
double d2T_dx2(const int n, const double x)
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
double dT_dx(const int n,const double x)
{
  return n*Cheb_Un(n-1,x);
}

