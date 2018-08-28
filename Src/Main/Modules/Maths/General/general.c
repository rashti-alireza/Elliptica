/*
// Alireza Rashti
// June 2018
*/

#include "general.h"

/* taking squre root of vector v2-v1 which has l double type components.
// ->return value: root mean square of v2-v1.
*/
double rms(const unsigned n, const double *const v2,const double *const v1)
{
  unsigned i;
  double sum;
  sum = 0;
  
  if (v2 == 0 && v1 == 0) return 0;
  else if (v1 == 0)    	  
    for(i = 0; i < n; i++)
      sum += SQR(v2[i]);
  else if (v2 == 0)       
    for(i = 0; i < n; i++)
      sum += SQR(v1[i]);
  else
    for(i = 0; i < n; i++)
      sum += SQR(v2[i]-v1[i]);
    
  sum = sqrt(sum);
  
  return sum;
}

/* taking root means squre of vector v2-v1 which has l double type components
// and l is of order of long unsigned.
// ->return value: root mean square of v2-v1
*/
long double rmsL(const long unsigned n, const double *const v2, const double *const v1)
{
  unsigned long i;
  long double sum;
  sum = 0;
  
  if (v2 == 0 && v1 == 0) return 0;
  else if (v1 == 0)    	  
    for(i = 0; i < n; i++)
      sum += SQR(v2[i]);
  else if (v2 == 0)       
    for(i = 0; i < n; i++)
      sum += SQR(v1[i]);
  else
    for(i = 0; i < n; i++)
      sum += SQR(v2[i]-v1[i]);
    
  sum = sqrtl(sum);
  
  return sum;
}

/* taking dot product of two v1 and v2 vector with n components
// ->return value : v2.v1
*/
double dot(const unsigned n, const double *const v2,const double *const v1)
{
  unsigned i;
  double d = 0;
  
  for (i = 0; i < n; i++)
    d += v2[i]*v1[i];
  
  return d;
}

/* taking absolute value of v
// ->return value: absolute(v)
*/
double ABS(const double v)
{
  return v > 0 ? v : -v;
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