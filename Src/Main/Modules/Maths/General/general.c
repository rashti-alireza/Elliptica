/*
// Alireza Rashti
// June 2018
*/

#include "general.h"

/* taking squre root of vector v2-v1 which has l components double version*/
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

/* taking squre root of vector v2-v1 which has l components 
// long double version
*/
long double rmsL(const unsigned long n, const double *const v2, const double *const v1)
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
// ->return value: U
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
