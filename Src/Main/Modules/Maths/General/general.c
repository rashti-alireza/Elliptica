/*
// Alireza Rashti
// June 2018
*/

#include "general.h"

/* taking squre root of vector v2-v1 which has l components double version*/
double rms(const int n, double *v2,double *v1)
{
  int i;
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
long double rmsL(const long int n, double *v2, double *v1)
{
  long int i;
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
double dot(const int n, double *v2,double *v1)
{
  int i;
  double d = 0;
  
  for (i = 0; i < n; i++)
    d += v2[i]*v1[i];
  
  return d;
}
