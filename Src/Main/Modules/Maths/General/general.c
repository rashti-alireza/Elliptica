/*
// Alireza Rashti
// June 2018
*/

#include "general.h"

/* taking squre root of vector v which has l components double version*/
double rms(const int n, double *v)
{
  long int i;
  double sum;
  sum = 0;
  
  for(i = 0; i < n; i++)
  {
    sum += SQR(v[i]);
  }
  
  sum = sqrt(sum);
  
  return sum;
}

/* taking squre root of vector v which has l components 
// long double version
*/
long double rmsL(const long int n, double *v)
{
  long int i;
  long double sum;
  sum = 0;
  
  for(i = 0; i < n; i++)
  {
    sum += SQR(v[i]);
  }
  
  sum = sqrt(sum);
  
  return sum;
}
