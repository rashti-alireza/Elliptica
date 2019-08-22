/*
// Alireza Rashti
// August 2019
*/

#include "special_functions.h"

/* calculate factorial for integer number/
// ->return value: n*(n-1)*...*1. */
int Factorial(const int n)
{
  if (n < 0)
    abortEr("Factorial argument is negative!\n");
  if (n == 0)
    return 1;
    
  return n*Factorial(n-1);
}
