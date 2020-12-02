/*
// Alireza Rashti
// June 2018
*/

#include "complex_numbers.h"

/* calloc memory for double complex */
void *alloc_double_complex(const Uint N)
{
  double complex *f = calloc(N,sizeof(*f));
  IsNull(f);
  
  return f;
}

