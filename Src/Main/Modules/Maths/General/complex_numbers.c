/*
// Alireza Rashti
// June 2018
*/

#include "complex_numbers.h"

/* calloc memory for double complex */
void *alloc_double_complex(const unsigned N)
{
  double complex *f = calloc(N,sizeof(*f));
  pointerEr(f);
  
  return f;
}

