/*
// Alireza Rashti
// May 2018
*/

#include "memory_alloc.h"

/* adding 2 block of memory for parameter data base 
// and puting the last block to null and 
// returning pointer to one before the last block
*/
void *alloc_parameter(Parameter_T **mem)
{
  int i;
  
  for (i = 0; mem != 0 && mem[i] != 0 ; i++);
  
  mem = realloc(mem,(i+2)*sizeof(*mem));
  pointerEr(mem);
  
  mem[i] = malloc(sizeof(*mem[i]));
  pointerEr(mem[i]);
  
  mem[i+1] = 0;
  
  return mem[i];
}
