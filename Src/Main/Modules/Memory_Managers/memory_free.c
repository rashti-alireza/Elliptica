/*
// Alireza Rashti
// May 2018
*/

#include "memory_free.h"

/* freeing 2 dimensions block of memory
// knowing the last column is pointing to null
*/
void free_2d(void *mem0)
{
  int i;
  void **mem = mem0;
    
  for (i = 0; mem[i] != 0; i++)
    free(mem[i]);
    
  free(mem);
    
}

/* freeing 2 dimensions block of memory
// knowing the number of columns is c
*/
void free_matrix(void *mem0, long int c)
{
  long int i;
  void **mem = mem0;
  
  for (i = 0; i < c; i++)
    free(mem[i]);
    
  free(mem);
    
}

