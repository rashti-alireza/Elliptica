/*
// Alireza Rashti
// May 2018
*/

#include "memory_alloc.h"

/* adding 2 block of memory for parameter data base 
// and puting the last block to null and 
// returning pointer to one before the last block
*/
void *alloc_parameter(Parameter_T ***mem)
{
  int i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  pointerEr((*mem));
  
  (*mem)[i] = malloc(sizeof(*(*mem)[i]));
  pointerEr((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}

/* adding 2 block of memory for project data base 
// and puting the last block to null and 
// returning pointer to one before the last block
*/
void *alloc_project(Project_T ***mem)
{
  int i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  pointerEr((*mem));
  
  (*mem)[i] = malloc(sizeof(*(*mem)[i]));
  pointerEr((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}

/* allocating memory for gird structure.
// there are flags which determine wheather or not the grid should 
// be newly allocated or there is already a gird which is recently 
// deleted and can be used readily. note: this function add grid to 
// grids_global; furthermore, the end of grids_global is determined by
// NULL.
*/
void *alloc_grid(void)
{
  extern Grid_T **grids_global;
  int i;
  
  for (i = 0; girds_global != 0 && grids_global[i] != 0; i++)
  {
    if (grids_global[i]->status == READY)
    {
      grids_global[i]->status = INUSED;
      return grids_global[i];
    }
  }
  
  grids_global = realloc(grids_global,(i+2)*sizeof(*grids_global));
  pointerEr(grids_global);
  
  grids_global[i] = malloc(sizeof(*grids_global[i]));
  pointerEr(grids_global[i]);
  
  grids_global[i+1] = 0;
  
  grids_global[i]->status = INUSED;
  return grids_global[i];
}