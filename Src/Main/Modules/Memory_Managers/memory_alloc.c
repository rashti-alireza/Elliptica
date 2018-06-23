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
  
  for (i = 0; grids_global != 0 && grids_global[i] != 0; i++)
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

/* allocating memory for patches based on type of coord sys */
void alloc_patches(Grid_T *grid)
{
  if (!strcmp(grid->kind,"Cartesian"))
    alloc_patches_cartesian(grid);
  
  //else if (grid->coord,"CubedSpherical")
    //alloc_patches_CubedSpherical(grid);
}

/* allocating memory for nodes based on type of coord sys */
void alloc_nodes(Grid_T *grid)
{
  if (!strcmp(grid->kind,"Cartesian"))
    alloc_nodes_cartesian(grid);
  
  //else if (grid->coord,"CubedSpherical")
    //alloc_nodes_CubedSpherical(grid);
}

/* memory alloc nodes for cartesian type */
static void alloc_nodes_cartesian(Grid_T *grid)
{
  int *n;
  int i;
  
  for_all_patches_macro(i,grid)
  {
    int j,U;
    Node_T **node;
    
    n = grid->patch[i]->n;
    grid->patch[i]->node = 
      malloc((n[0]*n[1]*n[2]+1)*sizeof(*grid->patch[i]->node));
    pointerEr(grid->patch[i]->node);
    
    node = grid->patch[i]->node;
    U = countf(node);
    for (j = 0; j < U; j++)
    {
      node[j] = malloc(sizeof(*node[j]));
      pointerEr(node[j]);
    }
    
  }
  
}

/* memory alloc patches for cartesian type */
static void alloc_patches_cartesian(Grid_T *grid)
{
  int Nboxes;// number of boxes
  int i;
  Flag_T flg;
  
  if (get_parameter("number_of_boxes") == 0)
    abortEr("number_of_boxes parameter is not defined!\n");
    
  Nboxes = get_parameter_value_I("number_of_boxes",&flg);
  parameterEr(flg);
  
  grid->patch = malloc((Nboxes+1)*sizeof(*grid->patch));
  pointerEr(grid->patch);
  grid->patch[Nboxes] = 0;
  
  for (i = 0; i < Nboxes; i++)
  {
    grid->patch[i] = malloc(sizeof(*grid->patch[i]));
    pointerEr(grid->patch[i]);
  }
  
}
