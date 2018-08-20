/*
// Alireza Rashti
// May 2018
*/

#include "memory_alloc.h"

/* adding 2 block of memory for parameter data base 
// and putting the last block to null and 
// returning pointer to one before the last block
*/
void *alloc_parameter(Parameter_T ***const mem)
{
  unsigned i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  pointerEr((*mem));
  
  (*mem)[i] = malloc(sizeof(*(*mem)[i]));
  pointerEr((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}

/* adding 2 block of memory for project data base 
// and putting the last block to null and 
// returning pointer to one before the last block
*/
void *alloc_project(Project_T ***const mem)
{
  unsigned i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  pointerEr((*mem));
  
  (*mem)[i] = malloc(sizeof(*(*mem)[i]));
  pointerEr((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}

/* allocating memory for grid structure.
// there are flags which determine whether or not the grid should 
// be newly allocated or there is already a gird which is recently 
// deleted and can be used readily. note: this function add grid to 
// grids_global; furthermore, the end of grids_global is determined by
// NULL.
// ->return value: a new grid
*/
void *alloc_grid(void)
{
  unsigned i;
  
  for (i = 0; grids_global != 0 && grids_global[i] != 0; i++)
  {
    if (grids_global[i]->status == READY)
    {
      grids_global[i]->status = INUSE;
      return grids_global[i];
    }
  }
  
  grids_global = realloc(grids_global,(i+2)*sizeof(*grids_global));
  pointerEr(grids_global);
  
  grids_global[i] = calloc(1,sizeof(*grids_global[i]));
  pointerEr(grids_global[i]);
  
  grids_global[i+1] = 0;
  
  grids_global[i]->status = INUSE;
  return grids_global[i];
}

/* allocating memory for patches based on type of coord sys */
void alloc_patches(Grid_T *const grid)
{
  if (strcmp_i(grid->kind,"Cartesian_grid"))
    alloc_patches_Cartesian_grid(grid);
  
  /*else if (strcmp_i(grid->kind,"CubedSpherical_grid"))
    //alloc_patches_CubedSpherical_grid(grid); */
  else
    abortEr_s("No such %s kind for grid.\n",grid->kind);
}

/* memory alloc nodes */
void alloc_nodes(Grid_T *const grid)
{
  unsigned *n;
  unsigned i;
  
  FOR_ALL(i,grid->patch)
  {
    unsigned j,U;
    Node_T **node;
    
    n = grid->patch[i]->n;
    grid->patch[i]->node = 
      malloc((n[0]*n[1]*n[2]+1)*sizeof(*grid->patch[i]->node));
    pointerEr(grid->patch[i]->node);
    
    node = grid->patch[i]->node;
    node[n[0]*n[1]*n[2]] = 0;
    
    U = n[0]*n[1]*n[2];
    for (j = 0; j < U; j++)
    {
      node[j] = calloc(1,sizeof(*node[j]));
      pointerEr(node[j]);
    }
    
  }
  
}

/* memory alloc patches for Cartesian type */
static void alloc_patches_Cartesian_grid(Grid_T *const grid)
{
  unsigned Nboxes;/* number of boxes */
  unsigned i;
  
  if (get_parameter("number_of_boxes") == 0)
    abortEr("\"number_of_boxes\" parameter is not defined!\n");
    
  Nboxes = (unsigned) GetParameterI_E("number_of_boxes");
  
  grid->patch = calloc((Nboxes+1),sizeof(*grid->patch));
  pointerEr(grid->patch);
  
  for (i = 0; i < Nboxes; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    pointerEr(grid->patch[i]);
  }
  
}

/* memory allocation for interface struct */
void alloc_interface(Patch_T *const patch)
{
  unsigned i;
  assert(patch);
  
  patch->interface = calloc(FACE_NUM+1,sizeof(*patch->interface));
  pointerEr(patch->interface);
  
  for (i = 0; i < FACE_NUM; i++)
  {
    patch->interface[i] = calloc(1,sizeof(*patch->interface[i]));
    pointerEr(patch->interface[i]);
  }
}

/*
// memory allocation for point struct;
// s is the number of point which is demanded 
// ->return value: pointer to new allocated memory
*/
void *alloc_point(const unsigned s)
{
  Point_T **point;
  unsigned i;
  
  point = calloc(s+1,sizeof(*point));
  pointerEr(point);
  
  for (i = 0; i < s; i++)
  {
    point[i] = calloc(1,sizeof(*point[i]));
    pointerEr(point[i]);
  }
  
  return point;
}

/* allocating 2 block of memory for sFunc_PtoV_T 
// and putting the last block to NULL and returning
// the new available pointer.
// ->return value: a pointer to a ready sFunc_PtoV
*/
void *alloc_sFunc_PtoV(sFunc_PtoV_T ***const mem)
{
  unsigned i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  pointerEr((*mem));
  
  (*mem)[i] = malloc(sizeof(*(*mem)[i]));
  pointerEr((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}

/* allocating 2 block of memory for sFunc_Patch2Pdouble_T 
// and putting the last block to NULL and returning
// the new available pointer.
// ->return value: a pointer to a ready sFunc_Patch2Pdouble
*/
void *alloc_sFunc_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const mem)
{
  unsigned i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  pointerEr((*mem));
  
  (*mem)[i] = calloc(1,sizeof(*(*mem)[i]));
  pointerEr((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}


/* make an empty needle 
// ->return value: pointer to new needle
*/
void *alloc_needle(void)
{
  Needle_T *needle;
  
  needle = calloc(1,sizeof(*needle));
  pointerEr(needle);

  return needle;
}

/* return value-> N*sizeof(double), using calloc*/
double *alloc_double(const unsigned N)
{
  double *d;
  
  d = calloc(N,sizeof(*d));
  pointerEr(d);
  
  return d;
}

/* return value-> M[R][C] double type memory using calloc */
double **alloc_matrix(const long unsigned R,const long unsigned C)
{
  double **M;
  long unsigned row;
  
  M = calloc(R,sizeof(*M));
  pointerEr(M);
  
  for (row = 0; row < R; ++row)
  {
    M[row] = calloc(C,sizeof(*M[row]));
    pointerEr(M[row]);
  }
  
  return M;
}

/* realloc patch->solution_man->solution struct by n units and
// return the pointer to solution[n-1]
// ->return value: patch->solution_man->solution[n-1]
*/
Solution_T *alloc_solution(Patch_T *const patch, const unsigned n)
{
  if (!patch->solution_man)
    return 0;
  
  Solution_T **sol = patch->solution_man->solution;

  sol = realloc(sol,n*sizeof(*sol));
  pointerEr(sol);
  sol[n-1] = calloc(1,sizeof(*sol[n-1]));
  pointerEr(sol[n-1]);
  patch->solution_man->solution = sol;
  
  return sol[n-1];
}
