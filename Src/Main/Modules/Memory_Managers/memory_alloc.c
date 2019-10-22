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
  
  (*mem)[i] = calloc(1,sizeof(*(*mem)[i]));
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
// this function add grid to grids_global; 
// furthermore, the end of grids_global is determined by NULL.
// ->return value: a new grid
*/
void *alloc_grid(void)
{
  unsigned i;
  
  for (i = 0; grids_global != 0 && grids_global[i] != 0; i++);
  
  grids_global = realloc(grids_global,(i+2)*sizeof(*grids_global));
  pointerEr(grids_global);
  
  /* allocate new grid */
  grids_global[i] = calloc(1,sizeof(*grids_global[i]));
  pointerEr(grids_global[i]);
  /* set grid number */
  grids_global[i]->gn = i;
  /* determine the last grid */
  grids_global[i+1] = 0;
  
  return grids_global[i];
}

/* allocating memory for patches based on type of coord sys */
void alloc_patches(Grid_T *const grid)
{
  if (strcmp_i(grid->kind,"Cartesian_grid"))
    alloc_patches_Cartesian_grid(grid);
  else if (strcmp_i(grid->kind,"BNS_Spherical_grid"))
    alloc_patches_BNS_Spherical_grid(grid);
  else if (strcmp_i(grid->kind,"BNS_CubedSpherical_grid"))
    alloc_patches_BNS_CubedSpherical_grid(grid);
  else if (strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
    alloc_patches_BBN_CubedSpherical_grid(grid);
  else if (strcmp_i(grid->kind,"SNS_CubedSpherical+Box_grid"))
    alloc_patches_SNS_CubedSpherical_Box_grid(grid);
  else if (strcmp_i(grid->kind,"SNS_CubedSpherical_grid"))
    alloc_patches_SNS_CubedSpherical_grid(grid);
  else if (strcmp_i(grid->kind,"SBH_CubedSpherical_grid"))
    alloc_patches_SBH_CubedSpherical_grid(grid);
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
double **alloc_2D_double(const long unsigned R,const long unsigned C)
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

/* given row and column of a asked matrix it allocates matrix 
// according to the given type using calloc and return the result.
// note : if either of col or row is zero, it returns null.
// ->return value: new empty matrix. null if something not defined correctly.
*/
Matrix_T *alloc_matrix(const Matrix_SF_T type_e,const long row,const long col)
{
  Matrix_T *m = 0;
  
  /* returns null if either of row or col is zero */
  if (row == 0 || col == 0)
    return m;
  
  if (type_e == UNDEF_SF)
    return m;
    
  m = calloc(1,sizeof(*m));
  pointerEr(m);
  m->row = row;
  m->col = col;
  
  switch(type_e)
  {
    case REG_SF:
      m->reg_f = 1;
      m->reg->A = alloc_2D_double((long unsigned)row,(long unsigned)col);
      break;
    case TRI_SF:
      m->tri_f = 1;
      break;
    case CCS_SF:
      m->ccs_f = 1;
      break;
    case CRS_SF:
      m->crs_f = 1;
      break;
    case TRI_L_SF:
      m->tri_l_f = 1;
      break;
    case CCS_L_SF:
      m->ccs_l_f = 1;
      break;
    case CRS_L_SF:
      m->crs_l_f = 1;
      break;
    default:
      abortEr("The specified type is undefined.\n");
  }
  
  return m;
}

/* calloc one sewing.
// ->return value = memory for 1 sewing struct.
*/
Sewing_T *alloc_sewing(void)
{
  Sewing_T *sewing = calloc(1, sizeof(*sewing));
  pointerEr(sewing);
  
  return sewing;
}

/* calloc memory for double complex */
void *alloc_double_complex(const unsigned N)
{
  double complex *f = calloc(N,sizeof(*f));
  pointerEr(f);
  
  return f;
}
