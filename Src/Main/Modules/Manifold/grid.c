/*
// Alireza Rashti
// June 2018
*/

#include "grid.h"

/* making the patches which cover the grid */
int make_patches(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Making the patches ...\n");
  
  unsigned p;
  
  /* allocating and filling basics of patches */
  alloc_patches(grid);
  fill_patches(grid);
  
  /* allocating and filling nodes */
  alloc_nodes(grid);
  make_nodes(grid);
  
  /* allocating and making Jacobian coordinate transformation */
  make_JacobianT(grid);
  
  /* filling grid->nn */
  grid->nn = total_nodes_grid(grid);
  /* filling grid->np */
  grid->np = countf(grid->patch);
  
  /* test printing coords */
  if (test_print(PRINT_COORDS))
    pr_coords(grid);
 
  /* print some info */
  for (p = 0; p < grid->np; ++p)
  {
    const Patch_T *patch = grid->patch[p];
    char str[1000];
    
    printf("|--> %s:\n",patch->name);
    printf("     |--> Resolution  = %ux%ux%u\n",patch->n[0],patch->n[1],patch->n[2]);
    printf("     |--> Coord. Sys. = %s\n",coord_sys_str(patch,str));
    printf("     |--> Collocation = %s\n",collocation_str(patch,str));
    printf("     |--> Bases       = %s\n",bases_str(patch,str));
    
  }
  printf("} Making the patches ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
  return EXIT_SUCCESS;
}

/* ->return value: string contains info about the coord sys of this patch */
static char *coord_sys_str(const Patch_T *const patch,char *const str)
{
  assert(patch);
  assert(str);
  
  if (patch->coordsys == Cartesian)
    sprintf(str,"Cartesian");
  else if (patch->coordsys == Spherical)
    sprintf(str,"Spherical");
  else if (patch->coordsys == CubedSpherical)
    sprintf(str,"Cubed Spherical");
  else
    sprintf(str,"Not defined!");
    
  return str;  
}

/* ->return value: string contains info about the bases of this patch */
static char *bases_str(const Patch_T *const patch,char *const str)
{
  assert(patch);
  assert(str);
  unsigned i,l;
  
  str[0] = '\0';
  for (i = 0; i < 3; ++i)
  {
    if (patch->basis[i] == Chebyshev_Tn_BASIS)
      strcat(str,"Cheb. Tn(x), ");
    else
      strcat(str,"Not defined!");
  }
  
  l = (unsigned)strlen(str);
  str[l]   = '\0';
  str[l-1] = '\0';
  str[l-2] = '\0';
  
  return str;  
}

/* ->return value: string contains info about collocation points of this patch */
static char *collocation_str(const Patch_T *const patch,char *const str)
{
  assert(patch);
  assert(str);
  unsigned i,l;
  
  str[0] = '\0';
  for (i = 0; i < 3; ++i)
  {
    if (patch->collocation[i] == EquiSpaced)
      strcat(str,"EquiSpaced, ");
    else if (patch->collocation[i] == Chebyshev_Extrema)
      strcat(str,"Cheb. Extrema, ");
    else if (patch->collocation[i] == Chebyshev_Nodes)
      strcat(str,"Cheb. Nodes, ");
    else
      strcat(str,"Not defined!");
  }
  
  l = (unsigned)strlen(str);
  str[l]   = '\0';
  str[l-1] = '\0';
  str[l-2] = '\0';
  
  return str;  
}

/* filling patch struct */
static void fill_patches(Grid_T *const grid)
{
  if (strcmp_i(grid->kind,"Cartesian_grid"))
    fill_patches_Cartesian_grid(grid);
  
  else if (strcmp_i(grid->kind,"BNS_Spherical_grid"))
    fill_patches_BNS_Spherical_grid(grid); 
    
  else if (strcmp_i(grid->kind,"BNS_CubedSpherical_grid"))
    fill_patches_BNS_CubedSpherical_grid(grid); 
    
  else if (strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
    fill_patches_BBN_CubedSpherical_grid(grid); 
    
  else if (strcmp_i(grid->kind,"SNS_CubedSpherical+Box_grid"))
    fill_patches_SNS_CubedSpherical_Box_grid(grid); 
    
  else if (strcmp_i(grid->kind,"SNS_CubedSpherical_grid"))
    fill_patches_SNS_CubedSpherical_grid(grid); 
    
  else if (strcmp_i(grid->kind,"SBH_CubedSpherical_grid"))
    fill_patches_SBH_CubedSpherical_grid(grid); 
    
  else
    abortEr_s("There is no such '%s' grid kind.\n",grid->kind);
}

/* check if all of houseKs have been marked */
void check_houseK(Patch_T *const patch)
{
  Interface_T **const interface = patch->interface;
  const unsigned nf = countf(interface);
  Node_T *node;
  unsigned i,f;
  
  for (f = 0; f < nf; f++)
    for (i = 0; i < interface[f]->np; i++)
      if (interface[f]->point[i]->houseK == 0)
      {
        node = patch->node[interface[f]->point[i]->ind];
        double *x = node->x;
        fprintf(stderr,"This point(%f,%f,%f) has not been found.\n",
                        x[0],x[1],x[2]);
        abortEr("Incomplete function.\n");
      }
}

/* setting all of houseK flags in Point_T to zero for given patch */
void flush_houseK(Patch_T *const patch)
{
  Interface_T **const interface = patch->interface;
  const unsigned nf = countf(interface);
  unsigned i,f;
  
  for (f = 0; f < nf; f++)
    for (i = 0; i < interface[f]->np; i++)
      interface[f]->point[i]->houseK = 0;
}

/* making a temporary patch for thread safety purposes */
Patch_T make_temp_patch(const Patch_T *const patch)
{
  Patch_T tmp_patch;

  tmp_patch = *patch;
  tmp_patch.nfld = 0;
  //tmp_patch.pn = UINT_MAX;
  tmp_patch.grid = 0;
  tmp_patch.name = 0;
  //tmp_patch.node = patch->node;
  tmp_patch.interface = 0;
  tmp_patch.pool = 0;
  tmp_patch.solving_man = 0;

  return tmp_patch;
}

/* freeing the temporary made patch */
void free_temp_patch(Patch_T *const patch)
{
  if (patch->pool != 0)
    free(patch->pool);
}

/* free the whole date base of grid */
void free_grid_db(void)
{
  unsigned i;
  
  i = 0;
  while (grids_global != 0 && grids_global[i] != 0)
  {
    Grid_T *grid = grids_global[i];
    free_grid(grid);/* since the last grid goes to i = 0, 
                    // don't increament i */          
  }
  
  free_2d(grids_global);
  grids_global = 0;
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
  if (i == 0)
    grids_global[i]->gn = i;
  else
    grids_global[i]->gn = grids_global[i-1]->gn+1;
    
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

/* free the given grid completely */
void free_grid(Grid_T *grid)
{
  unsigned p,ijk,nn,f,i,ng;
  
  if (!grid)
    return;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    nn             = patch->nn;
    
    _free(patch->name);
    
    if (patch->coordsys != Cartesian)
      for (ijk = 0; ijk < nn; ++ijk)
        _free(patch->node[ijk]->X);
    
    if (patch->node)    
      free_2d_mem(patch->node,nn);
    
    for (f = 0; f < patch->nfld; ++f)
    {
      Field_T *field = patch->pool[f];
      free_field(field);
    }
    _free(patch->pool);
    _free(patch->JacobianT);
    free_patch_interface(patch);
    if (patch->solving_man)
    {
      free_patch_SolMan_jacobian(patch);
      free_patch_SolMan_method_Schur(patch);
      _free(patch->solving_man->field_eq);
      _free(patch->solving_man->bc_eq);
      _free(patch->solving_man->jacobian_field_eq);
      _free(patch->solving_man->jacobian_bc_eq);
      free_2d_mem(patch->solving_man->field_name,patch->solving_man->nf);
      free(patch->solving_man);
    }
  }
  free_2d_mem(grid->patch,grid->np);
  _free(grid->kind);
  
  /* shrink the grids_global */
  Grid_T *last_grid = 0;
  
  for (ng = 0; grids_global != 0 && grids_global[ng] != 0; ng++);
  
  for (i = 0; i < ng; ++i)
  {
    if (grid == grids_global[i])
    {
      last_grid       = grids_global[ng-1];
      assert(last_grid);
      grids_global[i] = last_grid;
      
      grids_global = realloc(grids_global,ng*sizeof(*grids_global));
      pointerEr(grids_global);
      grids_global[ng-1] = 0;
      
      break;
    }
  }
  free(grid);
}

/* free the patch completely */
void free_patch(Patch_T *patch)
{
  unsigned ijk,nn,f;
  
  if (!patch)
    return;
  
  nn = patch->nn;
  
  _free(patch->name);
  
  if (patch->coordsys != Cartesian)
    for (ijk = 0; ijk < nn; ++ijk)
      _free(patch->node[ijk]->X);
  
  if (patch->node)    
    free_2d_mem(patch->node,nn);
  
  for (f = 0; f < patch->nfld; ++f)
  {
    Field_T *field = patch->pool[f];
    free_field(field);
  }
  _free(patch->pool);
  _free(patch->JacobianT);
  free_patch_interface(patch);
  if (patch->solving_man)
  {
    free_patch_SolMan_jacobian(patch);
    free_patch_SolMan_method_Schur(patch);
    _free(patch->solving_man->field_eq);
    _free(patch->solving_man->bc_eq);
    _free(patch->solving_man->jacobian_field_eq);
    _free(patch->solving_man->jacobian_bc_eq);
    free_2d_mem(patch->solving_man->field_name,patch->solving_man->nf);
    free(patch->solving_man);
  }
  
  free(patch);
}


