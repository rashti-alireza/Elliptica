/*
// Alireza Rashti
// June 2018
*/

#include "grid.h"

/* making the patches which cover the grid */
int make_patches(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Making the patches ...\n\n");
  
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
