/*
// Alireza Rashti
// June 2018
*/

#include "grid.h"

/* making the patches which cover the grid */
int make_patches(Grid_T *const grid)
{
  FUNC_TIC
  
  Uint p;
  
  /* allocate patches.
  // since split cubed spherical this is deprecated.
  // now, allocation of patches must be happened as they are filled.
  // also if you want to use previous grid kind, you must modify
  // their functions to meet this condition. */
  if(0) alloc_patches(grid);
  
  /* filling patches parameter */
  fill_patches(grid);
  
  /* allocating and filling nodes */
  alloc_nodes(grid);
  make_nodes(grid);
  
  /* allocating and making Jacobian coordinate transformation */
  make_JacobianT(grid);
  
  /* test printing coords */
  if (test_print(PRINT_COORDS))
  {
    pr_coords(grid);
    Pr_Field_T *pr  = init_PrField(grid);
    pr->folder = Pgets("top_directory");
    pr->cycle  = 0;
    pr_fields(pr);
    free_PrField(pr);
  }
 
  /* print some info */
  if(Pcmps("grid_verbose","yes"))
  for (p = 0; p < grid->np; ++p)
  {
    const Patch_T *patch = grid->patch[p];
    char str[1000];
    
    printf(Pretty0"%s:\n",patch->name);
    printf("     Resolution  = %ux%ux%u\n",patch->n[0],patch->n[1],patch->n[2]);
    printf("     Coord. Sys. = %s\n",coord_sys_str(patch,str));
    printf("     Collocation = %s\n",collocation_str(patch,str));
    printf("     Bases       = %s\n",bases_str(patch,str));
  }
  
  FUNC_TOC
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
  Uint i,l;
  
  str[0] = '\0';
  for (i = 0; i < 3; ++i)
  {
    if (patch->basis[i] == Chebyshev_Tn_BASIS)
      strcat(str,"Cheb. Tn(x), ");
    else
      strcat(str,"Not defined!");
  }
  
  l = (Uint)strlen(str);
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
  Uint i,l;
  
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
  
  l = (Uint)strlen(str);
  str[l]   = '\0';
  str[l-1] = '\0';
  str[l-2] = '\0';
  
  return str;  
}

/* allocating and filling patch struct for grid->patch */
static void fill_patches(Grid_T *const grid)
{
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_BHBH ||
      grid->kind == Grid_SplitCubedSpherical_SNS  ||
      grid->kind == Grid_SplitCubedSpherical_SBH
     )
    fill_patches_Split_CubedSpherical_grid(grid); 

  /* else if (strcmp_i(grid->kind,"Cartesian_grid"))
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
    fill_patches_SBH_CubedSpherical_grid(grid); */
    
  else
    Errors("There is no such '%s' grid kind.\n",Pgets("grid_kind"));
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
  tmp_patch.fields = 0;
  tmp_patch.solving_man = 0;

  return tmp_patch;
}

/* freeing the temporary made patch */
void free_temp_patch(Patch_T *const patch)
{
  if (patch->fields != 0)
    free(patch->fields);
}

/* free the whole date base of grid */
void free_grid_db(void)
{
  Uint i;
  
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
  Uint i;
  
  for (i = 0; grids_global != 0 && grids_global[i] != 0; i++);
  
  grids_global = realloc(grids_global,(i+2)*sizeof(*grids_global));
  IsNull(grids_global);
  
  /* allocate new grid */
  grids_global[i] = calloc(1,sizeof(*grids_global[i]));
  IsNull(grids_global[i]);
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
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_BHBH ||
      grid->kind == Grid_SplitCubedSpherical_SNS  ||
      grid->kind == Grid_SplitCubedSpherical_SBH
     )
    alloc_patches_Split_CubedSpherical_grid(grid);
  
  /*else if (strcmp_i(grid->kind,"Cartesian_grid"))
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
    alloc_patches_SBH_CubedSpherical_grid(grid);*/
  else
    Errors("No such %s kind for grid.\n",Pgets("grid_kind"));
}

/* free the given grid completely */
void free_grid(Grid_T *grid)
{
  Uint p,i,ng;
  
  if (!grid)
    return;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    free_patch(patch);
  }
  Free(grid->patch);
  
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
      IsNull(grids_global);
      grids_global[ng-1] = 0;
      
      break;
    }
  }
  free(grid);
}

/* free the patch completely */
void free_patch(Patch_T *patch)
{
  Uint ijk,nn,f;
  
  if (!patch)
    return;
  
  nn = patch->nn;
  
  Free(patch->name);
  
  if (patch->node)
  {
    if (patch->coordsys != Cartesian)
      for (ijk = 0; ijk < nn; ++ijk)
        Free(patch->node[ijk]->X);
    
    free_2d_mem(patch->node,nn);
    patch->node = 0;
  }
  for (f = 0; f < patch->nfld; ++f)
  {
    Field_T *field = patch->fields[f];
    free_field(field);
  }
  Free(patch->fields);
  Free(patch->JacobianT);
  free_patch_interface(patch);
  if (patch->solving_man)
  {
    free_patch_SolMan_jacobian(patch);
    free_patch_SolMan_method_Schur(patch);
    Free(patch->solving_man->field_eq);
    Free(patch->solving_man->bc_eq);
    Free(patch->solving_man->jacobian_field_eq);
    Free(patch->solving_man->jacobian_bc_eq);
    free_2d_mem(patch->solving_man->field_name,patch->solving_man->nf);
    free(patch->solving_man);
  }
  
  free(patch);
}

extern Parameter_T **parameters_global;
/* free all paramters that used for grid initialization
// which start with grid[0-9]+_. */
void free_grid_params(const Grid_T *const grid)
{
  if (!grid)
    return;
  
  char suffix[STR_LEN1] = {'\0'};
  Uint i,np;
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
  
  sprintf(suffix,PATCH_NAME_PRT_P_"",grid->gn);/* parameters related to this grid */
  for (i = 0; i < np;)/* no increment */
  {
    if (strstr(parameters_global[i]->lv,suffix))
    {
      /* note: the last par is put in palce of removed par
      // so don't increment i */
      free_parameter(parameters_global[i]->lv);
      np--;
    }
    else
      i++;
  }
}



