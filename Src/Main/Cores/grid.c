/*
// Alireza Rashti
// June 2018
*/

#include "grid.h"

/* making the patches which cover the grid */
int make_patches(Grid_T *const grid)
{
  const char *kind;
  
  /* finding the kind of grid */
  kind = GetParameterS_E("grid_kind");
  grid->kind = dup_s(kind);
  
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
  
  return EXIT_SUCCESS;
}

/* filling patch struct */
static void fill_patches(Grid_T *const grid)
{
  if (strcmp_i(grid->kind,"Cartesian_grid"))
    fill_patches_Cartesian_grid(grid);
  
  /*else if (grid->coord,"CubedSpherical")
    //fill_patches_CubedSpherical_grid(grid); 
    */
  else
    abortEr_s("There is no such %s grid kind.\n",grid->kind);
}

/*filling patch struct for Cartesian*/
static void fill_patches_Cartesian_grid(Grid_T *const grid)
{
  unsigned nc;
  char name[20] = {'\0'};
  Collocation_T c;
  Basis_T b;
  unsigned i;
  
  nc = 0;
  FOR_ALL(i,grid->patch)
  {
    struct Ret_S ret;
    Patch_T *const patch = grid->patch[i];
    c = UNDEFINED_COLLOCATION;
    b = UNDEFINED_BASIS;
    unsigned n;
    
    /* filling grid */
    patch->grid = grid;
    
    /* filling patch number */
    patch->pn = i;
    
    /* filling node counter */
    patch->nc = nc;
    
    /* filling inner boundary */
    patch->innerB = 0;
    
    /* filling name */
    sprintf(name,"box%d",i);
    patch->name = dup_s(name);
    
    /* filling n */
    n = (unsigned)GetParameterI_E("all_Nabc");
    if (n != INT_MAX)
      patch->n[0] = patch->n[1] = patch->n[2] = n;
    /* check for override */
    make_keyword_parameter(&ret,name,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
    
    assert(patch->n[0] && patch->n[1] && patch->n[2]);
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    /* filling center */
    make_keyword_parameter(&ret,name,"center");
    patch->c[0] = GetParameterD_E(ret.s0);
    patch->c[1] = GetParameterD_E(ret.s1);
    patch->c[2] = GetParameterD_E(ret.s2);
    
    /* filling size */
    make_keyword_parameter(&ret,name,"size");
    patch->s[0] = GetParameterD_E(ret.s0);
    patch->s[1] = GetParameterD_E(ret.s1);
    patch->s[2] = GetParameterD_E(ret.s2);
    
    /* filling min: min = center-l/2 */
    patch->min[0] = patch->c[0]-patch->s[0]/2;
    patch->min[1] = patch->c[1]-patch->s[1]/2;
    patch->min[2] = patch->c[2]-patch->s[2]/2;
    
    /* filling max: max = center+l/2 */
    patch->max[0] = patch->c[0]+patch->s[0]/2;
    patch->max[1] = patch->c[1]+patch->s[1]/2;
    patch->max[2] = patch->c[2]+patch->s[2]/2;
    
    /* filling flags */
    patch->coordsys = Cartesian;
    
    /* collocation */
    c = get_collocation(GetParameterS("all_collocation"));
    if (c != UNDEFINED_COLLOCATION)
      patch->collocation[0] = patch->collocation[1] = 
      patch->collocation[2] = c;
    /* check for override */
    make_keyword_parameter(&ret,name,"collocation");
    c = get_collocation(GetParameterS(ret.s0));
    if (c != UNDEFINED_COLLOCATION)
      patch->collocation[0] = c;
    c = get_collocation(GetParameterS(ret.s1));
    if (c != UNDEFINED_COLLOCATION)
      patch->collocation[1] = c;
    c = get_collocation(GetParameterS(ret.s2));
    if (c != UNDEFINED_COLLOCATION)
      patch->collocation[2] = c;
    
    assert(patch->collocation[0] != UNDEFINED_COLLOCATION);
    assert(patch->collocation[1] != UNDEFINED_COLLOCATION);
    assert(patch->collocation[2] != UNDEFINED_COLLOCATION);
    
    /* basis */
    b = get_basis(GetParameterS_E("all_basis"));
    if ( b != UNDEFINED_BASIS)
      patch->basis[0] = patch->basis[1] = patch->basis[2] = b;
    /* check for override */
    make_keyword_parameter(&ret,name,"basis");
    b = get_basis(GetParameterS(ret.s0));
    if ( b != UNDEFINED_BASIS)
      patch->basis[0] = b;
    b = get_basis(GetParameterS(ret.s1));
    if ( b != UNDEFINED_BASIS)
      patch->basis[1] = b;
    b = get_basis(GetParameterS(ret.s2));
    if ( b != UNDEFINED_BASIS)
      patch->basis[2] = b;
    
    assert(patch->basis[0] != UNDEFINED_BASIS);
    assert(patch->basis[1] != UNDEFINED_BASIS);
    assert(patch->basis[2] != UNDEFINED_BASIS);
    
    nc +=  patch->n[0]*patch->n[1]*patch->n[2];
  }
  
}

/* making keyword for searching of parameter value */
static void make_keyword_parameter(struct Ret_S *const ret,const char *const box,const char *const needle)
{
  /* for box?_n_? */
  if (strcmp_i(needle,"n"))
  {
    sprintf(ret->s0,"%s_n_a",box);
    sprintf(ret->s1,"%s_n_b",box);
    sprintf(ret->s2,"%s_n_c",box);
  }
  
  /* for box?_center_? */
  else if (strcmp_i(needle,"center"))
  {
    sprintf(ret->s0,"%s_center_a",box);
    sprintf(ret->s1,"%s_center_b",box);
    sprintf(ret->s2,"%s_center_c",box);
  }
  
  /* for box?_size_? */
  else if (strcmp_i(needle,"size"))
  {
    sprintf(ret->s0,"%s_size_a",box);
    sprintf(ret->s1,"%s_size_b",box);
    sprintf(ret->s2,"%s_size_c",box);
  }
  /* for box?_collocation_? */
  else if (strcmp_i(needle,"collocation"))
  {
    sprintf(ret->s0,"%s_collocation_a",box);
    sprintf(ret->s1,"%s_collocation_b",box);
    sprintf(ret->s2,"%s_collocation_c",box);
  }
  /* for box?_basis_? */
  else if (strcmp_i(needle,"basis"))
  {
    sprintf(ret->s0,"%s_basis_a",box);
    sprintf(ret->s1,"%s_basis_b",box);
    sprintf(ret->s2,"%s_basis_c",box);
  }
  else
  {
    abortEr_s("There is no such %s.\n",needle);
  }
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
  tmp_patch.node = 0;
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

