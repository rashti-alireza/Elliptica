/*
// Alireza Rashti
// June 2018
*/

#include "grid.h"

/* making the patches which cover the grid */
int make_patches(Grid_T *const grid)
{
  char *kind;
  Flag_T flg;
  
  /* finding the kind of grid */
  kind = get_parameter_value_S("grid_kind",&flg);
  parameterEr(flg);
  grid->kind = dup_s(kind);
  
  /* allocating and filling patches */
  alloc_patches(grid);
  fill_patches(grid);
  
  /* allocating and filling nodes */
  alloc_nodes(grid);
  make_nodes(grid);
  
  /* filling grid->nn */
  grid->nn = total_nodes_grid(grid);
  
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
  Flag_T flg;
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
    n = (unsigned)get_parameter_value_I("all_Nabc",&flg);
    if (flg != NONE)
      patch->n[0] = patch->n[1] = patch->n[2] = n;
    /* check for override */
    make_keyword_parameter(&ret,name,"n");
    n = (unsigned)get_parameter_value_I(ret.s0,&flg);
    if (flg != NONE)	patch->n[0] = n;
    n = (unsigned)get_parameter_value_I(ret.s1,&flg);
    if (flg != NONE)	patch->n[1] = n;
    n = (unsigned)get_parameter_value_I(ret.s2,&flg);
    if (flg != NONE)	patch->n[2] = n;
    
    assert(patch->n[0] && patch->n[1] && patch->n[2]);
    
    /* filling center */
    make_keyword_parameter(&ret,name,"center");
    patch->c[0] = get_parameter_value_D(ret.s0,&flg);
    parameterEr(flg);
    patch->c[1] = get_parameter_value_D(ret.s1,&flg);
    parameterEr(flg);
    patch->c[2] = get_parameter_value_D(ret.s2,&flg);
    parameterEr(flg);
    
    /* filling size */
    make_keyword_parameter(&ret,name,"size");
    patch->s[0] = get_parameter_value_D(ret.s0,&flg);
    parameterEr(flg);
    patch->s[1] = get_parameter_value_D(ret.s1,&flg);
    parameterEr(flg);
    patch->s[2] = get_parameter_value_D(ret.s2,&flg);
    parameterEr(flg);
    
    /* filling min: min = center-l/2 */
    patch->min[0] = patch->c[0]-patch->s[0]/2;
    patch->min[1] = patch->c[1]-patch->s[1]/2;
    patch->min[2] = patch->c[2]-patch->s[2]/2;
    
    /* filling max: max = center+l/2 */
    patch->max[0] = patch->c[0]+patch->s[0]/2;
    patch->max[1] = patch->c[1]+patch->s[1]/2;
    patch->max[2] = patch->c[2]+patch->s[2]/2;
    
    /* filling flags */
    patch->coordsys = dup_s("Cartesian");
    
    /* collocation */
    c = get_collocation(get_parameter_value_S("all_collocation",&flg));
    if (flg != NONE)
      patch->collocation[0] = patch->collocation[1] = 
          patch->collocation[2] = c;
    /* check for override */
    make_keyword_parameter(&ret,name,"collocation");
    c = get_collocation(get_parameter_value_S(ret.s0,&flg));
    if (flg != NONE)
      patch->collocation[0] = c;
    c = get_collocation(get_parameter_value_S(ret.s1,&flg));
    if (flg != NONE)
      patch->collocation[1] = c;
    c = get_collocation(get_parameter_value_S(ret.s2,&flg));
    if (flg != NONE)
      patch->collocation[2] = c;
    
    assert(patch->collocation[0] != UNDEFINED_COLLOCATION);
    assert(patch->collocation[1] != UNDEFINED_COLLOCATION);
    assert(patch->collocation[2] != UNDEFINED_COLLOCATION);
    
    /* basis */
    b = get_basis(get_parameter_value_S("all_basis",&flg));
    if (flg != NONE)
      patch->basis[0] = patch->basis[1] = patch->basis[2] = b;
    /* check for override */
    make_keyword_parameter(&ret,name,"basis");
    b = get_basis(get_parameter_value_S(ret.s0,&flg));
    if (flg != NONE)
      patch->basis[0] = b;
    b = get_basis(get_parameter_value_S(ret.s1,&flg));
    if (flg != NONE)
      patch->basis[1] = b;
    b = get_basis(get_parameter_value_S(ret.s2,&flg));
    if (flg != NONE)
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
