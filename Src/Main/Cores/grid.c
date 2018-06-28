/*
// Alireza Rashti
// June 2018
*/

#include "grid.h"

/* making the patches which cover the grid */
int make_patches(Grid_T *grid)
{
  char *kind;
  Flag_T flg;
  
  /* finding the kind of grid */
  kind = get_parameter_value_S("grid_kind",&flg);
  parameterEr(flg);
  grid->kind = strdup(kind);
  
  /* allocating and filling patches */
  alloc_patches(grid);
  fill_patches(grid);
  
  /* allocating and filling nodes */
  alloc_nodes(grid);
  fill_nodes(grid);
  
  return EXIT_SUCCESS;
}

/* filling patch struct */
static void fill_patches(Grid_T *grid)
{
  if (strcmp_i(grid->kind,"Cartesian_grid"))
    fill_patches_Cartesian_grid(grid);
  
  //else if (grid->coord,"CubedSpherical")
    //fill_patches_CubedSpherical_grid(grid);
  else
    abortEr_s("There is no such %s grid kind.\n",grid->kind);
}

/*filling patch struct for Cartesian*/
static void fill_patches_Cartesian_grid(Grid_T *grid)
{
  int i;
  char name[20] = {'\0'};
  Flag_T flg;
  
  FOR_ALL(i,grid->patch)
  {
    struct Ret_S ret;
    Patch_T *const patch = grid->patch[i];
    
    /* filling grid */
    patch->grid = grid;
    
    /* filling patch number */
    patch->pn = i;
    
    /* filling inner boundary */
    patch->innerB = 0;
    
    /* filling name */
    sprintf(name,"box%d",i);
    patch->name = strdup(name);
    
    /* filling n */
    make_keyword_parameter(&ret,name,"n");
    patch->n[0] = get_parameter_value_I(ret.s0,&flg);
    parameterEr(flg);
    patch->n[1] = get_parameter_value_I(ret.s1,&flg);
    parameterEr(flg);
    patch->n[2] = get_parameter_value_I(ret.s2,&flg);
    parameterEr(flg);
    
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
    sprintf(name,"box%d_collocation",i);
    patch->coordsys = strdup("Cartesian");
    patch->collocation = get_collocation(get_parameter_value_S(name,&flg));
    parameterEr(flg);
    
  }
  
}

/* making keyword for searching of parameter value */
static void make_keyword_parameter(struct Ret_S *ret,char *box,char *needle)
{
  /* for box?_n_? */
  if (!strcmp(needle,"n"))
  {
    sprintf(ret->s0,"%s_n_a",box);
    sprintf(ret->s1,"%s_n_b",box);
    sprintf(ret->s2,"%s_n_c",box);
  }
  
  /* for box?_center_? */
  else if (!strcmp(needle,"center"))
  {
    sprintf(ret->s0,"%s_center_a",box);
    sprintf(ret->s1,"%s_center_b",box);
    sprintf(ret->s2,"%s_center_c",box);
  }
  
  /* for box?_size_? */
  else if (!strcmp(needle,"size"))
  {
    sprintf(ret->s0,"%s_size_a",box);
    sprintf(ret->s1,"%s_size_b",box);
    sprintf(ret->s2,"%s_size_c",box);
  }
  else
  {
    abortEr_s("There is no such %s.\n",needle);
  }
}

/* setting all of houseK flags in Point_T to zero for given patch */
void flush_houseK(Patch_T *patch)
{
  Interface_T **const interface = patch->interface;
  const int nf = countf(interface);
  int i,f;
  
  for (f = 0; f < nf; f++)
    for (i = 0; i < interface[f]->np; i++)
      interface[f]->point[i]->houseK = 0;
}