/*
// Alireza Rashti
// June 2018
*/

#include "grid.h"

/* making the patches which cover the grid */
int make_patches(Grid_T *grid)
{
  char *coord;
  Flag_T flg;
  
  /* finding coord sys name in grid */
  coord = get_parameter_value_S("coord_sys",&flg);
  parameterEr(flg);
  grid->kind = strdup(coord);
  
  /* allocating and filling patches */
  alloc_patches(grid);
  fill_patches(grid);
  
  /* allocating and filling nodes */
  alloc_nodes(grid);
  fill_nodes(grid);
  
  return EXIT_SUCCESS;
}

/* filling node struc*/
static void fill_nodes(Grid_T *grid)
{
  if (!strcmp(grid->kind,"Cartesian"))
    fill_nodes_cartesian(grid);
  
  //else if (grid->coord,"CubedSpherical")
    //fill_nodes_CubedSpherical(grid);
}

/* filling node struc for Cartesian*/
static void fill_nodes_cartesian(Grid_T *grid)
{
  /* filling nodes based on coord sys type */
  make_coordinates(grid);
  
}

/* filling patch struct */
static void fill_patches(Grid_T *grid)
{
  if (!strcmp(grid->kind,"Cartesian"))
    fill_patches_cartesian(grid);
  
  //else if (grid->coord,"CubedSpherical")
    //fill_patches_CubedSpherical(grid);
}

/*filling patch struct for Cartesian*/
static void fill_patches_cartesian(Grid_T *grid)
{
  int i;
  char name[20] = {'\0'};
  Flag_T flg;
  
  for_all_patches_macro(i,grid)
  {
    struct Ret_S ret;
    Patch_T *const patch = grid->patch[i];
    
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
    patch->collocation = get_parameter_value_S(name,&flg);
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
