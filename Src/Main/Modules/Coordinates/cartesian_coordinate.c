/*
// Alireza Rashti
// June 2018
*/

#include "cartesian_coordinate.h"

/* making value of coords. it is a general function for Cartesian type */
void make_nodes_Cartesian_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const unsigned U = patch->nn;
  const unsigned *const n = patch->n;
  unsigned i,j,k,l;
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  initialize_collocation_struct(patch,&coll_s[2],2);
  
  for (l = 0; l < U; l++)
  {
    double *x = patch->node[l]->x;
    
    IJK(l,n,&i,&j,&k);
    x[0] = point_value(i,&coll_s[0]);
    x[1] = point_value(j,&coll_s[1]);
    x[2] = point_value(k,&coll_s[2]);
    
    /* since X and x are the same we have: */
    patch->node[l]->X = x;
  }
}

/* making Jacobian transformation for Cartesian coord. */
void make_JacobianT_Cartesian_coord(Patch_T *const patch)
{
  patch->JacobianT->j      = JT_Cartesian_patch;
}

/* Jacobian transformation for Cartesian patch.
// ->return value: dq2/dq1
*/
double JT_Cartesian_patch(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  double j;
  
  if (q2_e%3 == q1_e%3)/* e.g. _y_%3 = 1 and _b_%3 = 1 */
    j = 1;
  else
    j = 0;
    
  UNUSED(patch);
  UNUSED(p);
  
  return j;
}

/*filling patch struct for Cartesian*/
void fill_patches_Cartesian_grid(Grid_T *const grid)
{
  char name[20] = {'\0'};
  Collocation_T c;
  Basis_T b;
  unsigned i;
  
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
    
    /* filling inner boundary */
    patch->innerB = 0;
    
    /* filling name */
    sprintf(name,"box%d",i);
    patch->name = dup_s(name);
    
    /* filling n */
    patch->n[0] = (unsigned)GetParameterI("n_a");
    patch->n[1] = (unsigned)GetParameterI("n_b");
    patch->n[2] = (unsigned)GetParameterI("n_c");
    
    /* check for override */
    make_keyword_parameter(&ret,name,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;

    if(patch->n[0] == INT_MAX)
      abortEr("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      abortEr("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      abortEr("n_c could not be set.\n");
      
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
    patch->collocation[0] = get_collocation(GetParameterS("collocation_a"));
    patch->collocation[1] = get_collocation(GetParameterS("collocation_b"));
    patch->collocation[2] = get_collocation(GetParameterS("collocation_c"));
  
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
    patch->basis[0] = get_basis(GetParameterS("basis_a"));
    patch->basis[1] = get_basis(GetParameterS("basis_b"));
    patch->basis[2] = get_basis(GetParameterS("basis_c"));
  
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
    
  }
  
}

/* populating properties of the box at the middle of left NS */
void populate_left_NS_central_box(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_centeral_box",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  sprintf(var,"grid%u_left_centeral_box_n_a",grid->gn);
  patch->n[0] = (unsigned)GetParameterI_E(var);
  
  sprintf(var,"grid%u_left_centeral_box_n_b",grid->gn);
  patch->n[1] = (unsigned)GetParameterI_E(var);
  
  sprintf(var,"grid%u_left_centeral_box_n_c",grid->gn);
  patch->n[2] = (unsigned)GetParameterI_E(var);
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling size */
  sprintf(var,"grid%u_left_centeral_box_size_a",grid->gn);
  patch->s[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_centeral_box_size_b",grid->gn);
  patch->s[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_centeral_box_size_c",grid->gn);
  patch->s[2] = GetParameterDoubleF_E(var);
  
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
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of right box in for single neutron star */
void populate_right_box_sns(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_box",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  sprintf(var,"grid%u_right_box_n_a",grid->gn);
  patch->n[0] = (unsigned)GetParameterI_E(var);
  
  sprintf(var,"grid%u_right_box_n_b",grid->gn);
  patch->n[1] = (unsigned)GetParameterI_E(var);
  
  sprintf(var,"grid%u_right_box_n_c",grid->gn);
  patch->n[2] = (unsigned)GetParameterI_E(var);
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_right_box_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_box_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_box_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling size */
  sprintf(var,"grid%u_right_box_size_a",grid->gn);
  patch->s[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_box_size_b",grid->gn);
  patch->s[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_box_size_c",grid->gn);
  patch->s[2] = GetParameterDoubleF_E(var);
  
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
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}


/* populating properties of the filling box in cubed spherical grid */
void populate_filling_box_CubedSpherical(Grid_T *const grid,const unsigned pn,const Flag_T side)
{
  Patch_T *const patch = grid->patch[pn];
  unsigned n;
  double l;/* length */
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  
  switch(side)
  {
    case UP:
    /* filling name */
    sprintf(name,"grid%u_filling_box_up",grid->gn);
    patch->name = dup_s(name);
    
    sprintf(var,"grid%u_surrounding_box_length",grid->gn);
    l = GetParameterDoubleF_E(var);
    
    /* filling center */
    patch->c[0] = 0;
    patch->c[1] = 0;
    patch->c[2] = 3./4.*l;
    
    /* filling size */
    patch->s[0] = l;
    patch->s[1] = 2*l;
    patch->s[2] = 0.5*l;
  
    /* check for override n*/
    sprintf(var,"Outermost0");
    make_keyword_parameter(&ret,var,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)   patch->n[0] = n;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)   patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)   patch->n[2] = n/2;
    
    break;
    case DOWN:
    /* filling name */
    sprintf(name,"grid%u_filling_box_down",grid->gn);
    patch->name = dup_s(name);
    
    sprintf(var,"grid%u_surrounding_box_length",grid->gn);
    l = GetParameterDoubleF_E(var);
    
    /* filling center */
    patch->c[0] = 0;
    patch->c[1] = 0;
    patch->c[2] = -3./4.*l;
    
    /* filling size */
    patch->s[0] = l;
    patch->s[1] = 2*l;
    patch->s[2] = 0.5*l;
  
    /* check for override n*/
    sprintf(var,"Outermost0");
    make_keyword_parameter(&ret,var,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)   patch->n[0] = n;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)   patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)   patch->n[2] = n/2;
    
    break;
    case BACK:
    /* filling name */
    sprintf(name,"grid%u_filling_box_back",grid->gn);
    patch->name = dup_s(name);
    
    sprintf(var,"grid%u_surrounding_box_length",grid->gn);
    l = GetParameterDoubleF_E(var);
    
    /* filling center */
    patch->c[0] = -3./4.*l;
    patch->c[1] = 0;
    patch->c[2] = 0;
    
    /* filling size */
    patch->s[0] = 0.5*l;
    patch->s[1] = 2*l;
    patch->s[2] = 2*l;
  
    /* check for override n*/
    sprintf(var,"Outermost0");
    make_keyword_parameter(&ret,var,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)   patch->n[0] = n/2;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)   patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)   patch->n[2] = n;
    
    break;
    case FRONT:
    /* filling name */
    sprintf(name,"grid%u_filling_box_front",grid->gn);
    patch->name = dup_s(name);
    
    sprintf(var,"grid%u_surrounding_box_length",grid->gn);
    l = GetParameterDoubleF_E(var);
    
    /* filling center */
    patch->c[0] = 3./4.*l;
    patch->c[1] = 0;
    patch->c[2] = 0;
    
    /* filling size */
    patch->s[0] = 0.5*l;
    patch->s[1] = 2*l;
    patch->s[2] = 2*l;
    
    /* check for override n*/
    sprintf(var,"Outermost0");
    make_keyword_parameter(&ret,var,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)   patch->n[0] = n/2;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)   patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)   patch->n[2] = n;
    
    break;
    default:
      abortEr(NO_OPTION);
  }
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
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
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of the box at the middle of right NS  */
void populate_right_NS_central_box(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_centeral_box",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  sprintf(var,"grid%u_right_centeral_box_n_a",grid->gn);
  patch->n[0] = (unsigned)GetParameterI_E(var);
  
  sprintf(var,"grid%u_right_centeral_box_n_b",grid->gn);
  patch->n[1] = (unsigned)GetParameterI_E(var);
  
  sprintf(var,"grid%u_right_centeral_box_n_c",grid->gn);
  patch->n[2] = (unsigned)GetParameterI_E(var);
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling size */
  sprintf(var,"grid%u_right_centeral_box_size_a",grid->gn);
  patch->s[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_centeral_box_size_b",grid->gn);
  patch->s[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_centeral_box_size_c",grid->gn);
  patch->s[2] = GetParameterDoubleF_E(var);
  
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
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* memory alloc patches for Cartesian grid type */
void alloc_patches_Cartesian_grid(Grid_T *const grid)
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
