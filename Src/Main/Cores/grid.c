/*
// Alireza Rashti
// June 2018
*/

#include "grid.h"

/* making the patches which cover the grid */
int make_patches(Grid_T *const grid)
{
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
  
  else if (strcmp_i(grid->kind,"BNS_Projective_grid"))
    fill_patches_BNS_Projective_grid(grid); 
  
  else if (strcmp_i(grid->kind,"BNS_Spherical_grid"))
    fill_patches_BNS_Spherical_grid(grid); 
    
  else
    abortEr_s("There is no such %s grid kind.\n",grid->kind);
}

/* filling patch struct for BNS_Projective_grid */
static void fill_patches_BNS_Projective_grid(Grid_T *const grid)
{
  const unsigned N_outermost_split = (unsigned) GetParameterI_E("Number_of_Outermost_Split");
  unsigned pn,i;
  
  pn = 0;
  populate_left_NS_central_box(grid,pn++);
  populate_left_NS_hemisphere_up(grid,pn++);
  populate_left_NS_hemisphere_down(grid,pn++);
  populate_left_NS_surrounding_up(grid,pn++);
  populate_left_NS_surrounding_down(grid,pn++);
  for (i = 0; i < N_outermost_split; i++)
    populate_left_outermost(grid,pn++,i);
    
  populate_right_NS_central_box(grid,pn++);
  populate_right_NS_hemisphere_up(grid,pn++);
  populate_right_NS_hemisphere_down(grid,pn++);
  populate_right_NS_surrounding_up(grid,pn++);
  populate_right_NS_surrounding_down(grid,pn++);
  for (i = 0; i < N_outermost_split; i++)
    populate_right_outermost(grid,pn++,i);
  
}

/* filling patch struct for BNS_Spherical_grid */
static void fill_patches_BNS_Spherical_grid(Grid_T *const grid)
{
  const unsigned N_outermost_split = (unsigned) GetParameterI_E("Number_of_Outermost_Split");
  unsigned pn,i;
  
  pn = 0;
  populate_left_NS_sphere(grid,pn++);
  populate_left_NS_surrounding_sphere(grid,pn++);
  for (i = 0; i < N_outermost_split; i++)
    populate_left_outermost(grid,pn++,i);
    
  populate_right_NS_sphere(grid,pn++);
  populate_right_NS_surrounding_sphere(grid,pn++);
  for (i = 0; i < N_outermost_split; i++)
    populate_right_outermost(grid,pn++,i);
  
}

/*filling patch struct for Cartesian*/
static void fill_patches_Cartesian_grid(Grid_T *const grid)
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

/* populating properties of patch for outermost left */
static void populate_left_outermost(Grid_T *const grid,const unsigned pn,const unsigned outermost_n)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_outermost%u",grid->gn,outermost_n);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"Outermost%u",outermost_n);
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_outermost%u_R1",grid->gn,outermost_n);
  patch->CoordSysInfo->R1 = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_outermost%u_R2",grid->gn,outermost_n);
  patch->CoordSysInfo->R2 = GetParameterDoubleF_E(var);
  
  assert(GRT(patch->CoordSysInfo->R2,patch->CoordSysInfo->R1));
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = StereographicSphereLeft;
  
  /* collocation */
  patch->collocation[0] = Chebyshev_Nodes;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
  
}

/* populating properties of patch for right outermost */
static void populate_right_outermost(Grid_T *const grid,const unsigned pn,const unsigned outermost_n)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_outermost%u",grid->gn,outermost_n);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"Outermost%u",outermost_n);
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_outermost%u_R1",grid->gn,outermost_n);
  patch->CoordSysInfo->R1 = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_outermost%u_R2",grid->gn,outermost_n);
  patch->CoordSysInfo->R2 = GetParameterDoubleF_E(var);
  
  assert(GRT(patch->CoordSysInfo->R2,patch->CoordSysInfo->R1));
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = StereographicSphereRight;
  
  /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Nodes;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
  
}

/* populating properties of patch for NS left hemisphere up */
static void populate_left_NS_hemisphere_up(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double *R1_array,*R2_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_NS_hemisphere_up",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_NS_R1_up",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  sprintf(var,"grid%u_left_NS_R2_up",grid->gn);
  R2_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->R1_f = R1;
  patch->CoordSysInfo->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R2->v[ij0] = R2_array[ij0];
      R1->v[ij0] = R1_array[ij0];
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereUp;
  
  /* collocation */
  patch->collocation[0] = Chebyshev_Nodes;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
  
}

/* populating properties of patch for NS left hemisphere down */
static void populate_left_NS_hemisphere_down(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double *R1_array,*R2_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_NS_hemisphere_down",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_NS_R1_down",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  sprintf(var,"grid%u_left_NS_R2_down",grid->gn);
  R2_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->R1_f = R1;
  patch->CoordSysInfo->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R2->v[ij0] = R2_array[ij0];
      R1->v[ij0] = R1_array[ij0];
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereDown;
  
/* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Nodes;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for NS right hemisphere up */
static void populate_right_NS_hemisphere_up(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double *R1_array,*R2_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_NS_hemisphere_up",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_NS_R1_up",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  sprintf(var,"grid%u_right_NS_R2_up",grid->gn);
  R2_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->R1_f = R1;
  patch->CoordSysInfo->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R2->v[ij0] = R2_array[ij0];
      R1->v[ij0] = R1_array[ij0];
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereUp;
  
/* collocation */
  patch->collocation[0] = Chebyshev_Nodes;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for NS right hemisphere down */
static void populate_right_NS_hemisphere_down(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double *R1_array,*R2_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_NS_hemisphere_down",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_NS_R1_down",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  sprintf(var,"grid%u_right_NS_R2_down",grid->gn);
  R2_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->R1_f = R1;
  patch->CoordSysInfo->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R2->v[ij0] = R2_array[ij0];
      R1->v[ij0] = R1_array[ij0];
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereDown;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Nodes;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for left NS's surrounding up */
static void populate_left_NS_surrounding_up(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double R2_const,*R1_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_NS_surrounding_up",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_NS_Surrounding_R2",grid->gn);
  R2_const = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_R2_up",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->R1_f = R1;
  patch->CoordSysInfo->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R1->v[ij0] = R1_array[ij0];
      R2->v[ij0] = R2_const;
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereUp;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Nodes;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for left NS's surrounding down */
static void populate_left_NS_surrounding_down(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double R2_const,*R1_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_NS_surrounding_down",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_NS_Surrounding_R2",grid->gn);
  R2_const = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_R2_down",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->R1_f = R1;
  patch->CoordSysInfo->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R1->v[ij0] = R1_array[ij0];
      R2->v[ij0] = R2_const;
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereDown;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Nodes;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for right NS's surrounding up */
static void populate_right_NS_surrounding_up(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double R2_const,*R1_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_NS_surrounding_up",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_NS_Surrounding_R2",grid->gn);
  R2_const = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_R2_up",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->R1_f = R1;
  patch->CoordSysInfo->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R1->v[ij0] = R1_array[ij0];
      R2->v[ij0] = R2_const;
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereUp;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Nodes;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for right_NS's surrounding down */
static void populate_right_NS_surrounding_down(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double R2_const,*R1_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_NS_surrounding_down",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_NS_Surrounding_R2",grid->gn);
  R2_const = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_R2_down",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->R1_f = R1;
  patch->CoordSysInfo->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R1->v[ij0] = R1_array[ij0];
      R2->v[ij0] = R2_const;
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereDown;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Nodes;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of the box at the middle of left NS */
static void populate_left_NS_central_box(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  unsigned n;
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
  sprintf(var,"grid%u_left_centeral_box_n_abc",grid->gn);
  n = (unsigned)GetParameterI_E(var);
  patch->n[0] = patch->n[1] = patch->n[2] = n;
  
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
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

/* populating properties of the box at the middle of right NS */
static void populate_right_NS_central_box(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  unsigned n;
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
  sprintf(var,"grid%u_right_centeral_box_n_abc",grid->gn);
  n = (unsigned)GetParameterI_E(var);
  patch->n[0] = patch->n[1] = patch->n[2] = n;
  
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
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

/* populating properties of patch for left NS */
static void populate_left_NS_sphere(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_radius",0,patch,NO);
  Field_T *R2 = add_field("R2_radius",0,patch,NO);
  double *R1_array,*R2_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  Collocation_T c;
  Basis_T b;
  unsigned n,j,k;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_NS",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_NS_R1",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  sprintf(var,"grid%u_left_NS_R2",grid->gn);
  R2_array = GetParameterArrayF_E(var);
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (j = 0; j < patch->n[1]; ++j)
    for (k = 0; k < patch->n[2]; ++k)
    {
      unsigned ijk = L(patch->n,0,j,k);
      R2->v[ijk] = R2_array[ijk];
      R1->v[ijk] = R1_array[ijk];
    }
  
  /* filling min */
  patch->min[0] = 0;
  patch->min[1] = 0;
  patch->min[2] = 0;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = M_PI;
  patch->max[2] = 2*M_PI;
  
  /* filling flags */
  patch->coordsys = Spherical;
  
  /* collocation */
  patch->collocation[0] = get_collocation(GetParameterS_E("collocation_a"));
  patch->collocation[1] = get_collocation(GetParameterS_E("collocation_b"));
  patch->collocation[2] = get_collocation(GetParameterS_E("collocation_c"));

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
  patch->basis[0] = get_basis(GetParameterS_E("basis_a"));
  patch->basis[1] = get_basis(GetParameterS_E("basis_b"));
  patch->basis[2] = get_basis(GetParameterS_E("basis_c"));

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

/* populating properties of patch for right NS */
static void populate_right_NS_sphere(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_radius",0,patch,NO);
  Field_T *R2 = add_field("R2_radius",0,patch,NO);
  double *R1_array,*R2_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  Collocation_T c;
  Basis_T b;
  unsigned n,k,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_NS",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_NS_R1",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  sprintf(var,"grid%u_right_NS_R2",grid->gn);
  R2_array = GetParameterArrayF_E(var);
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (j = 0; j < patch->n[1]; ++j)
    for (k = 0; k < patch->n[2]; ++k)
    {
      unsigned ijk = L(patch->n,0,j,k);
      R2->v[ijk] = R2_array[ijk];
      R1->v[ijk] = R1_array[ijk];
    }
  
  /* filling min */
  patch->min[0] = 0;
  patch->min[1] = 0;
  patch->min[2] = 0;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = M_PI;
  patch->max[2] = 2*M_PI;
  
  /* filling flags */
  patch->coordsys = Spherical;
  
  /* collocation */
  patch->collocation[0] = get_collocation(GetParameterS_E("collocation_a"));
  patch->collocation[1] = get_collocation(GetParameterS_E("collocation_b"));
  patch->collocation[2] = get_collocation(GetParameterS_E("collocation_c"));

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
  patch->basis[0] = get_basis(GetParameterS_E("basis_a"));
  patch->basis[1] = get_basis(GetParameterS_E("basis_b"));
  patch->basis[2] = get_basis(GetParameterS_E("basis_c"));

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

/* populating properties of patch for left NS's surrounding */
static void populate_left_NS_surrounding_sphere(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_radius",0,patch,NO);
  Field_T *R2 = add_field("R2_radius",0,patch,NO);
  double R2_const,*R1_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  Collocation_T c;
  Basis_T b;
  unsigned n,k,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_NS_surrounding",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_NS_Surrounding_R2",grid->gn);
  R2_const = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_R2",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (j = 0; j < patch->n[1]; ++j)
    for (k = 0; k < patch->n[2]; ++k)
    {
      unsigned ijk = L(patch->n,0,j,k);
      R1->v[ijk] = R1_array[ijk];
      R2->v[ijk] = R2_const;
    }
  
  /* filling min */
  patch->min[0] = 0;
  patch->min[1] = 0;
  patch->min[2] = 0;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = M_PI;
  patch->max[2] = 2*M_PI;
  
  /* filling flags */
  patch->coordsys = Spherical;
  
  /* collocation */
  patch->collocation[0] = get_collocation(GetParameterS_E("collocation_a"));
  patch->collocation[1] = get_collocation(GetParameterS_E("collocation_b"));
  patch->collocation[2] = get_collocation(GetParameterS_E("collocation_c"));

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
  patch->basis[0] = get_basis(GetParameterS_E("basis_a"));
  patch->basis[1] = get_basis(GetParameterS_E("basis_b"));
  patch->basis[2] = get_basis(GetParameterS_E("basis_c"));

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

/* populating properties of patch for right NS's surrounding */
static void populate_right_NS_surrounding_sphere(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_radius",0,patch,NO);
  Field_T *R2 = add_field("R2_radius",0,patch,NO);
  double R2_const,*R1_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  Collocation_T c;
  Basis_T b;
  unsigned n,k,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_NS_surrounding",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
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
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_NS_Surrounding_R2",grid->gn);
  R2_const = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_R2",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (j = 0; j < patch->n[1]; ++j)
    for (k = 0; k < patch->n[2]; ++k)
    {
      unsigned ijk = L(patch->n,0,j,k);
      R1->v[ijk] = R1_array[ijk];
      R2->v[ijk] = R2_const;
    }
  
  /* filling min */
  patch->min[0] = 0;
  patch->min[1] = 0;
  patch->min[2] = 0;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = M_PI;
  patch->max[2] = 2*M_PI;
  
  /* filling flags */
  patch->coordsys = Spherical;
    
  /* collocation */
  patch->collocation[0] = get_collocation(GetParameterS_E("collocation_a"));
  patch->collocation[1] = get_collocation(GetParameterS_E("collocation_b"));
  patch->collocation[2] = get_collocation(GetParameterS_E("collocation_c"));

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
  patch->basis[0] = get_basis(GetParameterS_E("basis_a"));
  patch->basis[1] = get_basis(GetParameterS_E("basis_b"));
  patch->basis[2] = get_basis(GetParameterS_E("basis_c"));

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

