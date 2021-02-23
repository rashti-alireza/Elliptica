/*
// Alireza Rashti
// March 2018
*/

#include "spherical_coordinate.h"

/* making value of coords. it is a general function for Spherical type */
void make_nodes_Spherical_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const Uint U = patch->nn;
  const Uint *const n = patch->n;
  const Field_T *const R1_field = patch->fields[Ind("R1_radius")];
  const Field_T *const R2_field = patch->fields[Ind("R2_radius")];
  double R1,R2;
  const double *const c = patch->c;/* center of origine translated */
  Uint i,j,k,l;
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  initialize_collocation_struct(patch,&coll_s[2],2);
  
  for (l = 0; l < U; l++)
  {
    double *X = alloc_double(3);
    double *x = patch->node[l]->x;
    double r;
    
    ijk_to_i_j_k(l,n,&i,&j,&k);
    X[0] = point_value(i,&coll_s[0]);/* r */
    X[1] = point_value(j,&coll_s[1]);/* theta */
    X[2] = point_value(k,&coll_s[2]);/* phi */
    patch->node[l]->X = X;
    
    R1 = R1_field->v[i_j_k_to_ijk(n,0,j,k)];
    R2 = R2_field->v[i_j_k_to_ijk(n,0,j,k)];
    r = X[0]*(R2-R1)+R1;
    
    x[0] = r*sin(X[1])*cos(X[2]);
    x[1] = r*sin(X[1])*sin(X[2]);
    x[2] = r*cos(X[1]);
    
    x[0]+= c[0];
    x[1]+= c[1];
    x[2]+= c[2];
    
  }
}

/* filling patch struct for BNS_Spherical_grid */
void fill_patches_BNS_Spherical_grid(Grid_T *const grid)
{
  const Uint N_outermost_split = (Uint) Pgeti("Number_of_Outermost_Split");
  Uint pn,i;
  
  pn = 0;
  populate_left_NS_sphere(grid,pn++);
  populate_left_NS_around_sphere(grid,pn++);
  for (i = 0; i < N_outermost_split; i++)
    Error0(NO_JOB);
 
  populate_right_NS_sphere(grid,pn++);
  populate_right_NS_around_sphere(grid,pn++);
  for (i = 0; i < N_outermost_split; i++)
    Error0(NO_JOB);

  
}

/* populating properties of patch for left NS */
static void populate_left_NS_sphere(Grid_T *const grid,const Uint pn)
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
  Uint n,j,k;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,PATCH_NAME_PRT_P_"left_NS",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (Uint)PgetiEZ("n_a");
  patch->n[1] = (Uint)PgetiEZ("n_b");
  patch->n[2] = (Uint)PgetiEZ("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (Uint)PgetiEZ(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (Uint)PgetiEZ(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (Uint)PgetiEZ(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    Error0("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_a",grid->gn);
  patch->c[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_b",grid->gn);
  patch->c[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_c",grid->gn);
  patch->c[2] = Pgetd(var);
  
  /* filling Rs */
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_R1",grid->gn);
  R1_array = Pgetdd(var);
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_R2",grid->gn);
  R2_array = Pgetdd(var);
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (j = 0; j < patch->n[1]; ++j)
    for (k = 0; k < patch->n[2]; ++k)
    {
      Uint ijk = i_j_k_to_ijk(patch->n,0,j,k);
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
  patch->collocation[0] = get_collocation(Pgets("collocation_a"));
  patch->collocation[1] = get_collocation(Pgets("collocation_b"));
  patch->collocation[2] = get_collocation(Pgets("collocation_c"));

  /* check for override */
  make_keyword_parameter(&ret,name,"collocation");
  c = get_collocation(PgetsEZ(ret.s0));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[0] = c;
  c = get_collocation(PgetsEZ(ret.s1));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[1] = c;
  c = get_collocation(PgetsEZ(ret.s2));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[2] = c;
  
  assert(patch->collocation[0] != UNDEFINED_COLLOCATION);
  assert(patch->collocation[1] != UNDEFINED_COLLOCATION);
  assert(patch->collocation[2] != UNDEFINED_COLLOCATION);
  
  /* basis */
  patch->basis[0] = get_basis(Pgets("basis_a"));
  patch->basis[1] = get_basis(Pgets("basis_b"));
  patch->basis[2] = get_basis(Pgets("basis_c"));

  /* check for override */
  make_keyword_parameter(&ret,name,"basis");
  b = get_basis(PgetsEZ(ret.s0));
  if ( b != UNDEFINED_BASIS)
    patch->basis[0] = b;
  b = get_basis(PgetsEZ(ret.s1));
  if ( b != UNDEFINED_BASIS)
    patch->basis[1] = b;
  b = get_basis(PgetsEZ(ret.s2));
  if ( b != UNDEFINED_BASIS)
    patch->basis[2] = b;
  
  assert(patch->basis[0] != UNDEFINED_BASIS);
  assert(patch->basis[1] != UNDEFINED_BASIS);
  assert(patch->basis[2] != UNDEFINED_BASIS);
  
}

/* populating properties of patch for right NS */
static void populate_right_NS_sphere(Grid_T *const grid,const Uint pn)
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
  Uint n,k,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,PATCH_NAME_PRT_P_"right_NS",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (Uint)PgetiEZ("n_a");
  patch->n[1] = (Uint)PgetiEZ("n_b");
  patch->n[2] = (Uint)PgetiEZ("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (Uint)PgetiEZ(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (Uint)PgetiEZ(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (Uint)PgetiEZ(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    Error0("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_a",grid->gn);
  patch->c[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_b",grid->gn);
  patch->c[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_c",grid->gn);
  patch->c[2] = Pgetd(var);
  
  /* filling Rs */
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_R1",grid->gn);
  R1_array = Pgetdd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_R2",grid->gn);
  R2_array = Pgetdd(var);
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (j = 0; j < patch->n[1]; ++j)
    for (k = 0; k < patch->n[2]; ++k)
    {
      Uint ijk = i_j_k_to_ijk(patch->n,0,j,k);
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
  patch->collocation[0] = get_collocation(Pgets("collocation_a"));
  patch->collocation[1] = get_collocation(Pgets("collocation_b"));
  patch->collocation[2] = get_collocation(Pgets("collocation_c"));

  /* check for override */
  make_keyword_parameter(&ret,name,"collocation");
  c = get_collocation(PgetsEZ(ret.s0));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[0] = c;
  c = get_collocation(PgetsEZ(ret.s1));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[1] = c;
  c = get_collocation(PgetsEZ(ret.s2));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[2] = c;
  
  assert(patch->collocation[0] != UNDEFINED_COLLOCATION);
  assert(patch->collocation[1] != UNDEFINED_COLLOCATION);
  assert(patch->collocation[2] != UNDEFINED_COLLOCATION);
  
  /* basis */
  patch->basis[0] = get_basis(Pgets("basis_a"));
  patch->basis[1] = get_basis(Pgets("basis_b"));
  patch->basis[2] = get_basis(Pgets("basis_c"));

  /* check for override */
  make_keyword_parameter(&ret,name,"basis");
  b = get_basis(PgetsEZ(ret.s0));
  if ( b != UNDEFINED_BASIS)
    patch->basis[0] = b;
  b = get_basis(PgetsEZ(ret.s1));
  if ( b != UNDEFINED_BASIS)
    patch->basis[1] = b;
  b = get_basis(PgetsEZ(ret.s2));
  if ( b != UNDEFINED_BASIS)
    patch->basis[2] = b;
  
  assert(patch->basis[0] != UNDEFINED_BASIS);
  assert(patch->basis[1] != UNDEFINED_BASIS);
  assert(patch->basis[2] != UNDEFINED_BASIS);
    
}

/* populating properties of patch for left NS's around */
static void populate_left_NS_around_sphere(Grid_T *const grid,const Uint pn)
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
  Uint n,k,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,PATCH_NAME_PRT_P_"left_NS_around",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (Uint)PgetiEZ("n_a");
  patch->n[1] = (Uint)PgetiEZ("n_b");
  patch->n[2] = (Uint)PgetiEZ("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (Uint)PgetiEZ(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (Uint)PgetiEZ(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (Uint)PgetiEZ(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    Error0("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_a",grid->gn);
  patch->c[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_b",grid->gn);
  patch->c[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_c",grid->gn);
  patch->c[2] = Pgetd(var);
  
  /* filling Rs */
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_Surrounding_R2",grid->gn);
  R2_const = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_R2",grid->gn);
  R1_array = Pgetdd(var);
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (j = 0; j < patch->n[1]; ++j)
    for (k = 0; k < patch->n[2]; ++k)
    {
      Uint ijk = i_j_k_to_ijk(patch->n,0,j,k);
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
  patch->collocation[0] = get_collocation(Pgets("collocation_a"));
  patch->collocation[1] = get_collocation(Pgets("collocation_b"));
  patch->collocation[2] = get_collocation(Pgets("collocation_c"));

  /* check for override */
  make_keyword_parameter(&ret,name,"collocation");
  c = get_collocation(PgetsEZ(ret.s0));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[0] = c;
  c = get_collocation(PgetsEZ(ret.s1));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[1] = c;
  c = get_collocation(PgetsEZ(ret.s2));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[2] = c;
  
  assert(patch->collocation[0] != UNDEFINED_COLLOCATION);
  assert(patch->collocation[1] != UNDEFINED_COLLOCATION);
  assert(patch->collocation[2] != UNDEFINED_COLLOCATION);
  
  /* basis */
  patch->basis[0] = get_basis(Pgets("basis_a"));
  patch->basis[1] = get_basis(Pgets("basis_b"));
  patch->basis[2] = get_basis(Pgets("basis_c"));

  /* check for override */
  make_keyword_parameter(&ret,name,"basis");
  b = get_basis(PgetsEZ(ret.s0));
  if ( b != UNDEFINED_BASIS)
    patch->basis[0] = b;
  b = get_basis(PgetsEZ(ret.s1));
  if ( b != UNDEFINED_BASIS)
    patch->basis[1] = b;
  b = get_basis(PgetsEZ(ret.s2));
  if ( b != UNDEFINED_BASIS)
    patch->basis[2] = b;
  
  assert(patch->basis[0] != UNDEFINED_BASIS);
  assert(patch->basis[1] != UNDEFINED_BASIS);
  assert(patch->basis[2] != UNDEFINED_BASIS);
    
}

/* populating properties of patch for right NS's around */
static void populate_right_NS_around_sphere(Grid_T *const grid,const Uint pn)
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
  Uint n,k,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,PATCH_NAME_PRT_P_"right_NS_around",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (Uint)PgetiEZ("n_a");
  patch->n[1] = (Uint)PgetiEZ("n_b");
  patch->n[2] = (Uint)PgetiEZ("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (Uint)PgetiEZ(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (Uint)PgetiEZ(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (Uint)PgetiEZ(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    Error0("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_a",grid->gn);
  patch->c[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_b",grid->gn);
  patch->c[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_c",grid->gn);
  patch->c[2] = Pgetd(var);
  
  /* filling Rs */
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_Surrounding_R2",grid->gn);
  R2_const = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_R2",grid->gn);
  R1_array = Pgetdd(var);
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (j = 0; j < patch->n[1]; ++j)
    for (k = 0; k < patch->n[2]; ++k)
    {
      Uint ijk = i_j_k_to_ijk(patch->n,0,j,k);
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
  patch->collocation[0] = get_collocation(Pgets("collocation_a"));
  patch->collocation[1] = get_collocation(Pgets("collocation_b"));
  patch->collocation[2] = get_collocation(Pgets("collocation_c"));

  /* check for override */
  make_keyword_parameter(&ret,name,"collocation");
  c = get_collocation(PgetsEZ(ret.s0));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[0] = c;
  c = get_collocation(PgetsEZ(ret.s1));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[1] = c;
  c = get_collocation(PgetsEZ(ret.s2));
  if (c != UNDEFINED_COLLOCATION)
    patch->collocation[2] = c;
  
  assert(patch->collocation[0] != UNDEFINED_COLLOCATION);
  assert(patch->collocation[1] != UNDEFINED_COLLOCATION);
  assert(patch->collocation[2] != UNDEFINED_COLLOCATION);
  
  /* basis */
  patch->basis[0] = get_basis(Pgets("basis_a"));
  patch->basis[1] = get_basis(Pgets("basis_b"));
  patch->basis[2] = get_basis(Pgets("basis_c"));

  /* check for override */
  make_keyword_parameter(&ret,name,"basis");
  b = get_basis(PgetsEZ(ret.s0));
  if ( b != UNDEFINED_BASIS)
    patch->basis[0] = b;
  b = get_basis(PgetsEZ(ret.s1));
  if ( b != UNDEFINED_BASIS)
    patch->basis[1] = b;
  b = get_basis(PgetsEZ(ret.s2));
  if ( b != UNDEFINED_BASIS)
    patch->basis[2] = b;
  
  assert(patch->basis[0] != UNDEFINED_BASIS);
  assert(patch->basis[1] != UNDEFINED_BASIS);
  assert(patch->basis[2] != UNDEFINED_BASIS);
    
}

/* memory alloc patches for BNS_Spherical type */
void alloc_patches_BNS_Spherical_grid(Grid_T *const grid)
{
  Uint Np = 4;/* number of patches without outermost's*/
  Uint outermost;
  Uint i;
  
  outermost = (Uint) PgetiEZ("Number_of_Outermost_Split");
  if (outermost != (Uint)INT_MAX)
    Np += 2*outermost;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  IsNull(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    IsNull(grid->patch[i]);
  }
  
}
