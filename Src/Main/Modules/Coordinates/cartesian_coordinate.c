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
  patch->JacobianT->dN0_dx = dN0_dx_Cartesian_patch;
  patch->JacobianT->dN0_dy = dN0_dy_Cartesian_patch;
  patch->JacobianT->dN0_dz = dN0_dz_Cartesian_patch;
  patch->JacobianT->dN1_dx = dN1_dx_Cartesian_patch;
  patch->JacobianT->dN1_dy = dN1_dy_Cartesian_patch;
  patch->JacobianT->dN1_dz = dN1_dz_Cartesian_patch;
  patch->JacobianT->dN2_dx = dN2_dx_Cartesian_patch;
  patch->JacobianT->dN2_dy = dN2_dy_Cartesian_patch;
  patch->JacobianT->dN2_dz = dN2_dz_Cartesian_patch;
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

/* Calculating dN0/dx at arbitrary curvilinear point point X.
// used for interpolation.
// x = 0.5*(min-max)*N0 +0.5*(min+max)
// ->return value: dN0/dx */
double dN0_dx_Cartesian_patch(Patch_T *const patch,const double *const X)
{
  double dN0_dx = 0;
  
  if (patch->collocation[0] == Chebyshev_Extrema || 
      patch->collocation[0] == Chebyshev_Nodes)
      dN0_dx = 2./(-patch->max[0]+patch->min[0]);
  else
    abortEr(NO_JOB);
 
 UNUSED(X);
 
 return dN0_dx;   
}

/* Calculating dN0/dy at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN0/dy */
double dN0_dy_Cartesian_patch(Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN0/dz at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN0/dz */
double dN0_dz_Cartesian_patch(Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN1/dx at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN1/dx */
double dN1_dx_Cartesian_patch(Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN1/dy at arbitrary curvilinear point point X.
// used for interpolation.
// y = 0.5*(min-max)*N1 +0.5*(min+max)
// ->return value: dN1/dy */
double dN1_dy_Cartesian_patch(Patch_T *const patch,const double *const X)
{
  double dN1_dy = 0;
  
  if (patch->collocation[1] == Chebyshev_Extrema || 
      patch->collocation[1] == Chebyshev_Nodes)
      dN1_dy = 2./(-patch->max[1]+patch->min[1]);
  else
    abortEr(NO_JOB);
    
  UNUSED(X);
  
  return dN1_dy;
}
/* Calculating dN1/dz at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN1/dz */
double dN1_dz_Cartesian_patch(Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN0/dx at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN0/dx */
double dN2_dx_Cartesian_patch(Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN0/dx at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN0/dx */
double dN2_dy_Cartesian_patch(Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN0/dx at arbitrary curvilinear point point X.
// used for interpolation.
// z = 0.5*(min-max)*N2 +0.5*(min+max)
// ->return value: dN0/dx */
double dN2_dz_Cartesian_patch(Patch_T *const patch,const double *const X)
{
  double dN2_dz = 0;
  
  if (patch->collocation[2] == Chebyshev_Extrema || 
      patch->collocation[2] == Chebyshev_Nodes)
   dN2_dz = 2./(-patch->max[2]+patch->min[2]);
  else
    abortEr(NO_JOB);
    
  UNUSED(X);
  
  return dN2_dz;
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
