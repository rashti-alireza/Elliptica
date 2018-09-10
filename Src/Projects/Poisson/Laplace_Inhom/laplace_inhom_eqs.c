/*
// Alireza Rashti
// August 2018
*/

#include "laplace_inhom_eqs.h"

/* adding Laplace Inhomogeneous equations to data base */
void Laplace_Inhom_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_eq)
{
  /* adding field and boundary condition equations */
  *field_eq      = init_eq();
  *bc_eq         = init_eq();
  *jacobian_eq   = init_eq();
  
  add_eq(field_eq,eq_alpha,"eq_alpha");
  add_eq(bc_eq,bc_alpha,"bc_alpha");
  add_eq(jacobian_eq,jacobian_alpha,"jacobian_alpha");
}

/* laplace equation for alpha
// \nabla^2(alpha) = 6
// note: void *Structure is some general structure capable of embracing
// different types and situations for equations.
*/
static void *eq_alpha(void *vp1,void *vp2)
{
  Field_T *const alpha = vp1;
  double *const F      = vp2;
  double *alpha_xx = Partial_Derivative(alpha,"x,x");
  double *alpha_yy = Partial_Derivative(alpha,"y,y");
  double *alpha_zz = Partial_Derivative(alpha,"z,z");
  Patch_T *const patch = alpha->patch;
  const unsigned nn = patch->nn;
  unsigned i;
  
  for(i = 0; i < nn; ++i)
    F[i] = 
      alpha_xx[i]+alpha_yy[i]+alpha_zz[i]-6;//1/(1+pow(SQR(x_(i)),100)+pow(SQR(y_(i)),100)+pow(SQR(z_(i)),100);
    
  free(alpha_xx);
  free(alpha_yy);
  free(alpha_zz);
  return 0;
}

/* boundary condition for alpha */
static void *bc_alpha(void *vp1,void *vp2)
{
  Boundary_Condition_T *const bc = vp1;
  double *const F = vp2;
  Patch_T *const patch = bc->subface->patch;
  const unsigned *const Bnode = bc->node;/* nodes at boundary */
  const unsigned nn = bc->nn;/* number of nodes at boundary */
  double *const alpha = bc->field->v;
  unsigned i;
  
  /* alpha = 0 at outer boundary */
  for (i = 0; i < nn; ++i)
    F[Bnode[i]] = alpha[Bnode[i]]
      -SQR(x_(Bnode[i]))-SQR(y_(Bnode[i]))-SQR(z_(Bnode[i]));
      
  return 0;
}

/* jacobian equation for alpha.
// d(nabla^2 alpha(i))/d(alpha(j))
*/
static void *jacobian_alpha(void *vp1,void *vp2)
{
  Jacobian_Eq_T *const jac = vp1;
  Patch_T *const patch     = vp2;
  const unsigned nn = patch->nn;
  const char *types[] = {"j_xx","j_yy","j_zz",0};
  fJs_T *j_xx = 0,*j_yy = 0,*j_zz = 0;
  Matrix_T *j0 = 0,*j1 = 0,*j2 = 0;
  Matrix_T *alpha_jac;
  double **J;
  unsigned i;
  
  prepare_Js_jacobian_eq(patch,types);
  j0   = get_j_matrix(patch,"j_xx");
  j1   = get_j_matrix(patch,"j_yy");
  j2   = get_j_matrix(patch,"j_zz");
  j_xx = get_j_reader(j0);
  j_yy = get_j_reader(j1);
  j_zz = get_j_reader(j2);
  
  alpha_jac = alloc_matrix(REG_SF,nn,nn);
  J = alpha_jac->reg->A;
  
  /* fill up jacobian for alpha equation */
  for (i = 0; i < nn; ++i)
  {
    unsigned j;
    for (j = 0; j < nn; ++j)
      J[i][j] = j_xx(j0,i,j)+j_yy(j1,i,j)+j_zz(j2,i,j);
  }
  
  UNUSED(jac);
  
  return alpha_jac;
}

