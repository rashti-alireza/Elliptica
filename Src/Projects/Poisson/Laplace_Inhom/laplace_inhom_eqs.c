/*
// Alireza Rashti
// August 2018
*/

#include "laplace_inhom_eqs.h"

/* adding Laplace Inhomogeneous equations to data base */
void Laplace_Inhom_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq)
{
  /* adding field and boundary condition equations */
  *field_eq      = init_eq();
  *bc_eq         = init_eq();
  *jacobian_field_eq   = init_eq();
  *jacobian_bc_eq      = init_eq();
  
  add_eq(field_eq,eq_alpha,"eq_alpha");
  add_eq(bc_eq   ,bc_alpha,"bc_alpha");
  add_eq(jacobian_field_eq,jacobian_eq_alpha,"jacobian_eq_alpha");
  add_eq(jacobian_bc_eq   ,jacobian_bc_alpha,"jacobian_bc_alpha");
}

/* laplace equation for alpha
// \nabla^2(alpha) = 6
// note: void *Structure is some general structure capable of embracing
// different types and situations for equations.
*/
static void *eq_alpha(void *vp1,void *vp2)
{
  Patch_T *const patch = vp1;
  DDM_Schur_Complement_T *const S = vp2;
  Field_T *const alpha = patch->pool[Ind("alpha")];
  double *const F = S->f;
  const unsigned *const node  = S->inv;/* inverse map to node */
  double *alpha_xx = Partial_Derivative(alpha,"x,x");
  double *alpha_yy = Partial_Derivative(alpha,"y,y");
  double *alpha_zz = Partial_Derivative(alpha,"z,z");
  const unsigned N = S->Oi;/* number of inner mesh nodes */
  unsigned ijk;/* node */
  unsigned n;
  
  for (n = 0; n < N; ++n)
  {
    ijk  = node[n];
    F[n] = alpha_xx[ijk]+alpha_yy[ijk]+alpha_zz[ijk]-6;
  }
    
  free(alpha_xx);
  free(alpha_yy);
  free(alpha_zz);
  
  return 0;
}

/* boundary condition for alpha 
// alpha = (x^2+y^2+z^2) at boundary
*/
static void *bc_alpha(void *vp1,void *vp2)
{
  Boundary_Condition_T *const bc = vp1;
  DDM_Schur_Complement_T *const S = vp2;
  double *const F      = S->f;
  unsigned *const map  = S->map;
  Patch_T *const patch = bc->patch;
  const unsigned *const node = bc->node;/* nodes at boundary */
  const unsigned N = bc->nn;/* number of nodes at boundary */
  double *const alpha = patch->pool[Ind("alpha")]->v;
  unsigned n,ijk;
  
  /* alpha at outer boundary */
  for (n = 0; n < N; ++n)
  {
    ijk = node[n];
    F[map[ijk]] = alpha[ijk]-SQR(x_(ijk))-SQR(y_(ijk))-SQR(z_(ijk));
  }
      
  return 0;
}

/* jacobian equation for alpha.
// d(nabla^2 alpha(i))/d(alpha(j))
*/
static void *jacobian_eq_alpha(void *vp1,void *vp2)
{
  Patch_T *const patch  = vp1;
  DDM_Schur_Complement_T *const S = vp2;
  double **const B = S->B->reg->A;
  double **E_Trans;
  const unsigned *const inv = S->inv;
  const unsigned NInnerMesh = S->Oi;/* number of inner mesh nodes */
  const unsigned NInnerMeshAndOuterB = S->NS;/* number of inner mesh+outer-boundary nodes*/
  const unsigned NNode = patch->nn;/* total number of nodes */
  const unsigned Ref = NInnerMeshAndOuterB;/* for shorhand purposes */
  const char *types[] = {"j_xx","j_yy","j_zz",0};
  fJs_T *j_xx = 0,*j_yy = 0,*j_zz = 0;
  Matrix_T *j0 = 0,*j1 = 0,*j2 = 0;
  unsigned i,j,ii,jj;
  
  prepare_Js_jacobian_eq(patch,types);
  j0   = get_j_matrix(patch,"j_xx");
  j1   = get_j_matrix(patch,"j_yy");
  j2   = get_j_matrix(patch,"j_zz");
  j_xx = get_j_reader(j0);
  j_yy = get_j_reader(j1);
  j_zz = get_j_reader(j2);
  
  /* fill up jacobian for alpha equation: */
  //TIMER_ON(Jacobian_B);
  /* B part: */
  for (i = 0; i < NInnerMesh; ++i)
  {
    ii = inv[i];
    for (j = 0; j < NInnerMeshAndOuterB; ++j)
    {
      jj = inv[j];
      B[i][j] = j_xx(j0,ii,jj)+j_yy(j1,ii,jj)+j_zz(j2,ii,jj);
    }
  }
  //TIMER_OFF(Jacobian_B);
  /* E part: */
  if (S->NI)/* if there is any interface points then E is needed */
  {
    E_Trans = S->E_Trans->reg->A;
    for (j = NInnerMeshAndOuterB; j < NNode; ++j)
    {
      jj = inv[j];
      for (i = 0; i < NInnerMesh; ++i)
      {
        ii = inv[i];
        E_Trans[j-Ref][i] = j_xx(j0,ii,jj)+j_yy(j1,ii,jj)+j_zz(j2,ii,jj);
      }
    }
  }
  
  return 0;
}

/* jacobian alpha B.C. 
// alpha - (x^2+y^2+z^2) = 0
*/
static void *jacobian_bc_alpha(void *vp1,void *vp2)
{
  DDM_Schur_Complement_T *const S = vp2;
  double **const B = S->B->reg->A;
  const unsigned NInnerMesh = S->Oi;/* number of inner mesh nodes */
  const unsigned NInnerMeshAndOuterB = S->NS;/* number of inner mesh+outer-boundary nodes*/
  unsigned i;
  
  /* fill up jacobian for alpha equation: */
  
  /* B part: */
  for (i = NInnerMesh; i < NInnerMeshAndOuterB; ++i)
   B[i][i] = 1;
                
  /* E part: */
  /* since B.C. equation doesn't have any dependency on interface points,
  // E part is zero so nothing needs to be done. */
  
  UNUSED(vp1);
  
  return 0;
}


