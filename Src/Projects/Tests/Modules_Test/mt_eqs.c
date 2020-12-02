/*
// Alireza Rashti
// August 2018
*/

#include "mt_eqs.h"

/* adding mt equations to data base */
void mt_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq)
{
  /* adding field and boundary condition equations */
  *field_eq      = init_eq();
  *bc_eq         = init_eq();
  *jacobian_field_eq   = init_eq();
  *jacobian_bc_eq      = init_eq();

  add_eq(field_eq,eq_alpha,"eq_alpha");
  add_eq(bc_eq,bc_alpha,"bc_alpha");
  add_eq(jacobian_field_eq,jacobian_eq_alpha,"jacobian_eq_alpha");
  add_eq(jacobian_bc_eq   ,jacobian_bc_alpha,"jacobian_bc_alpha");
}

/* initial data for field alpha
// ->return value: EXIT_SUCCESS
*/
int mt_initial_data_alpha(Grid_T *const grid)
{
  Uint p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    double *alpha = patch->fields[Ind("alpha")]->v;
    Uint n;
    
    FOR_ALL_POINTS(n,patch)
      alpha[n] = Pow2(x_(n))+Pow2(y_(n))+Pow2(z_(n))+0.3+sin(x_(n));
  }
  
  return EXIT_SUCCESS;
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
  Field_T *const alpha = patch->fields[Ind("alpha")];
  double *const F = S->f;
  const Uint *const node  = S->inv;/* inverse map to node */
  double *alpha_xx = Partial_Derivative(alpha,"x,x");
  double *alpha_yy = Partial_Derivative(alpha,"y,y");
  double *alpha_zz = Partial_Derivative(alpha,"z,z");
  const Uint N = S->Oi;/* number of inner mesh nodes */
  Uint ijk;/* node */
  Uint n;
  
  for (n = 0; n < N; ++n)
  {
    ijk  = node[n];
    F[n] = alpha_xx[ijk]+alpha_yy[ijk]+alpha_zz[ijk]-6;/* note: for each new equation, only this line is changed */
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
  Uint *const map  = S->map;
  Patch_T *const patch = bc->patch;
  const Uint *const node = bc->node;/* nodes at boundary */
  const Uint N = bc->nn;/* number of nodes at boundary */
  double *const alpha = patch->fields[Ind("alpha")]->v;
  Uint n,ijk;
  
  /* alpha at outer boundary */
  for (n = 0; n < N; ++n)
  {
    ijk = node[n];
    F[map[ijk]] = alpha[ijk]-Pow2(x_(ijk))-Pow2(y_(ijk))-Pow2(z_(ijk));/* note: for each new equation, only this line is changed */
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
  const Uint *const node = S->inv;
  const Uint Ni = S->Oi;/* number of inner mesh nodes */
  const Uint Nj = S->NS;/* number of inner mesh+outer-boundary nodes*/
  const Uint K0 = S->NS;/* number of inner mesh+outer-boundary nodes*/
  const Uint Nk = patch->nn;/* total number of nodes */
  const Uint Ref = Nj;/* for shorhand purposes */
  const char *types[] = {"dfxx_df","dfyy_df","dfzz_df",0};
  fJs_T *dfxx_df = 0,*dfyy_df = 0,*dfzz_df = 0;
  Matrix_T *j0 = 0,*j1 = 0,*j2 = 0;
  Uint i,j,k,ijk,lmn;
  
  prepare_Js_jacobian_eq(patch,types);
  j0      = get_j_matrix(patch,"dfxx_df");
  j1      = get_j_matrix(patch,"dfyy_df");
  j2      = get_j_matrix(patch,"dfzz_df");
  dfxx_df = get_j_reader(j0);
  dfyy_df = get_j_reader(j1);
  dfzz_df = get_j_reader(j2);
  
  /* fill up jacobian for alpha equation: */
  
  /* B part: */
  for (i = 0; i < Ni; ++i)
  {
    ijk = node[i];
    for (j = 0; j < Nj; ++j)
    {
      lmn = node[j];
      B[i][j] = dfxx_df(j0,ijk,lmn)+dfyy_df(j1,ijk,lmn)+dfzz_df(j2,ijk,lmn);/* note: for each new equation, only this line is changed */
    }
  }
  
  /* E part: */
  if (S->NI)/* if there is any interface points then E is needed */
  {
    E_Trans = S->E_Trans->reg->A;
    
    for (k = K0; k < Nk; ++k)
    {
      lmn = node[k];
      j = k-Ref;
      for (i = 0; i < Ni; ++i)
      {
        ijk = node[i];
        E_Trans[j][i] = dfxx_df(j0,ijk,lmn)+dfyy_df(j1,ijk,lmn)+dfzz_df(j2,ijk,lmn);/* note: for each new equation, only this line is changed */
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
  const Uint I0 = S->Oi;/* number of inner mesh nodes */
  const Uint N = S->NS;/* number of inner mesh+outer-boundary nodes*/
  Uint i;
  
  /* fill up jacobian for alpha equation: */
  
  /* B part: */
  for (i = I0; i < N; ++i)
   B[i][i] = 1;/* note: for each new equation, only this line is changed, unless the b.c. is more complicated and has field in it */
                
  /* E part: */
  /* since B.C. equation doesn't have any dependency on interface points,
  // E part is zero so nothing needs to be done. */
  
  UNUSED(vp1);
  
  return 0;
}


