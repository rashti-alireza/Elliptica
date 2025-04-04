/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2021 Alireza Rashti.
*/


#include "eq_header.h"
#include "maths_equation_solvings_lib.h"

void *eq_XCTS_curve_T3_ddm_bc_phi(void *vp1,void *vp2);
void *eq_XCTS_curve_T3_ddm_bc_phi(void *vp1,void *vp2)
{
  DDM_SCHUR_BC_DECLARE
  Uint ijk;/* node index */

  /* declaring: */
  READ_v(dphi_D2)
  READ_v(dphi_D1)
  READ_v(dphi_D0)
  READ_v(W_U1)
  READ_v(W_U0)
  READ_v(W_U2)
  READ_v(u0)
  READ_v(psi)
  READ_v(beta_U1)
  READ_v(beta_U0)
  READ_v(beta_U2)
  READ_v(gConf_D0D2)
  READ_v(gConf_D0D0)
  READ_v(gConf_D0D1)
  READ_v(gConf_D1D2)
  READ_v(gConf_D1D1)
  READ_v(gConf_D2D2)


  if (patch->outerB)/* at outer boundary */
  {
  DDM_SCHUR_BC_OPEN

 double n_U0 = dq2_dq1(patch,_c_,_x_,ijk);
 double n_U1 = dq2_dq1(patch,_c_,_y_,ijk);
 double n_U2 = dq2_dq1(patch,_c_,_z_,ijk);
 double Norm = 
gConf_D0D0[ijk]*pow(n_U0, 2) + 2.0*gConf_D0D1[ijk]*n_U0*n_U1 + 2.0*
gConf_D0D2[ijk]*n_U0*n_U2 + gConf_D1D1[ijk]*pow(n_U1, 2) + 2.0*
gConf_D1D2[ijk]*n_U1*n_U2 + gConf_D2D2[ijk]*pow(n_U2, 2);

 double s_U1 = 
pow(Norm, -0.5)*n_U1;

 double s_U0 = 
pow(Norm, -0.5)*n_U0;

 double s_U2 = 
pow(Norm, -0.5)*n_U2;

 double psi4 = 
pow(psi[ijk], 4);

 double outerB_F = 
dphi_D0[ijk]*s_U0 + dphi_D1[ijk]*s_U1 + dphi_D2[ijk]*s_U2 + psi4*
(gConf_D0D0[ijk]*s_U0*(W_U0[ijk] - beta_U0[ijk]*u0[ijk]) +
gConf_D0D1[ijk]*s_U0*(W_U1[ijk] - beta_U1[ijk]*u0[ijk]) +
gConf_D0D1[ijk]*s_U1*(W_U0[ijk] - beta_U0[ijk]*u0[ijk]) +
gConf_D0D2[ijk]*s_U0*(W_U2[ijk] - beta_U2[ijk]*u0[ijk]) +
gConf_D0D2[ijk]*s_U2*(W_U0[ijk] - beta_U0[ijk]*u0[ijk]) +
gConf_D1D1[ijk]*s_U1*(W_U1[ijk] - beta_U1[ijk]*u0[ijk]) +
gConf_D1D2[ijk]*s_U1*(W_U2[ijk] - beta_U2[ijk]*u0[ijk]) +
gConf_D1D2[ijk]*s_U2*(W_U1[ijk] - beta_U1[ijk]*u0[ijk]) +
gConf_D2D2[ijk]*s_U2*(W_U2[ijk] - beta_U2[ijk]*u0[ijk]));

  schur_F[schur_ijk] = outerB_F;

  DDM_SCHUR_BC_CLOSE
  }/* end of if (patch->outerB) */
  else if (patch->innerB)/* at inner boundary */
  {
  DDM_SCHUR_BC_OPEN

  double innerB_F = 
0;

  schur_F[schur_ijk] = innerB_F;

  DDM_SCHUR_BC_CLOSE
  }/* end of else if (patch->innerB) */
  return 0;
}
