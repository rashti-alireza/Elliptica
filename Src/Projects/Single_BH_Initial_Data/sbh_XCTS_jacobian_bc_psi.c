/*
  These C codes generated by Cpi version 1.0
  Copyright (C) 2019 Alireza Rashti.
*/


#include "sbh_headers.h"
#include "maths_equation_solvings_lib.h"
#include "sbh_XCTS_equations_lib.h"


void *sbh_jacobian_bc_psi(void *vp1,void *vp2)
{
  DDM_SCHUR_JACOBIAN_BC_DECLARE
  unsigned ijk,lmn;/* for Jacobian entries J[ijk][lmn] */
  const double kd[2] = {0.,1.};/* Kronecker delta */

  /* declaring: */
  GET_FIELD(psi)
  GET_FIELD(K)
  GET_FIELD(_A_UiUj_U2U2)
  GET_FIELD(_A_UiUj_U1U2)
  GET_FIELD(_A_UiUj_U1U1)
  GET_FIELD(_A_UiUj_U0U2)
  GET_FIELD(_A_UiUj_U0U1)
  GET_FIELD(_A_UiUj_U0U0)
  GET_FIELD(dpsi_D0)
  GET_FIELD(dpsi_D1)
  GET_FIELD(dpsi_D2)
  GET_FIELD(_gamma_D2D2)
  GET_FIELD(_gamma_D0D2)
  GET_FIELD(_gamma_D0D0)
  GET_FIELD(_gamma_D0D1)
  GET_FIELD(_gamma_D1D2)
  GET_FIELD(_gamma_D1D1)
  GET_FIELD_IF_ON_HORIZON(_HS_U2)
  GET_FIELD_IF_ON_HORIZON(_HS_U0)
  GET_FIELD_IF_ON_HORIZON(_HS_U1)
  JACOBIAN_DERIVATIVE(Jpsi_D2)
  JACOBIAN_DERIVATIVE(Jpsi_D0)
  JACOBIAN_DERIVATIVE(Jpsi_D1)


  if (patch->outerB)
  {
  DDM_SCHUR_JACOBIAN_BC_Bpart_OPEN

  double outerB_Bpart = 
kd[ijk==lmn];

  B[i][j] = outerB_Bpart;

  DDM_SCHUR_JACOBIAN_BC_Bpart_CLOSE

  DDM_SCHUR_JACOBIAN_BC_Epart_OPEN

  double outerB_Epart = 
0;

  E_Trans[j][i] = outerB_Epart;

  DDM_SCHUR_JACOBIAN_BC_Epart_CLOSE
  }/* end of if (patch->outerB) */
  else if (patch->innerB)
  {
  DDM_SCHUR_JACOBIAN_BC_Bpart_OPEN

  double B_t1 = 
(_A_UiUj_U0U0[ijk]*pow(_HS_U0[ijk], 2)*pow(_gamma_D0D0[ijk], 2) + 2.0*
_A_UiUj_U0U0[ijk]*_HS_U0[ijk]*_HS_U1[ijk]*_gamma_D0D0[ijk]*
_gamma_D0D1[ijk] + 2.0*_A_UiUj_U0U0[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*
_gamma_D0D0[ijk]*_gamma_D0D2[ijk] + _A_UiUj_U0U0[ijk]*
pow(_HS_U1[ijk], 2)*pow(_gamma_D0D1[ijk], 2) + 2.0*_A_UiUj_U0U0[ijk]*
_HS_U1[ijk]*_HS_U2[ijk]*_gamma_D0D1[ijk]*_gamma_D0D2[ijk] +
_A_UiUj_U0U0[ijk]*pow(_HS_U2[ijk], 2)*pow(_gamma_D0D2[ijk], 2) + 2.0*
_A_UiUj_U0U1[ijk]*pow(_HS_U0[ijk], 2)*_gamma_D0D0[ijk]*
_gamma_D0D1[ijk] + 2.0*_A_UiUj_U0U1[ijk]*_HS_U0[ijk]*_HS_U1[ijk]*
_gamma_D0D0[ijk]*_gamma_D1D1[ijk] + 2.0*_A_UiUj_U0U1[ijk]*_HS_U0[ijk]*
_HS_U1[ijk]*pow(_gamma_D0D1[ijk], 2) + 2.0*_A_UiUj_U0U1[ijk]*
_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D0[ijk]*_gamma_D1D2[ijk] + 2.0*
_A_UiUj_U0U1[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D1[ijk]*
_gamma_D0D2[ijk] + 2.0*_A_UiUj_U0U1[ijk]*pow(_HS_U1[ijk], 2)*
_gamma_D0D1[ijk]*_gamma_D1D1[ijk] + 2.0*_A_UiUj_U0U1[ijk]*_HS_U1[ijk]*
_HS_U2[ijk]*_gamma_D0D1[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UiUj_U0U1[ijk]*
_HS_U1[ijk]*_HS_U2[ijk]*_gamma_D0D2[ijk]*_gamma_D1D1[ijk] + 2.0*
_A_UiUj_U0U1[ijk]*pow(_HS_U2[ijk], 2)*_gamma_D0D2[ijk]*
_gamma_D1D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*pow(_HS_U0[ijk], 2)*
_gamma_D0D0[ijk]*_gamma_D0D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*_HS_U0[ijk]*
_HS_U1[ijk]*_gamma_D0D0[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*
_HS_U0[ijk]*_HS_U1[ijk]*_gamma_D0D1[ijk]*_gamma_D0D2[ijk] + 2.0*
_A_UiUj_U0U2[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D0[ijk]*
_gamma_D2D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*
pow(_gamma_D0D2[ijk], 2) + 2.0*_A_UiUj_U0U2[ijk]*pow(_HS_U1[ijk], 2)*
_gamma_D0D1[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*_HS_U1[ijk]*
_HS_U2[ijk]*_gamma_D0D1[ijk]*_gamma_D2D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*
_HS_U1[ijk]*_HS_U2[ijk]*_gamma_D0D2[ijk]*_gamma_D1D2[ijk] + 2.0*
_A_UiUj_U0U2[ijk]*pow(_HS_U2[ijk], 2)*_gamma_D0D2[ijk]*
_gamma_D2D2[ijk] + _A_UiUj_U1U1[ijk]*pow(_HS_U0[ijk], 2)*
pow(_gamma_D0D1[ijk], 2) + 2.0*_A_UiUj_U1U1[ijk]*_HS_U0[ijk]*
_HS_U1[ijk]*_gamma_D0D1[ijk]*_gamma_D1D1[ijk] + 2.0*_A_UiUj_U1U1[ijk]*
_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D1[ijk]*_gamma_D1D2[ijk] +
_A_UiUj_U1U1[ijk]*pow(_HS_U1[ijk], 2)*pow(_gamma_D1D1[ijk], 2) + 2.0*
_A_UiUj_U1U1[ijk]*_HS_U1[ijk]*_HS_U2[ijk]*_gamma_D1D1[ijk]*
_gamma_D1D2[ijk] + _A_UiUj_U1U1[ijk]*pow(_HS_U2[ijk], 2)*
pow(_gamma_D1D2[ijk], 2) + 2.0*_A_UiUj_U1U2[ijk]*pow(_HS_U0[ijk], 2)*
_gamma_D0D1[ijk]*_gamma_D0D2[ijk] + 2.0*_A_UiUj_U1U2[ijk]*_HS_U0[ijk]*
_HS_U1[ijk]*_gamma_D0D1[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UiUj_U1U2[ijk]*
_HS_U0[ijk]*_HS_U1[ijk]*_gamma_D0D2[ijk]*_gamma_D1D1[ijk] + 2.0*
_A_UiUj_U1U2[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D1[ijk]*
_gamma_D2D2[ijk] + 2.0*_A_UiUj_U1U2[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*
_gamma_D0D2[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UiUj_U1U2[ijk]*
pow(_HS_U1[ijk], 2)*_gamma_D1D1[ijk]*_gamma_D1D2[ijk] + 2.0*
_A_UiUj_U1U2[ijk]*_HS_U1[ijk]*_HS_U2[ijk]*_gamma_D1D1[ijk]*
_gamma_D2D2[ijk] + 2.0*_A_UiUj_U1U2[ijk]*_HS_U1[ijk]*_HS_U2[ijk]*
pow(_gamma_D1D2[ijk], 2) + 2.0*_A_UiUj_U1U2[ijk]*pow(_HS_U2[ijk], 2)*
_gamma_D1D2[ijk]*_gamma_D2D2[ijk] + _A_UiUj_U2U2[ijk]*
pow(_HS_U0[ijk], 2)*pow(_gamma_D0D2[ijk], 2) + 2.0*_A_UiUj_U2U2[ijk]*
_HS_U0[ijk]*_HS_U1[ijk]*_gamma_D0D2[ijk]*_gamma_D1D2[ijk] + 2.0*
_A_UiUj_U2U2[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D2[ijk]*
_gamma_D2D2[ijk] + _A_UiUj_U2U2[ijk]*pow(_HS_U1[ijk], 2)*
pow(_gamma_D1D2[ijk], 2) + 2.0*_A_UiUj_U2U2[ijk]*_HS_U1[ijk]*
_HS_U2[ijk]*_gamma_D1D2[ijk]*_gamma_D2D2[ijk] + _A_UiUj_U2U2[ijk]*
pow(_HS_U2[ijk], 2)*pow(_gamma_D2D2[ijk], 2))/pow(psi[ijk], 5);

  double innerB_Bpart = 
(-1.0/3.0*kd[ijk==lmn]*pow(psi[ijk], 2)*(3*B_t1 + K[ijk]*psi[ijk]) -
kd[ijk==lmn]*(_HS_U0[ijk]*dpsi_D0[ijk] + _HS_U1[ijk]*dpsi_D1[ijk] +
_HS_U2[ijk]*dpsi_D2[ijk]) + psi[ijk]*(Jpsi_D0(j_Jpsi_D0,ijk,lmn)*
_HS_U0[ijk] + Jpsi_D1(j_Jpsi_D1,ijk,lmn)*_HS_U1[ijk] +
Jpsi_D2(j_Jpsi_D2,ijk,lmn)*_HS_U2[ijk]))/pow(psi[ijk], 2);

  B[i][j] = innerB_Bpart;

  DDM_SCHUR_JACOBIAN_BC_Bpart_CLOSE

  DDM_SCHUR_JACOBIAN_BC_Epart_OPEN

  double E_t1 = 
(_A_UiUj_U0U0[ijk]*pow(_HS_U0[ijk], 2)*pow(_gamma_D0D0[ijk], 2) + 2.0*
_A_UiUj_U0U0[ijk]*_HS_U0[ijk]*_HS_U1[ijk]*_gamma_D0D0[ijk]*
_gamma_D0D1[ijk] + 2.0*_A_UiUj_U0U0[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*
_gamma_D0D0[ijk]*_gamma_D0D2[ijk] + _A_UiUj_U0U0[ijk]*
pow(_HS_U1[ijk], 2)*pow(_gamma_D0D1[ijk], 2) + 2.0*_A_UiUj_U0U0[ijk]*
_HS_U1[ijk]*_HS_U2[ijk]*_gamma_D0D1[ijk]*_gamma_D0D2[ijk] +
_A_UiUj_U0U0[ijk]*pow(_HS_U2[ijk], 2)*pow(_gamma_D0D2[ijk], 2) + 2.0*
_A_UiUj_U0U1[ijk]*pow(_HS_U0[ijk], 2)*_gamma_D0D0[ijk]*
_gamma_D0D1[ijk] + 2.0*_A_UiUj_U0U1[ijk]*_HS_U0[ijk]*_HS_U1[ijk]*
_gamma_D0D0[ijk]*_gamma_D1D1[ijk] + 2.0*_A_UiUj_U0U1[ijk]*_HS_U0[ijk]*
_HS_U1[ijk]*pow(_gamma_D0D1[ijk], 2) + 2.0*_A_UiUj_U0U1[ijk]*
_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D0[ijk]*_gamma_D1D2[ijk] + 2.0*
_A_UiUj_U0U1[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D1[ijk]*
_gamma_D0D2[ijk] + 2.0*_A_UiUj_U0U1[ijk]*pow(_HS_U1[ijk], 2)*
_gamma_D0D1[ijk]*_gamma_D1D1[ijk] + 2.0*_A_UiUj_U0U1[ijk]*_HS_U1[ijk]*
_HS_U2[ijk]*_gamma_D0D1[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UiUj_U0U1[ijk]*
_HS_U1[ijk]*_HS_U2[ijk]*_gamma_D0D2[ijk]*_gamma_D1D1[ijk] + 2.0*
_A_UiUj_U0U1[ijk]*pow(_HS_U2[ijk], 2)*_gamma_D0D2[ijk]*
_gamma_D1D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*pow(_HS_U0[ijk], 2)*
_gamma_D0D0[ijk]*_gamma_D0D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*_HS_U0[ijk]*
_HS_U1[ijk]*_gamma_D0D0[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*
_HS_U0[ijk]*_HS_U1[ijk]*_gamma_D0D1[ijk]*_gamma_D0D2[ijk] + 2.0*
_A_UiUj_U0U2[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D0[ijk]*
_gamma_D2D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*
pow(_gamma_D0D2[ijk], 2) + 2.0*_A_UiUj_U0U2[ijk]*pow(_HS_U1[ijk], 2)*
_gamma_D0D1[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*_HS_U1[ijk]*
_HS_U2[ijk]*_gamma_D0D1[ijk]*_gamma_D2D2[ijk] + 2.0*_A_UiUj_U0U2[ijk]*
_HS_U1[ijk]*_HS_U2[ijk]*_gamma_D0D2[ijk]*_gamma_D1D2[ijk] + 2.0*
_A_UiUj_U0U2[ijk]*pow(_HS_U2[ijk], 2)*_gamma_D0D2[ijk]*
_gamma_D2D2[ijk] + _A_UiUj_U1U1[ijk]*pow(_HS_U0[ijk], 2)*
pow(_gamma_D0D1[ijk], 2) + 2.0*_A_UiUj_U1U1[ijk]*_HS_U0[ijk]*
_HS_U1[ijk]*_gamma_D0D1[ijk]*_gamma_D1D1[ijk] + 2.0*_A_UiUj_U1U1[ijk]*
_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D1[ijk]*_gamma_D1D2[ijk] +
_A_UiUj_U1U1[ijk]*pow(_HS_U1[ijk], 2)*pow(_gamma_D1D1[ijk], 2) + 2.0*
_A_UiUj_U1U1[ijk]*_HS_U1[ijk]*_HS_U2[ijk]*_gamma_D1D1[ijk]*
_gamma_D1D2[ijk] + _A_UiUj_U1U1[ijk]*pow(_HS_U2[ijk], 2)*
pow(_gamma_D1D2[ijk], 2) + 2.0*_A_UiUj_U1U2[ijk]*pow(_HS_U0[ijk], 2)*
_gamma_D0D1[ijk]*_gamma_D0D2[ijk] + 2.0*_A_UiUj_U1U2[ijk]*_HS_U0[ijk]*
_HS_U1[ijk]*_gamma_D0D1[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UiUj_U1U2[ijk]*
_HS_U0[ijk]*_HS_U1[ijk]*_gamma_D0D2[ijk]*_gamma_D1D1[ijk] + 2.0*
_A_UiUj_U1U2[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D1[ijk]*
_gamma_D2D2[ijk] + 2.0*_A_UiUj_U1U2[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*
_gamma_D0D2[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UiUj_U1U2[ijk]*
pow(_HS_U1[ijk], 2)*_gamma_D1D1[ijk]*_gamma_D1D2[ijk] + 2.0*
_A_UiUj_U1U2[ijk]*_HS_U1[ijk]*_HS_U2[ijk]*_gamma_D1D1[ijk]*
_gamma_D2D2[ijk] + 2.0*_A_UiUj_U1U2[ijk]*_HS_U1[ijk]*_HS_U2[ijk]*
pow(_gamma_D1D2[ijk], 2) + 2.0*_A_UiUj_U1U2[ijk]*pow(_HS_U2[ijk], 2)*
_gamma_D1D2[ijk]*_gamma_D2D2[ijk] + _A_UiUj_U2U2[ijk]*
pow(_HS_U0[ijk], 2)*pow(_gamma_D0D2[ijk], 2) + 2.0*_A_UiUj_U2U2[ijk]*
_HS_U0[ijk]*_HS_U1[ijk]*_gamma_D0D2[ijk]*_gamma_D1D2[ijk] + 2.0*
_A_UiUj_U2U2[ijk]*_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D2[ijk]*
_gamma_D2D2[ijk] + _A_UiUj_U2U2[ijk]*pow(_HS_U1[ijk], 2)*
pow(_gamma_D1D2[ijk], 2) + 2.0*_A_UiUj_U2U2[ijk]*_HS_U1[ijk]*
_HS_U2[ijk]*_gamma_D1D2[ijk]*_gamma_D2D2[ijk] + _A_UiUj_U2U2[ijk]*
pow(_HS_U2[ijk], 2)*pow(_gamma_D2D2[ijk], 2))/pow(psi[ijk], 5);

  double innerB_Epart = 
(-1.0/3.0*kd[ijk==lmn]*pow(psi[ijk], 2)*(3*E_t1 + K[ijk]*psi[ijk]) -
kd[ijk==lmn]*(_HS_U0[ijk]*dpsi_D0[ijk] + _HS_U1[ijk]*dpsi_D1[ijk] +
_HS_U2[ijk]*dpsi_D2[ijk]) + psi[ijk]*(Jpsi_D0(j_Jpsi_D0,ijk,lmn)*
_HS_U0[ijk] + Jpsi_D1(j_Jpsi_D1,ijk,lmn)*_HS_U1[ijk] +
Jpsi_D2(j_Jpsi_D2,ijk,lmn)*_HS_U2[ijk]))/pow(psi[ijk], 2);

  E_Trans[j][i] = innerB_Epart;

  DDM_SCHUR_JACOBIAN_BC_Epart_CLOSE

  }/* end of else if (patch->innerB) */

  return 0;
}