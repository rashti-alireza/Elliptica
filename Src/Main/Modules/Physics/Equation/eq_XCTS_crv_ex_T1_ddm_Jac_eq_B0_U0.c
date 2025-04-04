/*
  These C codes generated by Cpi version 3.0
  Copyright (C) 2019-2022 Alireza Rashti.
*/


#include "eq_header.h"
#include "maths_equation_solvings_lib.h"

void *eq_XCTS_curve_exc_T1_ddm_jacobian_eq_B0_U0(void *vp1,void *vp2);
void *eq_XCTS_curve_exc_T1_ddm_jacobian_eq_B0_U0(void *vp1,void *vp2)
{
  DDM_SCHUR_JACOBIAN_EQ_DECLARE
  Uint ijk,lmn;/* for Jacobian entries J[ijk][lmn] */
  const double kd[2]   = {0.,1.};/* Kronecker delta */

  Header_Jacobian
  Init_Jacobian(JB0_D0)
  Init_Jacobian(JB0_D1)
  Init_Jacobian(JB0_D2)
  Init_Jacobian(JJB0_D0D0)
  Init_Jacobian(JJB0_D0D1)
  Init_Jacobian(JJB0_D0D2)
  Init_Jacobian(JJB0_D1D1)
  Init_Jacobian(JJB0_D1D2)
  Init_Jacobian(JJB0_D2D2)
  READ_v(psi)
  READ_v_UNUSED(dpsi_D0)
  READ_v_UNUSED(dpsi_D1)
  READ_v_UNUSED(dpsi_D2)
  READ_v(alphaPsi)
  READ_v_UNUSED(dalphaPsi_D0)
  READ_v_UNUSED(dalphaPsi_D1)
  READ_v_UNUSED(dalphaPsi_D2)
  READ_v_UNUSED(RicciConf_D0D0)
  READ_v_UNUSED(RicciConf_D0D1)
  READ_v_UNUSED(RicciConf_D0D2)
  READ_v_UNUSED(RicciConf_D1D1)
  READ_v_UNUSED(RicciConf_D1D2)
  READ_v_UNUSED(RicciConf_D2D2)
  READ_v_UNUSED(igConf_U0U0)
  READ_v_UNUSED(igConf_U0U1)
  READ_v_UNUSED(igConf_U0U2)
  READ_v_UNUSED(igConf_U1U1)
  READ_v_UNUSED(igConf_U1U2)
  READ_v_UNUSED(igConf_U2U2)
  READ_v_UNUSED(ChrisConf_U0D0D0)
  READ_v_UNUSED(ChrisConf_U0D0D1)
  READ_v_UNUSED(ChrisConf_U0D0D2)
  READ_v_UNUSED(ChrisConf_U0D1D1)
  READ_v_UNUSED(ChrisConf_U0D1D2)
  READ_v_UNUSED(ChrisConf_U0D2D2)
  READ_v_UNUSED(ChrisConf_U1D0D0)
  READ_v_UNUSED(ChrisConf_U1D0D1)
  READ_v_UNUSED(ChrisConf_U1D0D2)
  READ_v_UNUSED(ChrisConf_U1D1D1)
  READ_v_UNUSED(ChrisConf_U1D1D2)
  READ_v_UNUSED(ChrisConf_U1D2D2)
  READ_v_UNUSED(ChrisConf_U2D0D0)
  READ_v_UNUSED(ChrisConf_U2D0D1)
  READ_v_UNUSED(ChrisConf_U2D0D2)
  READ_v_UNUSED(ChrisConf_U2D1D1)
  READ_v_UNUSED(ChrisConf_U2D1D2)
  READ_v_UNUSED(ChrisConf_U2D2D2)
  READ_v_UNUSED(dChrisConf_U0D0D0D0)
  READ_v_UNUSED(dChrisConf_U0D0D0D1)
  READ_v_UNUSED(dChrisConf_U0D0D0D2)
  READ_v_UNUSED(dChrisConf_U0D0D1D0)
  READ_v_UNUSED(dChrisConf_U0D0D1D1)
  READ_v_UNUSED(dChrisConf_U0D0D1D2)
  READ_v_UNUSED(dChrisConf_U0D0D2D0)
  READ_v_UNUSED(dChrisConf_U0D0D2D1)
  READ_v_UNUSED(dChrisConf_U0D0D2D2)
  READ_v_UNUSED(dChrisConf_U0D1D1D0)
  READ_v_UNUSED(dChrisConf_U0D1D1D1)
  READ_v_UNUSED(dChrisConf_U0D1D1D2)
  READ_v_UNUSED(dChrisConf_U0D1D2D0)
  READ_v_UNUSED(dChrisConf_U0D1D2D1)
  READ_v_UNUSED(dChrisConf_U0D1D2D2)
  READ_v_UNUSED(dChrisConf_U0D2D2D0)
  READ_v_UNUSED(dChrisConf_U0D2D2D1)
  READ_v_UNUSED(dChrisConf_U0D2D2D2)
  READ_v_UNUSED(dChrisConf_U1D0D0D0)
  READ_v_UNUSED(dChrisConf_U1D0D0D1)
  READ_v_UNUSED(dChrisConf_U1D0D0D2)
  READ_v_UNUSED(dChrisConf_U1D0D1D0)
  READ_v_UNUSED(dChrisConf_U1D0D1D1)
  READ_v_UNUSED(dChrisConf_U1D0D1D2)
  READ_v_UNUSED(dChrisConf_U1D0D2D0)
  READ_v_UNUSED(dChrisConf_U1D0D2D1)
  READ_v_UNUSED(dChrisConf_U1D0D2D2)
  READ_v_UNUSED(dChrisConf_U1D1D1D0)
  READ_v_UNUSED(dChrisConf_U1D1D1D1)
  READ_v_UNUSED(dChrisConf_U1D1D1D2)
  READ_v_UNUSED(dChrisConf_U1D1D2D0)
  READ_v_UNUSED(dChrisConf_U1D1D2D1)
  READ_v_UNUSED(dChrisConf_U1D1D2D2)
  READ_v_UNUSED(dChrisConf_U1D2D2D0)
  READ_v_UNUSED(dChrisConf_U1D2D2D1)
  READ_v_UNUSED(dChrisConf_U1D2D2D2)
  READ_v_UNUSED(dChrisConf_U2D0D0D0)
  READ_v_UNUSED(dChrisConf_U2D0D0D1)
  READ_v_UNUSED(dChrisConf_U2D0D0D2)
  READ_v_UNUSED(dChrisConf_U2D0D1D0)
  READ_v_UNUSED(dChrisConf_U2D0D1D1)
  READ_v_UNUSED(dChrisConf_U2D0D1D2)
  READ_v_UNUSED(dChrisConf_U2D0D2D0)
  READ_v_UNUSED(dChrisConf_U2D0D2D1)
  READ_v_UNUSED(dChrisConf_U2D0D2D2)
  READ_v_UNUSED(dChrisConf_U2D1D1D0)
  READ_v_UNUSED(dChrisConf_U2D1D1D1)
  READ_v_UNUSED(dChrisConf_U2D1D1D2)
  READ_v_UNUSED(dChrisConf_U2D1D2D0)
  READ_v_UNUSED(dChrisConf_U2D1D2D1)
  READ_v_UNUSED(dChrisConf_U2D1D2D2)
  READ_v_UNUSED(dChrisConf_U2D2D2D0)
  READ_v_UNUSED(dChrisConf_U2D2D2D1)
  READ_v_UNUSED(dChrisConf_U2D2D2D2)
  READ_v(JConfC)
  DDM_SCHUR_JACOBIAN_EQ_Bpart_OPEN

  double JB0_D0 = d2f_dxdu_Jacobian(patch,0,ijk,lmn,JB0_D0);
  double JB0_D1 = d2f_dxdu_Jacobian(patch,1,ijk,lmn,JB0_D1);
  double JB0_D2 = d2f_dxdu_Jacobian(patch,2,ijk,lmn,JB0_D2);
  double JJB0_D0D0 = d3f_dx2du_Jacobian(patch,0,ijk,lmn,JJB0_D0D0);
  double JJB0_D0D1 = d3f_dx2du_Jacobian(patch,1,ijk,lmn,JJB0_D0D1);
  double JJB0_D0D2 = d3f_dx2du_Jacobian(patch,2,ijk,lmn,JJB0_D0D2);
  double JJB0_D1D1 = d3f_dx2du_Jacobian(patch,3,ijk,lmn,JJB0_D1D1);
  double JJB0_D1D2 = d3f_dx2du_Jacobian(patch,4,ijk,lmn,JJB0_D1D2);
  double JJB0_D2D2 = d3f_dx2du_Jacobian(patch,5,ijk,lmn,JJB0_D2D2);
  double dLnOf_alpha_B_U0 =
-7*dpsi_D0[ijk]/psi[ijk] + dalphaPsi_D0[ijk]/alphaPsi[ijk];

  double dLnOf_alpha_B_U1 =
-7*dpsi_D1[ijk]/psi[ijk] + dalphaPsi_D1[ijk]/alphaPsi[ijk];

  double dLnOf_alpha_B_U2 =
-7*dpsi_D2[ijk]/psi[ijk] + dalphaPsi_D2[ijk]/alphaPsi[ijk];

  double t1_B_U0 =
igConf_U0U0[ijk]*(JJB0_D0D0 + dChrisConf_U0D0D0D0[ijk]*kd[ijk==lmn]) +
igConf_U0U1[ijk]*(JJB0_D0D1 + dChrisConf_U0D0D0D1[ijk]*kd[ijk==lmn]) +
igConf_U0U1[ijk]*(JJB0_D0D1 + dChrisConf_U0D0D1D0[ijk]*kd[ijk==lmn]) +
igConf_U0U2[ijk]*(JJB0_D0D2 + dChrisConf_U0D0D0D2[ijk]*kd[ijk==lmn]) +
igConf_U0U2[ijk]*(JJB0_D0D2 + dChrisConf_U0D0D2D0[ijk]*kd[ijk==lmn]) +
igConf_U1U1[ijk]*(JJB0_D1D1 + dChrisConf_U0D0D1D1[ijk]*kd[ijk==lmn]) +
igConf_U1U2[ijk]*(JJB0_D1D2 + dChrisConf_U0D0D1D2[ijk]*kd[ijk==lmn]) +
igConf_U1U2[ijk]*(JJB0_D1D2 + dChrisConf_U0D0D2D1[ijk]*kd[ijk==lmn]) +
igConf_U2U2[ijk]*(JJB0_D2D2 + dChrisConf_U0D0D2D2[ijk]*
kd[ijk==lmn]);



  double t2_B_U0 =
2.0*ChrisConf_U0D0D0[ijk]*JB0_D0*igConf_U0U0[ijk] + 2.0*
ChrisConf_U0D0D1[ijk]*JB0_D1*igConf_U1U1[ijk] + 2.0*
ChrisConf_U0D0D2[ijk]*JB0_D2*igConf_U2U2[ijk] + 2.0*igConf_U0U1[ijk]*
(ChrisConf_U0D0D0[ijk]*JB0_D1 + ChrisConf_U0D0D1[ijk]*JB0_D0) + 2.0*
igConf_U0U2[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D2 + ChrisConf_U0D0D2[ijk]*
JB0_D0) + 2.0*igConf_U1U2[ijk]*(ChrisConf_U0D0D1[ijk]*JB0_D2 +
ChrisConf_U0D0D2[ijk]*JB0_D1);



  double t3_B_U0 =
kd[ijk==lmn]*(pow(ChrisConf_U0D0D0[ijk], 2)*igConf_U0U0[ijk] + 2.0*
ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D1[ijk]*igConf_U0U1[ijk] + 2.0*
ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D2[ijk]*igConf_U0U2[ijk] +
pow(ChrisConf_U0D0D1[ijk], 2)*igConf_U1U1[ijk] + 2.0*
ChrisConf_U0D0D1[ijk]*ChrisConf_U0D0D2[ijk]*igConf_U1U2[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D0[ijk]*igConf_U0U0[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk]*igConf_U0U1[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D2[ijk]*igConf_U0U2[ijk] +
pow(ChrisConf_U0D0D2[ijk], 2)*igConf_U2U2[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D0[ijk]*igConf_U0U0[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U1[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D2[ijk]*igConf_U0U2[ijk] + ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D0[ijk]*igConf_U0U1[ijk] + ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D1[ijk]*igConf_U1U1[ijk] + ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D2[ijk]*igConf_U1U2[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U1D0D0[ijk]*igConf_U0U2[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U1D0D1[ijk]*igConf_U1U2[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U1D0D2[ijk]*igConf_U2U2[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U2D0D0[ijk]*igConf_U0U1[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U1U1[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U2D0D2[ijk]*igConf_U1U2[ijk] + ChrisConf_U0D2D2[ijk]*
ChrisConf_U2D0D0[ijk]*igConf_U0U2[ijk] + ChrisConf_U0D2D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U1U2[ijk] + ChrisConf_U0D2D2[ijk]*
ChrisConf_U2D0D2[ijk]*igConf_U2U2[ijk]);



  double t4_B_U0 =
-igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D0 + ChrisConf_U1D0D0[ijk]*
JB0_D1 + ChrisConf_U2D0D0[ijk]*JB0_D2 + kd[ijk==lmn]*
(pow(ChrisConf_U0D0D0[ijk], 2) + ChrisConf_U0D0D1[ijk]*
ChrisConf_U1D0D0[ijk] + ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D0[ijk])) -
2.0*igConf_U0U1[ijk]*(ChrisConf_U0D0D1[ijk]*JB0_D0 +
ChrisConf_U1D0D1[ijk]*JB0_D1 + ChrisConf_U2D0D1[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D1[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D1[ijk])) - 2.0*igConf_U0U2[ijk]*(ChrisConf_U0D0D2[ijk]*
JB0_D0 + ChrisConf_U1D0D2[ijk]*JB0_D1 + ChrisConf_U2D0D2[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D2[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D2[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D2[ijk])) - igConf_U1U1[ijk]*(ChrisConf_U0D1D1[ijk]*
JB0_D0 + ChrisConf_U1D1D1[ijk]*JB0_D1 + ChrisConf_U2D1D1[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D1D1[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D1D1[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D1D1[ijk])) - 2.0*igConf_U1U2[ijk]*(ChrisConf_U0D1D2[ijk]*
JB0_D0 + ChrisConf_U1D1D2[ijk]*JB0_D1 + ChrisConf_U2D1D2[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D1D2[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D1D2[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D1D2[ijk])) - igConf_U2U2[ijk]*(ChrisConf_U0D2D2[ijk]*
JB0_D0 + ChrisConf_U1D2D2[ijk]*JB0_D1 + ChrisConf_U2D2D2[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D2D2[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D2D2[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D2D2[ijk]));



  double t5_B_U0 =
0.33333333333333331*igConf_U0U0[ijk]*(JJB0_D0D0 + kd[ijk==lmn]*
(dChrisConf_U0D0D0D0[ijk] + dChrisConf_U1D0D1D0[ijk] +
dChrisConf_U2D0D2D0[ijk])) + 0.33333333333333331*igConf_U0U1[ijk]*
(JJB0_D0D1 + kd[ijk==lmn]*(dChrisConf_U0D0D0D1[ijk] +
dChrisConf_U1D0D1D1[ijk] + dChrisConf_U2D0D2D1[ijk])) +
0.33333333333333331*igConf_U0U2[ijk]*(JJB0_D0D2 + kd[ijk==lmn]*
(dChrisConf_U0D0D0D2[ijk] + dChrisConf_U1D0D1D2[ijk] +
dChrisConf_U2D0D2D2[ijk]));



  double t6_B_U0 =
0.33333333333333331*igConf_U0U0[ijk]*(2.0*ChrisConf_U0D0D0[ijk]*JB0_D0 +
ChrisConf_U1D0D0[ijk]*JB0_D1 + ChrisConf_U1D0D1[ijk]*JB0_D0 +
ChrisConf_U2D0D0[ijk]*JB0_D2 + ChrisConf_U2D0D2[ijk]*JB0_D0) +
0.33333333333333331*igConf_U0U1[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D1 +
ChrisConf_U0D0D1[ijk]*JB0_D0 + 2.0*ChrisConf_U1D0D1[ijk]*JB0_D1 +
ChrisConf_U2D0D1[ijk]*JB0_D2 + ChrisConf_U2D0D2[ijk]*JB0_D1) +
0.33333333333333331*igConf_U0U2[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D2 +
ChrisConf_U0D0D2[ijk]*JB0_D0 + ChrisConf_U1D0D1[ijk]*JB0_D2 +
ChrisConf_U1D0D2[ijk]*JB0_D1 + 2.0*ChrisConf_U2D0D2[ijk]*
JB0_D2);



  double t7_B_U0 =
kd[ijk==lmn]*(0.33333333333333331*pow(ChrisConf_U0D0D0[ijk], 2)*
igConf_U0U0[ijk] + 0.33333333333333331*ChrisConf_U0D0D0[ijk]*
ChrisConf_U0D0D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D2[ijk]*igConf_U0U2[ijk] +
0.66666666666666663*ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D0[ijk]*
igConf_U0U0[ijk] + 0.33333333333333331*ChrisConf_U0D0D1[ijk]*
ChrisConf_U1D0D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D2[ijk]*igConf_U0U2[ijk] +
0.66666666666666663*ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D0[ijk]*
igConf_U0U0[ijk] + 0.33333333333333331*ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D2[ijk]*igConf_U0U2[ijk] +
0.33333333333333331*ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D0[ijk]*
igConf_U0U1[ijk] + 0.33333333333333331*ChrisConf_U0D1D2[ijk]*
ChrisConf_U1D0D0[ijk]*igConf_U0U2[ijk] + 0.33333333333333331*
ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D0[ijk]*igConf_U0U1[ijk] +
0.33333333333333331*ChrisConf_U0D2D2[ijk]*ChrisConf_U2D0D0[ijk]*
igConf_U0U2[ijk] + 0.33333333333333331*pow(ChrisConf_U1D0D1[ijk], 2)*
igConf_U0U0[ijk] + 0.33333333333333331*ChrisConf_U1D0D1[ijk]*
ChrisConf_U1D1D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D2[ijk]*igConf_U0U2[ijk] +
0.66666666666666663*ChrisConf_U1D0D2[ijk]*ChrisConf_U2D0D1[ijk]*
igConf_U0U0[ijk] + 0.33333333333333331*ChrisConf_U1D0D2[ijk]*
ChrisConf_U2D1D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D2[ijk]*igConf_U0U2[ijk] +
0.33333333333333331*ChrisConf_U1D1D2[ijk]*ChrisConf_U2D0D1[ijk]*
igConf_U0U1[ijk] + 0.33333333333333331*ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U2[ijk] + 0.33333333333333331*
pow(ChrisConf_U2D0D2[ijk], 2)*igConf_U0U0[ijk] + 0.33333333333333331*
ChrisConf_U2D0D2[ijk]*ChrisConf_U2D1D2[ijk]*igConf_U0U1[ijk] +
0.33333333333333331*ChrisConf_U2D0D2[ijk]*ChrisConf_U2D2D2[ijk]*
igConf_U0U2[ijk]);



  double t8_B_U0 =
-0.33333333333333331*igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D0 +
ChrisConf_U1D0D0[ijk]*JB0_D1 + ChrisConf_U2D0D0[ijk]*JB0_D2 +
kd[ijk==lmn]*(pow(ChrisConf_U0D0D0[ijk], 2) + 2.0*ChrisConf_U0D0D1[ijk]*
ChrisConf_U1D0D0[ijk] + 2.0*ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D0[ijk] +
pow(ChrisConf_U1D0D1[ijk], 2) + 2.0*ChrisConf_U1D0D2[ijk]*
ChrisConf_U2D0D1[ijk] + pow(ChrisConf_U2D0D2[ijk], 2))) -
0.33333333333333331*igConf_U0U1[ijk]*(ChrisConf_U0D0D1[ijk]*JB0_D0 +
ChrisConf_U1D0D1[ijk]*JB0_D1 + ChrisConf_U2D0D1[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D1[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D1[ijk] + ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D0[ijk] +
ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D0[ijk] + ChrisConf_U1D0D1[ijk]*
ChrisConf_U1D1D1[ijk] + ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D1[ijk] +
ChrisConf_U1D1D2[ijk]*ChrisConf_U2D0D1[ijk] + ChrisConf_U2D0D2[ijk]*
ChrisConf_U2D1D2[ijk])) - 0.33333333333333331*igConf_U0U2[ijk]*
(ChrisConf_U0D0D2[ijk]*JB0_D0 + ChrisConf_U1D0D2[ijk]*JB0_D1 +
ChrisConf_U2D0D2[ijk]*JB0_D2 + kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*
ChrisConf_U0D0D2[ijk] + ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D2[ijk] +
ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D2[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U1D0D0[ijk] + ChrisConf_U0D2D2[ijk]*ChrisConf_U2D0D0[ijk] +
ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D2[ijk] + ChrisConf_U1D0D2[ijk]*
ChrisConf_U2D1D2[ijk] + ChrisConf_U1D2D2[ijk]*ChrisConf_U2D0D1[ijk] +
ChrisConf_U2D0D2[ijk]*ChrisConf_U2D2D2[ijk]));



  double t9_B_U0 =
kd[ijk==lmn]*(RicciConf_D0D0[ijk]*igConf_U0U0[ijk] +
RicciConf_D0D1[ijk]*igConf_U0U1[ijk] + RicciConf_D0D2[ijk]*
igConf_U0U2[ijk]);



  double contB =
JB0_D0 + kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk] + ChrisConf_U1D0D1[ijk] +
ChrisConf_U2D0D2[ijk]);

  double t1p10_B_U0 =
0.66666666666666663*contB*(dLnOf_alpha_B_U0*igConf_U0U0[ijk] +
dLnOf_alpha_B_U1*igConf_U0U1[ijk] + dLnOf_alpha_B_U2*
igConf_U0U2[ijk]);



  double t2p10_B_U0 =
-dLnOf_alpha_B_U0*igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*kd[ijk==lmn] +
JB0_D0) - dLnOf_alpha_B_U0*igConf_U0U1[ijk]*(ChrisConf_U0D0D1[ijk]*
kd[ijk==lmn] + JB0_D1) - dLnOf_alpha_B_U0*igConf_U0U2[ijk]*
(ChrisConf_U0D0D2[ijk]*kd[ijk==lmn] + JB0_D2) - dLnOf_alpha_B_U1*
igConf_U0U1[ijk]*(ChrisConf_U0D0D0[ijk]*kd[ijk==lmn] + JB0_D0) -
dLnOf_alpha_B_U1*igConf_U1U1[ijk]*(ChrisConf_U0D0D1[ijk]*kd[ijk==lmn] +
JB0_D1) - dLnOf_alpha_B_U1*igConf_U1U2[ijk]*(ChrisConf_U0D0D2[ijk]*
kd[ijk==lmn] + JB0_D2) - dLnOf_alpha_B_U2*igConf_U0U2[ijk]*
(ChrisConf_U0D0D0[ijk]*kd[ijk==lmn] + JB0_D0) - dLnOf_alpha_B_U2*
igConf_U1U2[ijk]*(ChrisConf_U0D0D1[ijk]*kd[ijk==lmn] + JB0_D1) -
dLnOf_alpha_B_U2*igConf_U2U2[ijk]*(ChrisConf_U0D0D2[ijk]*kd[ijk==lmn] +
JB0_D2);



  double t3p10_B_U0 =
-ChrisConf_U1D0D0[ijk]*dLnOf_alpha_B_U1*igConf_U0U0[ijk]*kd[ijk==lmn] -
ChrisConf_U1D0D1[ijk]*dLnOf_alpha_B_U1*igConf_U0U1[ijk]*kd[ijk==lmn] -
ChrisConf_U1D0D2[ijk]*dLnOf_alpha_B_U1*igConf_U0U2[ijk]*kd[ijk==lmn] -
ChrisConf_U2D0D0[ijk]*dLnOf_alpha_B_U2*igConf_U0U0[ijk]*kd[ijk==lmn] -
ChrisConf_U2D0D1[ijk]*dLnOf_alpha_B_U2*igConf_U0U1[ijk]*kd[ijk==lmn] -
ChrisConf_U2D0D2[ijk]*dLnOf_alpha_B_U2*igConf_U0U2[ijk]*kd[ijk==lmn] -
dLnOf_alpha_B_U0*igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*kd[ijk==lmn] +
JB0_D0) - dLnOf_alpha_B_U0*igConf_U0U1[ijk]*(ChrisConf_U0D0D1[ijk]*
kd[ijk==lmn] + JB0_D1) - dLnOf_alpha_B_U0*igConf_U0U2[ijk]*
(ChrisConf_U0D0D2[ijk]*kd[ijk==lmn] + JB0_D2);



  double t10_B_U0 =
t1p10_B_U0 + t2p10_B_U0 + t3p10_B_U0;



  double Bpart_U0 =
t10_B_U0 + t1_B_U0 + t2_B_U0 + t3_B_U0 + t4_B_U0 + t5_B_U0 + t6_B_U0 +
t7_B_U0 + t8_B_U0 + t9_B_U0;



  Bpart_U0 *= JConfC[ijk];
  schur_B[schur_ijk][schur_lmn] = Bpart_U0;

  DDM_SCHUR_JACOBIAN_EQ_Bpart_CLOSE

  DDM_SCHUR_JACOBIAN_EQ_Epart_OPEN

  double JB0_D0 = d2f_dxdu_Jacobian(patch,0,ijk,lmn,JB0_D0);
  double JB0_D1 = d2f_dxdu_Jacobian(patch,1,ijk,lmn,JB0_D1);
  double JB0_D2 = d2f_dxdu_Jacobian(patch,2,ijk,lmn,JB0_D2);
  double JJB0_D0D0 = d3f_dx2du_Jacobian(patch,0,ijk,lmn,JJB0_D0D0);
  double JJB0_D0D1 = d3f_dx2du_Jacobian(patch,1,ijk,lmn,JJB0_D0D1);
  double JJB0_D0D2 = d3f_dx2du_Jacobian(patch,2,ijk,lmn,JJB0_D0D2);
  double JJB0_D1D1 = d3f_dx2du_Jacobian(patch,3,ijk,lmn,JJB0_D1D1);
  double JJB0_D1D2 = d3f_dx2du_Jacobian(patch,4,ijk,lmn,JJB0_D1D2);
  double JJB0_D2D2 = d3f_dx2du_Jacobian(patch,5,ijk,lmn,JJB0_D2D2);
  double dLnOf_alpha_E_U0 =
-7*dpsi_D0[ijk]/psi[ijk] + dalphaPsi_D0[ijk]/alphaPsi[ijk];

  double dLnOf_alpha_E_U1 =
-7*dpsi_D1[ijk]/psi[ijk] + dalphaPsi_D1[ijk]/alphaPsi[ijk];

  double dLnOf_alpha_E_U2 =
-7*dpsi_D2[ijk]/psi[ijk] + dalphaPsi_D2[ijk]/alphaPsi[ijk];

  double t1_E_U0 =
igConf_U0U0[ijk]*(JJB0_D0D0 + dChrisConf_U0D0D0D0[ijk]*kd[ijk==lmn]) +
igConf_U0U1[ijk]*(JJB0_D0D1 + dChrisConf_U0D0D0D1[ijk]*kd[ijk==lmn]) +
igConf_U0U1[ijk]*(JJB0_D0D1 + dChrisConf_U0D0D1D0[ijk]*kd[ijk==lmn]) +
igConf_U0U2[ijk]*(JJB0_D0D2 + dChrisConf_U0D0D0D2[ijk]*kd[ijk==lmn]) +
igConf_U0U2[ijk]*(JJB0_D0D2 + dChrisConf_U0D0D2D0[ijk]*kd[ijk==lmn]) +
igConf_U1U1[ijk]*(JJB0_D1D1 + dChrisConf_U0D0D1D1[ijk]*kd[ijk==lmn]) +
igConf_U1U2[ijk]*(JJB0_D1D2 + dChrisConf_U0D0D1D2[ijk]*kd[ijk==lmn]) +
igConf_U1U2[ijk]*(JJB0_D1D2 + dChrisConf_U0D0D2D1[ijk]*kd[ijk==lmn]) +
igConf_U2U2[ijk]*(JJB0_D2D2 + dChrisConf_U0D0D2D2[ijk]*
kd[ijk==lmn]);



  double t2_E_U0 =
2.0*ChrisConf_U0D0D0[ijk]*JB0_D0*igConf_U0U0[ijk] + 2.0*
ChrisConf_U0D0D1[ijk]*JB0_D1*igConf_U1U1[ijk] + 2.0*
ChrisConf_U0D0D2[ijk]*JB0_D2*igConf_U2U2[ijk] + 2.0*igConf_U0U1[ijk]*
(ChrisConf_U0D0D0[ijk]*JB0_D1 + ChrisConf_U0D0D1[ijk]*JB0_D0) + 2.0*
igConf_U0U2[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D2 + ChrisConf_U0D0D2[ijk]*
JB0_D0) + 2.0*igConf_U1U2[ijk]*(ChrisConf_U0D0D1[ijk]*JB0_D2 +
ChrisConf_U0D0D2[ijk]*JB0_D1);



  double t3_E_U0 =
kd[ijk==lmn]*(pow(ChrisConf_U0D0D0[ijk], 2)*igConf_U0U0[ijk] + 2.0*
ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D1[ijk]*igConf_U0U1[ijk] + 2.0*
ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D2[ijk]*igConf_U0U2[ijk] +
pow(ChrisConf_U0D0D1[ijk], 2)*igConf_U1U1[ijk] + 2.0*
ChrisConf_U0D0D1[ijk]*ChrisConf_U0D0D2[ijk]*igConf_U1U2[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D0[ijk]*igConf_U0U0[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk]*igConf_U0U1[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D2[ijk]*igConf_U0U2[ijk] +
pow(ChrisConf_U0D0D2[ijk], 2)*igConf_U2U2[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D0[ijk]*igConf_U0U0[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U1[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D2[ijk]*igConf_U0U2[ijk] + ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D0[ijk]*igConf_U0U1[ijk] + ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D1[ijk]*igConf_U1U1[ijk] + ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D2[ijk]*igConf_U1U2[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U1D0D0[ijk]*igConf_U0U2[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U1D0D1[ijk]*igConf_U1U2[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U1D0D2[ijk]*igConf_U2U2[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U2D0D0[ijk]*igConf_U0U1[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U1U1[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U2D0D2[ijk]*igConf_U1U2[ijk] + ChrisConf_U0D2D2[ijk]*
ChrisConf_U2D0D0[ijk]*igConf_U0U2[ijk] + ChrisConf_U0D2D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U1U2[ijk] + ChrisConf_U0D2D2[ijk]*
ChrisConf_U2D0D2[ijk]*igConf_U2U2[ijk]);



  double t4_E_U0 =
-igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D0 + ChrisConf_U1D0D0[ijk]*
JB0_D1 + ChrisConf_U2D0D0[ijk]*JB0_D2 + kd[ijk==lmn]*
(pow(ChrisConf_U0D0D0[ijk], 2) + ChrisConf_U0D0D1[ijk]*
ChrisConf_U1D0D0[ijk] + ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D0[ijk])) -
2.0*igConf_U0U1[ijk]*(ChrisConf_U0D0D1[ijk]*JB0_D0 +
ChrisConf_U1D0D1[ijk]*JB0_D1 + ChrisConf_U2D0D1[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D1[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D1[ijk])) - 2.0*igConf_U0U2[ijk]*(ChrisConf_U0D0D2[ijk]*
JB0_D0 + ChrisConf_U1D0D2[ijk]*JB0_D1 + ChrisConf_U2D0D2[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D2[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D2[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D2[ijk])) - igConf_U1U1[ijk]*(ChrisConf_U0D1D1[ijk]*
JB0_D0 + ChrisConf_U1D1D1[ijk]*JB0_D1 + ChrisConf_U2D1D1[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D1D1[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D1D1[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D1D1[ijk])) - 2.0*igConf_U1U2[ijk]*(ChrisConf_U0D1D2[ijk]*
JB0_D0 + ChrisConf_U1D1D2[ijk]*JB0_D1 + ChrisConf_U2D1D2[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D1D2[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D1D2[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D1D2[ijk])) - igConf_U2U2[ijk]*(ChrisConf_U0D2D2[ijk]*
JB0_D0 + ChrisConf_U1D2D2[ijk]*JB0_D1 + ChrisConf_U2D2D2[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D2D2[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D2D2[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D2D2[ijk]));



  double t5_E_U0 =
0.33333333333333331*igConf_U0U0[ijk]*(JJB0_D0D0 + kd[ijk==lmn]*
(dChrisConf_U0D0D0D0[ijk] + dChrisConf_U1D0D1D0[ijk] +
dChrisConf_U2D0D2D0[ijk])) + 0.33333333333333331*igConf_U0U1[ijk]*
(JJB0_D0D1 + kd[ijk==lmn]*(dChrisConf_U0D0D0D1[ijk] +
dChrisConf_U1D0D1D1[ijk] + dChrisConf_U2D0D2D1[ijk])) +
0.33333333333333331*igConf_U0U2[ijk]*(JJB0_D0D2 + kd[ijk==lmn]*
(dChrisConf_U0D0D0D2[ijk] + dChrisConf_U1D0D1D2[ijk] +
dChrisConf_U2D0D2D2[ijk]));



  double t6_E_U0 =
0.33333333333333331*igConf_U0U0[ijk]*(2.0*ChrisConf_U0D0D0[ijk]*JB0_D0 +
ChrisConf_U1D0D0[ijk]*JB0_D1 + ChrisConf_U1D0D1[ijk]*JB0_D0 +
ChrisConf_U2D0D0[ijk]*JB0_D2 + ChrisConf_U2D0D2[ijk]*JB0_D0) +
0.33333333333333331*igConf_U0U1[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D1 +
ChrisConf_U0D0D1[ijk]*JB0_D0 + 2.0*ChrisConf_U1D0D1[ijk]*JB0_D1 +
ChrisConf_U2D0D1[ijk]*JB0_D2 + ChrisConf_U2D0D2[ijk]*JB0_D1) +
0.33333333333333331*igConf_U0U2[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D2 +
ChrisConf_U0D0D2[ijk]*JB0_D0 + ChrisConf_U1D0D1[ijk]*JB0_D2 +
ChrisConf_U1D0D2[ijk]*JB0_D1 + 2.0*ChrisConf_U2D0D2[ijk]*
JB0_D2);



  double t7_E_U0 =
kd[ijk==lmn]*(0.33333333333333331*pow(ChrisConf_U0D0D0[ijk], 2)*
igConf_U0U0[ijk] + 0.33333333333333331*ChrisConf_U0D0D0[ijk]*
ChrisConf_U0D0D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D2[ijk]*igConf_U0U2[ijk] +
0.66666666666666663*ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D0[ijk]*
igConf_U0U0[ijk] + 0.33333333333333331*ChrisConf_U0D0D1[ijk]*
ChrisConf_U1D0D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D2[ijk]*igConf_U0U2[ijk] +
0.66666666666666663*ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D0[ijk]*
igConf_U0U0[ijk] + 0.33333333333333331*ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D2[ijk]*igConf_U0U2[ijk] +
0.33333333333333331*ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D0[ijk]*
igConf_U0U1[ijk] + 0.33333333333333331*ChrisConf_U0D1D2[ijk]*
ChrisConf_U1D0D0[ijk]*igConf_U0U2[ijk] + 0.33333333333333331*
ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D0[ijk]*igConf_U0U1[ijk] +
0.33333333333333331*ChrisConf_U0D2D2[ijk]*ChrisConf_U2D0D0[ijk]*
igConf_U0U2[ijk] + 0.33333333333333331*pow(ChrisConf_U1D0D1[ijk], 2)*
igConf_U0U0[ijk] + 0.33333333333333331*ChrisConf_U1D0D1[ijk]*
ChrisConf_U1D1D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D2[ijk]*igConf_U0U2[ijk] +
0.66666666666666663*ChrisConf_U1D0D2[ijk]*ChrisConf_U2D0D1[ijk]*
igConf_U0U0[ijk] + 0.33333333333333331*ChrisConf_U1D0D2[ijk]*
ChrisConf_U2D1D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D2[ijk]*igConf_U0U2[ijk] +
0.33333333333333331*ChrisConf_U1D1D2[ijk]*ChrisConf_U2D0D1[ijk]*
igConf_U0U1[ijk] + 0.33333333333333331*ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U2[ijk] + 0.33333333333333331*
pow(ChrisConf_U2D0D2[ijk], 2)*igConf_U0U0[ijk] + 0.33333333333333331*
ChrisConf_U2D0D2[ijk]*ChrisConf_U2D1D2[ijk]*igConf_U0U1[ijk] +
0.33333333333333331*ChrisConf_U2D0D2[ijk]*ChrisConf_U2D2D2[ijk]*
igConf_U0U2[ijk]);



  double t8_E_U0 =
-0.33333333333333331*igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D0 +
ChrisConf_U1D0D0[ijk]*JB0_D1 + ChrisConf_U2D0D0[ijk]*JB0_D2 +
kd[ijk==lmn]*(pow(ChrisConf_U0D0D0[ijk], 2) + 2.0*ChrisConf_U0D0D1[ijk]*
ChrisConf_U1D0D0[ijk] + 2.0*ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D0[ijk] +
pow(ChrisConf_U1D0D1[ijk], 2) + 2.0*ChrisConf_U1D0D2[ijk]*
ChrisConf_U2D0D1[ijk] + pow(ChrisConf_U2D0D2[ijk], 2))) -
0.33333333333333331*igConf_U0U1[ijk]*(ChrisConf_U0D0D1[ijk]*JB0_D0 +
ChrisConf_U1D0D1[ijk]*JB0_D1 + ChrisConf_U2D0D1[ijk]*JB0_D2 +
kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D1[ijk] +
ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk] + ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D1[ijk] + ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D0[ijk] +
ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D0[ijk] + ChrisConf_U1D0D1[ijk]*
ChrisConf_U1D1D1[ijk] + ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D1[ijk] +
ChrisConf_U1D1D2[ijk]*ChrisConf_U2D0D1[ijk] + ChrisConf_U2D0D2[ijk]*
ChrisConf_U2D1D2[ijk])) - 0.33333333333333331*igConf_U0U2[ijk]*
(ChrisConf_U0D0D2[ijk]*JB0_D0 + ChrisConf_U1D0D2[ijk]*JB0_D1 +
ChrisConf_U2D0D2[ijk]*JB0_D2 + kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*
ChrisConf_U0D0D2[ijk] + ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D2[ijk] +
ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D2[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U1D0D0[ijk] + ChrisConf_U0D2D2[ijk]*ChrisConf_U2D0D0[ijk] +
ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D2[ijk] + ChrisConf_U1D0D2[ijk]*
ChrisConf_U2D1D2[ijk] + ChrisConf_U1D2D2[ijk]*ChrisConf_U2D0D1[ijk] +
ChrisConf_U2D0D2[ijk]*ChrisConf_U2D2D2[ijk]));



  double t9_E_U0 =
kd[ijk==lmn]*(RicciConf_D0D0[ijk]*igConf_U0U0[ijk] +
RicciConf_D0D1[ijk]*igConf_U0U1[ijk] + RicciConf_D0D2[ijk]*
igConf_U0U2[ijk]);



  double contE =
JB0_D0 + kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk] + ChrisConf_U1D0D1[ijk] +
ChrisConf_U2D0D2[ijk]);

  double t1p10_E_U0 =
0.66666666666666663*contE*(dLnOf_alpha_E_U0*igConf_U0U0[ijk] +
dLnOf_alpha_E_U1*igConf_U0U1[ijk] + dLnOf_alpha_E_U2*
igConf_U0U2[ijk]);



  double t2p10_E_U0 =
-dLnOf_alpha_E_U0*igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*kd[ijk==lmn] +
JB0_D0) - dLnOf_alpha_E_U0*igConf_U0U1[ijk]*(ChrisConf_U0D0D1[ijk]*
kd[ijk==lmn] + JB0_D1) - dLnOf_alpha_E_U0*igConf_U0U2[ijk]*
(ChrisConf_U0D0D2[ijk]*kd[ijk==lmn] + JB0_D2) - dLnOf_alpha_E_U1*
igConf_U0U1[ijk]*(ChrisConf_U0D0D0[ijk]*kd[ijk==lmn] + JB0_D0) -
dLnOf_alpha_E_U1*igConf_U1U1[ijk]*(ChrisConf_U0D0D1[ijk]*kd[ijk==lmn] +
JB0_D1) - dLnOf_alpha_E_U1*igConf_U1U2[ijk]*(ChrisConf_U0D0D2[ijk]*
kd[ijk==lmn] + JB0_D2) - dLnOf_alpha_E_U2*igConf_U0U2[ijk]*
(ChrisConf_U0D0D0[ijk]*kd[ijk==lmn] + JB0_D0) - dLnOf_alpha_E_U2*
igConf_U1U2[ijk]*(ChrisConf_U0D0D1[ijk]*kd[ijk==lmn] + JB0_D1) -
dLnOf_alpha_E_U2*igConf_U2U2[ijk]*(ChrisConf_U0D0D2[ijk]*kd[ijk==lmn] +
JB0_D2);



  double t3p10_E_U0 =
-ChrisConf_U1D0D0[ijk]*dLnOf_alpha_E_U1*igConf_U0U0[ijk]*kd[ijk==lmn] -
ChrisConf_U1D0D1[ijk]*dLnOf_alpha_E_U1*igConf_U0U1[ijk]*kd[ijk==lmn] -
ChrisConf_U1D0D2[ijk]*dLnOf_alpha_E_U1*igConf_U0U2[ijk]*kd[ijk==lmn] -
ChrisConf_U2D0D0[ijk]*dLnOf_alpha_E_U2*igConf_U0U0[ijk]*kd[ijk==lmn] -
ChrisConf_U2D0D1[ijk]*dLnOf_alpha_E_U2*igConf_U0U1[ijk]*kd[ijk==lmn] -
ChrisConf_U2D0D2[ijk]*dLnOf_alpha_E_U2*igConf_U0U2[ijk]*kd[ijk==lmn] -
dLnOf_alpha_E_U0*igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*kd[ijk==lmn] +
JB0_D0) - dLnOf_alpha_E_U0*igConf_U0U1[ijk]*(ChrisConf_U0D0D1[ijk]*
kd[ijk==lmn] + JB0_D1) - dLnOf_alpha_E_U0*igConf_U0U2[ijk]*
(ChrisConf_U0D0D2[ijk]*kd[ijk==lmn] + JB0_D2);



  double t10_E_U0 =
t1p10_E_U0 + t2p10_E_U0 + t3p10_E_U0;



  double Epart_U0 =
t10_E_U0 + t1_E_U0 + t2_E_U0 + t3_E_U0 + t4_E_U0 + t5_E_U0 + t6_E_U0 +
t7_E_U0 + t8_E_U0 + t9_E_U0;



  Epart_U0 *= JConfC[ijk];
  schur_Et[schur_lmn][schur_ijk] = Epart_U0;

  DDM_SCHUR_JACOBIAN_EQ_Epart_CLOSE


  Free_Jacobian(JB0_D0)
  Free_Jacobian(JB0_D1)
  Free_Jacobian(JB0_D2)
  Free_Jacobian(JJB0_D0D0)
  Free_Jacobian(JJB0_D0D1)
  Free_Jacobian(JJB0_D0D2)
  Free_Jacobian(JJB0_D1D1)
  Free_Jacobian(JJB0_D1D2)
  Free_Jacobian(JJB0_D2D2)
  Footer_Jacobian

  return 0;
}
