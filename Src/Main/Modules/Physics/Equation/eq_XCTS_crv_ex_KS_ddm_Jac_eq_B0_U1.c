/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2021 Alireza Rashti.
*/


#include "eq_header.h"
#include "maths_equation_solvings_lib.h"

void *eq_XCTS_curve_exc_KS_ddm_jacobian_eq_B0_U1(void *vp1,void *vp2);
void *eq_XCTS_curve_exc_KS_ddm_jacobian_eq_B0_U1(void *vp1,void *vp2)
{
  DDM_SCHUR_JACOBIAN_EQ_DECLARE
  Uint ijk,lmn;/* for Jacobian entries J[ijk][lmn] */
  const double kd[2]   = {0.,1.};/* Kronecker delta */

  /* declaring: */
  JACOBIAN_DERIVATIVE(JB0_D1)
  JACOBIAN_DERIVATIVE(JB0_D0)
  JACOBIAN_DERIVATIVE(JB0_D2)
  JACOBIAN_DERIVATIVE(JJB0_D1D1)
  JACOBIAN_DERIVATIVE(JJB0_D1D2)
  JACOBIAN_DERIVATIVE(JJB0_D2D2)
  JACOBIAN_DERIVATIVE(JJB0_D0D0)
  JACOBIAN_DERIVATIVE(JJB0_D0D1)
  JACOBIAN_DERIVATIVE(JJB0_D0D2)
  READ_v(psi)
  READ_v_UNUSED(dpsi_D0)
  READ_v_UNUSED(dpsi_D1)
  READ_v_UNUSED(dpsi_D2)
  READ_v(alphaPsi)
  READ_v_UNUSED(dalphaPsi_D2)
  READ_v_UNUSED(dalphaPsi_D1)
  READ_v_UNUSED(dalphaPsi_D0)
  READ_v_UNUSED(RicciConf_D0D1)
  READ_v_UNUSED(RicciConf_D0D0)
  READ_v_UNUSED(RicciConf_D0D2)
  READ_v_UNUSED(RicciConf_D1D1)
  READ_v_UNUSED(RicciConf_D1D2)
  READ_v_UNUSED(RicciConf_D2D2)
  READ_v_UNUSED(igConf_U2U2)
  READ_v_UNUSED(igConf_U1U2)
  READ_v_UNUSED(igConf_U1U1)
  READ_v_UNUSED(igConf_U0U2)
  READ_v_UNUSED(igConf_U0U0)
  READ_v_UNUSED(igConf_U0U1)
  READ_v_UNUSED(ChrisConf_U1D0D0)
  READ_v_UNUSED(ChrisConf_U0D2D2)
  READ_v_UNUSED(ChrisConf_U1D0D2)
  READ_v_UNUSED(ChrisConf_U2D2D2)
  READ_v_UNUSED(ChrisConf_U2D1D1)
  READ_v_UNUSED(ChrisConf_U2D0D0)
  READ_v_UNUSED(ChrisConf_U1D1D1)
  READ_v_UNUSED(ChrisConf_U0D1D1)
  READ_v_UNUSED(ChrisConf_U0D1D2)
  READ_v_UNUSED(ChrisConf_U1D1D2)
  READ_v_UNUSED(ChrisConf_U0D0D1)
  READ_v_UNUSED(ChrisConf_U0D0D0)
  READ_v_UNUSED(ChrisConf_U2D0D1)
  READ_v_UNUSED(ChrisConf_U0D0D2)
  READ_v_UNUSED(ChrisConf_U2D1D2)
  READ_v_UNUSED(ChrisConf_U1D2D2)
  READ_v_UNUSED(ChrisConf_U2D0D2)
  READ_v_UNUSED(ChrisConf_U1D0D1)
  READ_v_UNUSED(dChrisConf_U1D1D1D1)
  READ_v_UNUSED(dChrisConf_U2D2D2D0)
  READ_v_UNUSED(dChrisConf_U1D1D1D0)
  READ_v_UNUSED(dChrisConf_U2D1D2D2)
  READ_v_UNUSED(dChrisConf_U2D1D2D1)
  READ_v_UNUSED(dChrisConf_U2D1D2D0)
  READ_v_UNUSED(dChrisConf_U2D2D2D1)
  READ_v_UNUSED(dChrisConf_U0D2D2D2)
  READ_v_UNUSED(dChrisConf_U1D1D1D2)
  READ_v_UNUSED(dChrisConf_U0D2D2D0)
  READ_v_UNUSED(dChrisConf_U0D2D2D1)
  READ_v_UNUSED(dChrisConf_U2D1D1D0)
  READ_v_UNUSED(dChrisConf_U0D0D1D2)
  READ_v_UNUSED(dChrisConf_U0D0D1D1)
  READ_v_UNUSED(dChrisConf_U0D0D1D0)
  READ_v_UNUSED(dChrisConf_U0D1D2D1)
  READ_v_UNUSED(dChrisConf_U0D1D2D0)
  READ_v_UNUSED(dChrisConf_U0D1D2D2)
  READ_v_UNUSED(dChrisConf_U1D2D2D2)
  READ_v_UNUSED(dChrisConf_U1D2D2D1)
  READ_v_UNUSED(dChrisConf_U1D2D2D0)
  READ_v_UNUSED(dChrisConf_U2D0D0D0)
  READ_v_UNUSED(dChrisConf_U0D1D1D2)
  READ_v_UNUSED(dChrisConf_U0D1D1D0)
  READ_v_UNUSED(dChrisConf_U0D1D1D1)
  READ_v_UNUSED(dChrisConf_U1D0D2D1)
  READ_v_UNUSED(dChrisConf_U1D0D2D0)
  READ_v_UNUSED(dChrisConf_U2D2D2D2)
  READ_v_UNUSED(dChrisConf_U1D0D2D2)
  READ_v_UNUSED(dChrisConf_U2D0D1D2)
  READ_v_UNUSED(dChrisConf_U2D0D0D1)
  READ_v_UNUSED(dChrisConf_U1D0D0D2)
  READ_v_UNUSED(dChrisConf_U1D0D0D1)
  READ_v_UNUSED(dChrisConf_U1D0D0D0)
  READ_v_UNUSED(dChrisConf_U1D0D1D2)
  READ_v_UNUSED(dChrisConf_U1D0D1D0)
  READ_v_UNUSED(dChrisConf_U1D0D1D1)
  READ_v_UNUSED(dChrisConf_U2D0D2D2)
  READ_v_UNUSED(dChrisConf_U2D1D1D1)
  READ_v_UNUSED(dChrisConf_U0D0D0D2)
  READ_v_UNUSED(dChrisConf_U2D0D2D0)
  READ_v_UNUSED(dChrisConf_U0D0D0D0)
  READ_v_UNUSED(dChrisConf_U0D0D0D1)
  READ_v_UNUSED(dChrisConf_U2D0D0D2)
  READ_v_UNUSED(dChrisConf_U2D0D2D1)
  READ_v_UNUSED(dChrisConf_U1D1D2D0)
  READ_v_UNUSED(dChrisConf_U1D1D2D1)
  READ_v_UNUSED(dChrisConf_U1D1D2D2)
  READ_v_UNUSED(dChrisConf_U2D1D1D2)
  READ_v_UNUSED(dChrisConf_U0D0D2D0)
  READ_v_UNUSED(dChrisConf_U0D0D2D1)
  READ_v_UNUSED(dChrisConf_U0D0D2D2)
  READ_v_UNUSED(dChrisConf_U2D0D1D1)
  READ_v_UNUSED(dChrisConf_U2D0D1D0)


  DDM_SCHUR_JACOBIAN_EQ_Bpart_OPEN

  double dLnOf_alpha_B_U1 = 
-7*dpsi_D1[ijk]/psi[ijk] + dalphaPsi_D1[ijk]/alphaPsi[ijk];

  double dLnOf_alpha_B_U0 = 
-7*dpsi_D0[ijk]/psi[ijk] + dalphaPsi_D0[ijk]/alphaPsi[ijk];

  double dLnOf_alpha_B_U2 = 
-7*dpsi_D2[ijk]/psi[ijk] + dalphaPsi_D2[ijk]/alphaPsi[ijk];

  double t1_B_U1 = 
igConf_U0U0[ijk]*(JJB0_D0D0(j_JJB0_D0D0,ijk,lmn) + dChrisConf_U1D0D1D0[ijk]*
kd[ijk==lmn]) + igConf_U0U1[ijk]*(JJB0_D0D1(j_JJB0_D0D1,ijk,lmn) + 
dChrisConf_U1D0D1D1[ijk]*kd[ijk==lmn]) + igConf_U0U1[ijk]*
(JJB0_D0D1(j_JJB0_D0D1,ijk,lmn) + dChrisConf_U1D1D1D0[ijk]*
kd[ijk==lmn]) + igConf_U0U2[ijk]*(JJB0_D0D2(j_JJB0_D0D2,ijk,lmn) + 
dChrisConf_U1D0D1D2[ijk]*kd[ijk==lmn]) + igConf_U0U2[ijk]*
(JJB0_D0D2(j_JJB0_D0D2,ijk,lmn) + dChrisConf_U1D1D2D0[ijk]*
kd[ijk==lmn]) + igConf_U1U1[ijk]*(JJB0_D1D1(j_JJB0_D1D1,ijk,lmn) + 
dChrisConf_U1D1D1D1[ijk]*kd[ijk==lmn]) + igConf_U1U2[ijk]*
(JJB0_D1D2(j_JJB0_D1D2,ijk,lmn) + dChrisConf_U1D1D1D2[ijk]*
kd[ijk==lmn]) + igConf_U1U2[ijk]*(JJB0_D1D2(j_JJB0_D1D2,ijk,lmn) + 
dChrisConf_U1D1D2D1[ijk]*kd[ijk==lmn]) + igConf_U2U2[ijk]*
(JJB0_D2D2(j_JJB0_D2D2,ijk,lmn) + dChrisConf_U1D1D2D2[ijk]*
kd[ijk==lmn]);




  double t2_B_U1 = 
2.0*ChrisConf_U1D0D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn)*igConf_U0U0[ijk] + 
2.0*ChrisConf_U1D1D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn)*igConf_U1U1[ijk] + 
2.0*ChrisConf_U1D1D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn)*igConf_U2U2[ijk] + 
2.0*igConf_U0U1[ijk]*(ChrisConf_U1D0D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U1D1D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn)) + 2.0*igConf_U0U2[ijk]*
(ChrisConf_U1D0D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + ChrisConf_U1D1D2[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn)) + 2.0*igConf_U1U2[ijk]*(ChrisConf_U1D1D1[ijk]*
JB0_D2(j_JB0_D2,ijk,lmn) + ChrisConf_U1D1D2[ijk]*JB0_D1(j_JB0_D1,ijk,lmn));



  double t3_B_U1 = 
kd[ijk==lmn]*(ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D0[ijk]*
igConf_U0U0[ijk] + ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk]*
igConf_U0U1[ijk] + ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D2[ijk]*
igConf_U0U2[ijk] + ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D0[ijk]*
igConf_U0U1[ijk] + ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D1[ijk]*
igConf_U1U1[ijk] + ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D2[ijk]*
igConf_U1U2[ijk] + ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D0[ijk]*
igConf_U0U2[ijk] + ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D1[ijk]*
igConf_U1U2[ijk] + ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D2[ijk]*
igConf_U2U2[ijk] + pow(ChrisConf_U1D0D1[ijk], 2)*igConf_U0U0[ijk] + 2.0*
ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D1[ijk]*igConf_U0U1[ijk] + 2.0*
ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D2[ijk]*igConf_U0U2[ijk] + 
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D0D1[ijk]*igConf_U0U0[ijk] + 
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D1[ijk]*igConf_U0U1[ijk] + 
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D2[ijk]*igConf_U0U2[ijk] + 
pow(ChrisConf_U1D1D1[ijk], 2)*igConf_U1U1[ijk] + 2.0*
ChrisConf_U1D1D1[ijk]*ChrisConf_U1D1D2[ijk]*igConf_U1U2[ijk] + 
pow(ChrisConf_U1D1D2[ijk], 2)*igConf_U2U2[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U1[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D1D1[ijk]*igConf_U1U1[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D1D2[ijk]*igConf_U1U2[ijk] + ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U2[ijk] + ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D1D1[ijk]*igConf_U1U2[ijk] + ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D1D2[ijk]*igConf_U2U2[ijk]);




  double t4_B_U1 = 
-igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + 
ChrisConf_U1D0D0[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + ChrisConf_U2D0D0[ijk]*
JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*
ChrisConf_U1D0D1[ijk] + ChrisConf_U1D0D0[ijk]*ChrisConf_U1D1D1[ijk] + 
ChrisConf_U1D1D2[ijk]*ChrisConf_U2D0D0[ijk])) - 2.0*igConf_U0U1[ijk]*
(ChrisConf_U0D0D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D0D1[ijk]*
JB0_D1(j_JB0_D1,ijk,lmn) + ChrisConf_U2D0D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + 
kd[ijk==lmn]*(ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk] + 
ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D1[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D0D1[ijk])) - 2.0*igConf_U0U2[ijk]*(ChrisConf_U0D0D2[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D0D2[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D0D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*
(ChrisConf_U0D0D2[ijk]*ChrisConf_U1D0D1[ijk] + ChrisConf_U1D0D2[ijk]*
ChrisConf_U1D1D1[ijk] + ChrisConf_U1D1D2[ijk]*ChrisConf_U2D0D2[ijk])) - 
igConf_U1U1[ijk]*(ChrisConf_U0D1D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + 
ChrisConf_U1D1D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + ChrisConf_U2D1D1[ijk]*
JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*(ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D1[ijk] + pow(ChrisConf_U1D1D1[ijk], 2) + 
ChrisConf_U1D1D2[ijk]*ChrisConf_U2D1D1[ijk])) - 2.0*igConf_U1U2[ijk]*
(ChrisConf_U0D1D2[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D1D2[ijk]*
JB0_D1(j_JB0_D1,ijk,lmn) + ChrisConf_U2D1D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + 
kd[ijk==lmn]*(ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D1[ijk] + 
ChrisConf_U1D1D1[ijk]*ChrisConf_U1D1D2[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D1D2[ijk])) - igConf_U2U2[ijk]*(ChrisConf_U0D2D2[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D2D2[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D2D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*
(ChrisConf_U0D2D2[ijk]*ChrisConf_U1D0D1[ijk] + ChrisConf_U1D1D1[ijk]*
ChrisConf_U1D2D2[ijk] + ChrisConf_U1D1D2[ijk]*ChrisConf_U2D2D2[ijk]));

  double t5_B_U1 = 
0.33333333333333331*igConf_U0U1[ijk]*(JJB0_D0D1(j_JJB0_D0D1,ijk,lmn) + 
kd[ijk==lmn]*(dChrisConf_U0D0D1D0[ijk] + dChrisConf_U1D1D1D0[ijk] + 
dChrisConf_U2D1D2D0[ijk])) + 0.33333333333333331*igConf_U1U1[ijk]*
(JJB0_D1D1(j_JJB0_D1D1,ijk,lmn) + kd[ijk==lmn]*(dChrisConf_U0D0D1D1[ijk] + 
dChrisConf_U1D1D1D1[ijk] + dChrisConf_U2D1D2D1[ijk])) + 
0.33333333333333331*igConf_U1U2[ijk]*(JJB0_D1D2(j_JJB0_D1D2,ijk,lmn) + 
kd[ijk==lmn]*(dChrisConf_U0D0D1D2[ijk] + dChrisConf_U1D1D1D2[ijk] + 
dChrisConf_U2D1D2D2[ijk]));




  double t6_B_U1 = 
0.33333333333333331*igConf_U0U1[ijk]*(2.0*ChrisConf_U0D0D1[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D0D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U1D1D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U2D0D1[ijk]*
JB0_D2(j_JB0_D2,ijk,lmn) + ChrisConf_U2D1D2[ijk]*JB0_D0(j_JB0_D0,ijk,lmn)) + 
0.33333333333333331*igConf_U1U1[ijk]*(ChrisConf_U0D0D1[ijk]*
JB0_D1(j_JB0_D1,ijk,lmn) + ChrisConf_U0D1D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + 
2.0*ChrisConf_U1D1D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D1D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + ChrisConf_U2D1D2[ijk]*
JB0_D1(j_JB0_D1,ijk,lmn)) + 0.33333333333333331*igConf_U1U2[ijk]*
(ChrisConf_U0D0D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + ChrisConf_U0D1D2[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D1D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + 
ChrisConf_U1D1D2[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 2.0*
ChrisConf_U2D1D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn));



  double t7_B_U1 = 
kd[ijk==lmn]*(0.33333333333333331*ChrisConf_U0D0D0[ijk]*
ChrisConf_U0D0D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
pow(ChrisConf_U0D0D1[ijk], 2)*igConf_U1U1[ijk] + 0.33333333333333331*
ChrisConf_U0D0D1[ijk]*ChrisConf_U0D0D2[ijk]*igConf_U1U2[ijk] + 
0.33333333333333331*ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk]*
igConf_U0U1[ijk] + 0.33333333333333331*ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D0[ijk]*igConf_U0U1[ijk] + 
0.66666666666666663*ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D1[ijk]*
igConf_U1U1[ijk] + 0.33333333333333331*ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D2[ijk]*igConf_U1U2[ijk] + 0.33333333333333331*
ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D1[ijk]*igConf_U1U2[ijk] + 
0.33333333333333331*ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D0[ijk]*
igConf_U0U1[ijk] + 0.66666666666666663*ChrisConf_U0D1D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U1U1[ijk] + 0.33333333333333331*
ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D2[ijk]*igConf_U1U2[ijk] + 
0.33333333333333331*ChrisConf_U0D2D2[ijk]*ChrisConf_U2D0D1[ijk]*
igConf_U1U2[ijk] + 0.33333333333333331*ChrisConf_U1D0D1[ijk]*
ChrisConf_U1D1D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D1[ijk]*igConf_U0U1[ijk] + 
0.33333333333333331*pow(ChrisConf_U1D1D1[ijk], 2)*igConf_U1U1[ijk] + 
0.33333333333333331*ChrisConf_U1D1D1[ijk]*ChrisConf_U1D1D2[ijk]*
igConf_U1U2[ijk] + 0.33333333333333331*ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U1[ijk] + 0.66666666666666663*
ChrisConf_U1D1D2[ijk]*ChrisConf_U2D1D1[ijk]*igConf_U1U1[ijk] + 
0.33333333333333331*ChrisConf_U1D1D2[ijk]*ChrisConf_U2D1D2[ijk]*
igConf_U1U2[ijk] + 0.33333333333333331*ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D1D1[ijk]*igConf_U1U2[ijk] + 0.33333333333333331*
ChrisConf_U2D0D2[ijk]*ChrisConf_U2D1D2[ijk]*igConf_U0U1[ijk] + 
0.33333333333333331*pow(ChrisConf_U2D1D2[ijk], 2)*igConf_U1U1[ijk] + 
0.33333333333333331*ChrisConf_U2D1D2[ijk]*ChrisConf_U2D2D2[ijk]*
igConf_U1U2[ijk]);




  double t8_B_U1 = 
-0.33333333333333331*igConf_U0U1[ijk]*(ChrisConf_U0D0D1[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D0D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D0D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*
(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D1[ijk] + ChrisConf_U0D0D1[ijk]*
ChrisConf_U1D0D1[ijk] + ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D1[ijk] + 
ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D0[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U2D0D0[ijk] + ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D1[ijk] + 
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D1[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D0D1[ijk] + ChrisConf_U2D0D2[ijk]*ChrisConf_U2D1D2[ijk])) - 
0.33333333333333331*igConf_U1U1[ijk]*(ChrisConf_U0D1D1[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D1D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D1D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*
(pow(ChrisConf_U0D0D1[ijk], 2) + 2.0*ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D1[ijk] + 2.0*ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D1[ijk] + 
pow(ChrisConf_U1D1D1[ijk], 2) + 2.0*ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D1D1[ijk] + pow(ChrisConf_U2D1D2[ijk], 2))) - 
0.33333333333333331*igConf_U1U2[ijk]*(ChrisConf_U0D1D2[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D1D2[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D1D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*
(ChrisConf_U0D0D1[ijk]*ChrisConf_U0D0D2[ijk] + ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D2[ijk] + ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D1[ijk] + 
ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D2[ijk] + ChrisConf_U0D2D2[ijk]*
ChrisConf_U2D0D1[ijk] + ChrisConf_U1D1D1[ijk]*ChrisConf_U1D1D2[ijk] + 
ChrisConf_U1D1D2[ijk]*ChrisConf_U2D1D2[ijk] + ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D1D1[ijk] + ChrisConf_U2D1D2[ijk]*ChrisConf_U2D2D2[ijk]));

  double t9_B_U1 = 
kd[ijk==lmn]*(RicciConf_D0D1[ijk]*igConf_U0U1[ijk] + 
RicciConf_D1D1[ijk]*igConf_U1U1[ijk] + RicciConf_D1D2[ijk]*
igConf_U1U2[ijk]);



  double contB = 
JB0_D1(j_JB0_D1,ijk,lmn) + kd[ijk==lmn]*(ChrisConf_U0D0D1[ijk] +
ChrisConf_U1D1D1[ijk] + ChrisConf_U2D1D2[ijk]);


  double t1p10_B_U1 = 
0.66666666666666663*contB*(dLnOf_alpha_B_U0*igConf_U0U1[ijk] + 
dLnOf_alpha_B_U1*igConf_U1U1[ijk] + dLnOf_alpha_B_U2*
igConf_U1U2[ijk]);



  double t2p10_B_U1 = 
-dLnOf_alpha_B_U0*igConf_U0U0[ijk]*(ChrisConf_U1D0D1[ijk]*kd[ijk==lmn] + 
JB0_D0(j_JB0_D0,ijk,lmn)) - dLnOf_alpha_B_U0*igConf_U0U1[ijk]*
(ChrisConf_U1D1D1[ijk]*kd[ijk==lmn] + JB0_D1(j_JB0_D1,ijk,lmn)) - 
dLnOf_alpha_B_U0*igConf_U0U2[ijk]*(ChrisConf_U1D1D2[ijk]*kd[ijk==lmn] + 
JB0_D2(j_JB0_D2,ijk,lmn)) - dLnOf_alpha_B_U1*igConf_U0U1[ijk]*
(ChrisConf_U1D0D1[ijk]*kd[ijk==lmn] + JB0_D0(j_JB0_D0,ijk,lmn)) - 
dLnOf_alpha_B_U1*igConf_U1U1[ijk]*(ChrisConf_U1D1D1[ijk]*kd[ijk==lmn] + 
JB0_D1(j_JB0_D1,ijk,lmn)) - dLnOf_alpha_B_U1*igConf_U1U2[ijk]*
(ChrisConf_U1D1D2[ijk]*kd[ijk==lmn] + JB0_D2(j_JB0_D2,ijk,lmn)) - 
dLnOf_alpha_B_U2*igConf_U0U2[ijk]*(ChrisConf_U1D0D1[ijk]*kd[ijk==lmn] + 
JB0_D0(j_JB0_D0,ijk,lmn)) - dLnOf_alpha_B_U2*igConf_U1U2[ijk]*
(ChrisConf_U1D1D1[ijk]*kd[ijk==lmn] + JB0_D1(j_JB0_D1,ijk,lmn)) - 
dLnOf_alpha_B_U2*igConf_U2U2[ijk]*(ChrisConf_U1D1D2[ijk]*kd[ijk==lmn] + 
JB0_D2(j_JB0_D2,ijk,lmn));


  double t3p10_B_U1 = 
-ChrisConf_U0D0D1[ijk]*dLnOf_alpha_B_U0*igConf_U0U1[ijk]*kd[ijk==lmn] - 
ChrisConf_U0D1D1[ijk]*dLnOf_alpha_B_U0*igConf_U1U1[ijk]*kd[ijk==lmn] - 
ChrisConf_U0D1D2[ijk]*dLnOf_alpha_B_U0*igConf_U1U2[ijk]*kd[ijk==lmn] - 
ChrisConf_U2D0D1[ijk]*dLnOf_alpha_B_U2*igConf_U0U1[ijk]*kd[ijk==lmn] - 
ChrisConf_U2D1D1[ijk]*dLnOf_alpha_B_U2*igConf_U1U1[ijk]*kd[ijk==lmn] - 
ChrisConf_U2D1D2[ijk]*dLnOf_alpha_B_U2*igConf_U1U2[ijk]*kd[ijk==lmn] - 
dLnOf_alpha_B_U1*igConf_U0U1[ijk]*(ChrisConf_U1D0D1[ijk]*kd[ijk==lmn] + 
JB0_D0(j_JB0_D0,ijk,lmn)) - dLnOf_alpha_B_U1*igConf_U1U1[ijk]*
(ChrisConf_U1D1D1[ijk]*kd[ijk==lmn] + JB0_D1(j_JB0_D1,ijk,lmn)) - 
dLnOf_alpha_B_U1*igConf_U1U2[ijk]*(ChrisConf_U1D1D2[ijk]*kd[ijk==lmn] + 
JB0_D2(j_JB0_D2,ijk,lmn));




  double t10_B_U1 = 
t1p10_B_U1 + t2p10_B_U1 + t3p10_B_U1;



  double Bpart_U1 = 
t10_B_U1 + t1_B_U1 + t2_B_U1 + t3_B_U1 + t4_B_U1 + t5_B_U1 + t6_B_U1 + 
t7_B_U1 + t8_B_U1 + t9_B_U1;


  B[i][j] = Bpart_U1;

  DDM_SCHUR_JACOBIAN_EQ_Bpart_CLOSE

  DDM_SCHUR_JACOBIAN_EQ_Epart_OPEN

  double dLnOf_alpha_E_U0 = 
-7*dpsi_D0[ijk]/psi[ijk] + dalphaPsi_D0[ijk]/alphaPsi[ijk];

  double dLnOf_alpha_E_U1 = 
-7*dpsi_D1[ijk]/psi[ijk] + dalphaPsi_D1[ijk]/alphaPsi[ijk];

  double dLnOf_alpha_E_U2 = 
-7*dpsi_D2[ijk]/psi[ijk] + dalphaPsi_D2[ijk]/alphaPsi[ijk];


  double t1_E_U1 = 
igConf_U0U0[ijk]*(JJB0_D0D0(j_JJB0_D0D0,ijk,lmn) + dChrisConf_U1D0D1D0[ijk]*
kd[ijk==lmn]) + igConf_U0U1[ijk]*(JJB0_D0D1(j_JJB0_D0D1,ijk,lmn) + 
dChrisConf_U1D0D1D1[ijk]*kd[ijk==lmn]) + igConf_U0U1[ijk]*
(JJB0_D0D1(j_JJB0_D0D1,ijk,lmn) + dChrisConf_U1D1D1D0[ijk]*
kd[ijk==lmn]) + igConf_U0U2[ijk]*(JJB0_D0D2(j_JJB0_D0D2,ijk,lmn) + 
dChrisConf_U1D0D1D2[ijk]*kd[ijk==lmn]) + igConf_U0U2[ijk]*
(JJB0_D0D2(j_JJB0_D0D2,ijk,lmn) + dChrisConf_U1D1D2D0[ijk]*
kd[ijk==lmn]) + igConf_U1U1[ijk]*(JJB0_D1D1(j_JJB0_D1D1,ijk,lmn) + 
dChrisConf_U1D1D1D1[ijk]*kd[ijk==lmn]) + igConf_U1U2[ijk]*
(JJB0_D1D2(j_JJB0_D1D2,ijk,lmn) + dChrisConf_U1D1D1D2[ijk]*
kd[ijk==lmn]) + igConf_U1U2[ijk]*(JJB0_D1D2(j_JJB0_D1D2,ijk,lmn) + 
dChrisConf_U1D1D2D1[ijk]*kd[ijk==lmn]) + igConf_U2U2[ijk]*
(JJB0_D2D2(j_JJB0_D2D2,ijk,lmn) + dChrisConf_U1D1D2D2[ijk]*
kd[ijk==lmn]);


  double t2_E_U1 = 
2.0*ChrisConf_U1D0D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn)*igConf_U0U0[ijk] + 
2.0*ChrisConf_U1D1D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn)*igConf_U1U1[ijk] + 
2.0*ChrisConf_U1D1D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn)*igConf_U2U2[ijk] + 
2.0*igConf_U0U1[ijk]*(ChrisConf_U1D0D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U1D1D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn)) + 2.0*igConf_U0U2[ijk]*
(ChrisConf_U1D0D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + ChrisConf_U1D1D2[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn)) + 2.0*igConf_U1U2[ijk]*(ChrisConf_U1D1D1[ijk]*
JB0_D2(j_JB0_D2,ijk,lmn) + ChrisConf_U1D1D2[ijk]*JB0_D1(j_JB0_D1,ijk,lmn));





  double t3_E_U1 = 
kd[ijk==lmn]*(ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D0[ijk]*
igConf_U0U0[ijk] + ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk]*
igConf_U0U1[ijk] + ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D2[ijk]*
igConf_U0U2[ijk] + ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D0[ijk]*
igConf_U0U1[ijk] + ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D1[ijk]*
igConf_U1U1[ijk] + ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D2[ijk]*
igConf_U1U2[ijk] + ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D0[ijk]*
igConf_U0U2[ijk] + ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D1[ijk]*
igConf_U1U2[ijk] + ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D2[ijk]*
igConf_U2U2[ijk] + pow(ChrisConf_U1D0D1[ijk], 2)*igConf_U0U0[ijk] + 2.0*
ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D1[ijk]*igConf_U0U1[ijk] + 2.0*
ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D2[ijk]*igConf_U0U2[ijk] + 
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D0D1[ijk]*igConf_U0U0[ijk] + 
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D1[ijk]*igConf_U0U1[ijk] + 
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D2[ijk]*igConf_U0U2[ijk] + 
pow(ChrisConf_U1D1D1[ijk], 2)*igConf_U1U1[ijk] + 2.0*
ChrisConf_U1D1D1[ijk]*ChrisConf_U1D1D2[ijk]*igConf_U1U2[ijk] + 
pow(ChrisConf_U1D1D2[ijk], 2)*igConf_U2U2[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U1[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D1D1[ijk]*igConf_U1U1[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D1D2[ijk]*igConf_U1U2[ijk] + ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U2[ijk] + ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D1D1[ijk]*igConf_U1U2[ijk] + ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D1D2[ijk]*igConf_U2U2[ijk]);


  double t4_E_U1 = 
-igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + 
ChrisConf_U1D0D0[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + ChrisConf_U2D0D0[ijk]*
JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*(ChrisConf_U0D0D0[ijk]*
ChrisConf_U1D0D1[ijk] + ChrisConf_U1D0D0[ijk]*ChrisConf_U1D1D1[ijk] + 
ChrisConf_U1D1D2[ijk]*ChrisConf_U2D0D0[ijk])) - 2.0*igConf_U0U1[ijk]*
(ChrisConf_U0D0D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D0D1[ijk]*
JB0_D1(j_JB0_D1,ijk,lmn) + ChrisConf_U2D0D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + 
kd[ijk==lmn]*(ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk] + 
ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D1[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D0D1[ijk])) - 2.0*igConf_U0U2[ijk]*(ChrisConf_U0D0D2[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D0D2[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D0D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*
(ChrisConf_U0D0D2[ijk]*ChrisConf_U1D0D1[ijk] + ChrisConf_U1D0D2[ijk]*
ChrisConf_U1D1D1[ijk] + ChrisConf_U1D1D2[ijk]*ChrisConf_U2D0D2[ijk])) - 
igConf_U1U1[ijk]*(ChrisConf_U0D1D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + 
ChrisConf_U1D1D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + ChrisConf_U2D1D1[ijk]*
JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*(ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D1[ijk] + pow(ChrisConf_U1D1D1[ijk], 2) + 
ChrisConf_U1D1D2[ijk]*ChrisConf_U2D1D1[ijk])) - 2.0*igConf_U1U2[ijk]*
(ChrisConf_U0D1D2[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D1D2[ijk]*
JB0_D1(j_JB0_D1,ijk,lmn) + ChrisConf_U2D1D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + 
kd[ijk==lmn]*(ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D1[ijk] + 
ChrisConf_U1D1D1[ijk]*ChrisConf_U1D1D2[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D1D2[ijk])) - igConf_U2U2[ijk]*(ChrisConf_U0D2D2[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D2D2[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D2D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*
(ChrisConf_U0D2D2[ijk]*ChrisConf_U1D0D1[ijk] + ChrisConf_U1D1D1[ijk]*
ChrisConf_U1D2D2[ijk] + ChrisConf_U1D1D2[ijk]*ChrisConf_U2D2D2[ijk]));



  double t5_E_U1 = 
0.33333333333333331*igConf_U0U1[ijk]*(JJB0_D0D1(j_JJB0_D0D1,ijk,lmn) + 
kd[ijk==lmn]*(dChrisConf_U0D0D1D0[ijk] + dChrisConf_U1D1D1D0[ijk] + 
dChrisConf_U2D1D2D0[ijk])) + 0.33333333333333331*igConf_U1U1[ijk]*
(JJB0_D1D1(j_JJB0_D1D1,ijk,lmn) + kd[ijk==lmn]*(dChrisConf_U0D0D1D1[ijk] + 
dChrisConf_U1D1D1D1[ijk] + dChrisConf_U2D1D2D1[ijk])) + 
0.33333333333333331*igConf_U1U2[ijk]*(JJB0_D1D2(j_JJB0_D1D2,ijk,lmn) + 
kd[ijk==lmn]*(dChrisConf_U0D0D1D2[ijk] + dChrisConf_U1D1D1D2[ijk] + 
dChrisConf_U2D1D2D2[ijk]));


  double t6_E_U1 = 
0.33333333333333331*igConf_U0U1[ijk]*(2.0*ChrisConf_U0D0D1[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D0D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U1D1D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U2D0D1[ijk]*
JB0_D2(j_JB0_D2,ijk,lmn) + ChrisConf_U2D1D2[ijk]*JB0_D0(j_JB0_D0,ijk,lmn)) + 
0.33333333333333331*igConf_U1U1[ijk]*(ChrisConf_U0D0D1[ijk]*
JB0_D1(j_JB0_D1,ijk,lmn) + ChrisConf_U0D1D1[ijk]*JB0_D0(j_JB0_D0,ijk,lmn) + 
2.0*ChrisConf_U1D1D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D1D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + ChrisConf_U2D1D2[ijk]*
JB0_D1(j_JB0_D1,ijk,lmn)) + 0.33333333333333331*igConf_U1U2[ijk]*
(ChrisConf_U0D0D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + ChrisConf_U0D1D2[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D1D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + 
ChrisConf_U1D1D2[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 2.0*
ChrisConf_U2D1D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn));





  double t7_E_U1 = 
kd[ijk==lmn]*(0.33333333333333331*ChrisConf_U0D0D0[ijk]*
ChrisConf_U0D0D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
pow(ChrisConf_U0D0D1[ijk], 2)*igConf_U1U1[ijk] + 0.33333333333333331*
ChrisConf_U0D0D1[ijk]*ChrisConf_U0D0D2[ijk]*igConf_U1U2[ijk] + 
0.33333333333333331*ChrisConf_U0D0D1[ijk]*ChrisConf_U1D0D1[ijk]*
igConf_U0U1[ijk] + 0.33333333333333331*ChrisConf_U0D0D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D0[ijk]*igConf_U0U1[ijk] + 
0.66666666666666663*ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D1[ijk]*
igConf_U1U1[ijk] + 0.33333333333333331*ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D2[ijk]*igConf_U1U2[ijk] + 0.33333333333333331*
ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D1[ijk]*igConf_U1U2[ijk] + 
0.33333333333333331*ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D0[ijk]*
igConf_U0U1[ijk] + 0.66666666666666663*ChrisConf_U0D1D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U1U1[ijk] + 0.33333333333333331*
ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D2[ijk]*igConf_U1U2[ijk] + 
0.33333333333333331*ChrisConf_U0D2D2[ijk]*ChrisConf_U2D0D1[ijk]*
igConf_U1U2[ijk] + 0.33333333333333331*ChrisConf_U1D0D1[ijk]*
ChrisConf_U1D1D1[ijk]*igConf_U0U1[ijk] + 0.33333333333333331*
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D1[ijk]*igConf_U0U1[ijk] + 
0.33333333333333331*pow(ChrisConf_U1D1D1[ijk], 2)*igConf_U1U1[ijk] + 
0.33333333333333331*ChrisConf_U1D1D1[ijk]*ChrisConf_U1D1D2[ijk]*
igConf_U1U2[ijk] + 0.33333333333333331*ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D0D1[ijk]*igConf_U0U1[ijk] + 0.66666666666666663*
ChrisConf_U1D1D2[ijk]*ChrisConf_U2D1D1[ijk]*igConf_U1U1[ijk] + 
0.33333333333333331*ChrisConf_U1D1D2[ijk]*ChrisConf_U2D1D2[ijk]*
igConf_U1U2[ijk] + 0.33333333333333331*ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D1D1[ijk]*igConf_U1U2[ijk] + 0.33333333333333331*
ChrisConf_U2D0D2[ijk]*ChrisConf_U2D1D2[ijk]*igConf_U0U1[ijk] + 
0.33333333333333331*pow(ChrisConf_U2D1D2[ijk], 2)*igConf_U1U1[ijk] + 
0.33333333333333331*ChrisConf_U2D1D2[ijk]*ChrisConf_U2D2D2[ijk]*
igConf_U1U2[ijk]);


  double t8_E_U1 = 
-0.33333333333333331*igConf_U0U1[ijk]*(ChrisConf_U0D0D1[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D0D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D0D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*
(ChrisConf_U0D0D0[ijk]*ChrisConf_U0D0D1[ijk] + ChrisConf_U0D0D1[ijk]*
ChrisConf_U1D0D1[ijk] + ChrisConf_U0D0D2[ijk]*ChrisConf_U2D0D1[ijk] + 
ChrisConf_U0D1D1[ijk]*ChrisConf_U1D0D0[ijk] + ChrisConf_U0D1D2[ijk]*
ChrisConf_U2D0D0[ijk] + ChrisConf_U1D0D1[ijk]*ChrisConf_U1D1D1[ijk] + 
ChrisConf_U1D0D2[ijk]*ChrisConf_U2D1D1[ijk] + ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D0D1[ijk] + ChrisConf_U2D0D2[ijk]*ChrisConf_U2D1D2[ijk])) - 
0.33333333333333331*igConf_U1U1[ijk]*(ChrisConf_U0D1D1[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D1D1[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D1D1[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*
(pow(ChrisConf_U0D0D1[ijk], 2) + 2.0*ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D1[ijk] + 2.0*ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D1[ijk] + 
pow(ChrisConf_U1D1D1[ijk], 2) + 2.0*ChrisConf_U1D1D2[ijk]*
ChrisConf_U2D1D1[ijk] + pow(ChrisConf_U2D1D2[ijk], 2))) - 
0.33333333333333331*igConf_U1U2[ijk]*(ChrisConf_U0D1D2[ijk]*
JB0_D0(j_JB0_D0,ijk,lmn) + ChrisConf_U1D1D2[ijk]*JB0_D1(j_JB0_D1,ijk,lmn) + 
ChrisConf_U2D1D2[ijk]*JB0_D2(j_JB0_D2,ijk,lmn) + kd[ijk==lmn]*
(ChrisConf_U0D0D1[ijk]*ChrisConf_U0D0D2[ijk] + ChrisConf_U0D1D1[ijk]*
ChrisConf_U1D0D2[ijk] + ChrisConf_U0D1D2[ijk]*ChrisConf_U1D0D1[ijk] + 
ChrisConf_U0D1D2[ijk]*ChrisConf_U2D0D2[ijk] + ChrisConf_U0D2D2[ijk]*
ChrisConf_U2D0D1[ijk] + ChrisConf_U1D1D1[ijk]*ChrisConf_U1D1D2[ijk] + 
ChrisConf_U1D1D2[ijk]*ChrisConf_U2D1D2[ijk] + ChrisConf_U1D2D2[ijk]*
ChrisConf_U2D1D1[ijk] + ChrisConf_U2D1D2[ijk]*ChrisConf_U2D2D2[ijk]));



  double t9_E_U1 = 
kd[ijk==lmn]*(RicciConf_D0D1[ijk]*igConf_U0U1[ijk] + 
RicciConf_D1D1[ijk]*igConf_U1U1[ijk] + RicciConf_D1D2[ijk]*
igConf_U1U2[ijk]);


  double contE = 
JB0_D1(j_JB0_D1,ijk,lmn) + kd[ijk==lmn]*(ChrisConf_U0D0D1[ijk] +
ChrisConf_U1D1D1[ijk] + ChrisConf_U2D1D2[ijk]);



  double t1p10_E_U1 = 
0.66666666666666663*contE*(dLnOf_alpha_E_U0*igConf_U0U1[ijk] + 
dLnOf_alpha_E_U1*igConf_U1U1[ijk] + dLnOf_alpha_E_U2*
igConf_U1U2[ijk]);

  double t2p10_E_U1 = 
-dLnOf_alpha_E_U0*igConf_U0U0[ijk]*(ChrisConf_U1D0D1[ijk]*kd[ijk==lmn] + 
JB0_D0(j_JB0_D0,ijk,lmn)) - dLnOf_alpha_E_U0*igConf_U0U1[ijk]*
(ChrisConf_U1D1D1[ijk]*kd[ijk==lmn] + JB0_D1(j_JB0_D1,ijk,lmn)) - 
dLnOf_alpha_E_U0*igConf_U0U2[ijk]*(ChrisConf_U1D1D2[ijk]*kd[ijk==lmn] + 
JB0_D2(j_JB0_D2,ijk,lmn)) - dLnOf_alpha_E_U1*igConf_U0U1[ijk]*
(ChrisConf_U1D0D1[ijk]*kd[ijk==lmn] + JB0_D0(j_JB0_D0,ijk,lmn)) - 
dLnOf_alpha_E_U1*igConf_U1U1[ijk]*(ChrisConf_U1D1D1[ijk]*kd[ijk==lmn] + 
JB0_D1(j_JB0_D1,ijk,lmn)) - dLnOf_alpha_E_U1*igConf_U1U2[ijk]*
(ChrisConf_U1D1D2[ijk]*kd[ijk==lmn] + JB0_D2(j_JB0_D2,ijk,lmn)) - 
dLnOf_alpha_E_U2*igConf_U0U2[ijk]*(ChrisConf_U1D0D1[ijk]*kd[ijk==lmn] + 
JB0_D0(j_JB0_D0,ijk,lmn)) - dLnOf_alpha_E_U2*igConf_U1U2[ijk]*
(ChrisConf_U1D1D1[ijk]*kd[ijk==lmn] + JB0_D1(j_JB0_D1,ijk,lmn)) - 
dLnOf_alpha_E_U2*igConf_U2U2[ijk]*(ChrisConf_U1D1D2[ijk]*kd[ijk==lmn] + 
JB0_D2(j_JB0_D2,ijk,lmn));




  double t3p10_E_U1 = 
-ChrisConf_U0D0D1[ijk]*dLnOf_alpha_E_U0*igConf_U0U1[ijk]*kd[ijk==lmn] - 
ChrisConf_U0D1D1[ijk]*dLnOf_alpha_E_U0*igConf_U1U1[ijk]*kd[ijk==lmn] - 
ChrisConf_U0D1D2[ijk]*dLnOf_alpha_E_U0*igConf_U1U2[ijk]*kd[ijk==lmn] - 
ChrisConf_U2D0D1[ijk]*dLnOf_alpha_E_U2*igConf_U0U1[ijk]*kd[ijk==lmn] - 
ChrisConf_U2D1D1[ijk]*dLnOf_alpha_E_U2*igConf_U1U1[ijk]*kd[ijk==lmn] - 
ChrisConf_U2D1D2[ijk]*dLnOf_alpha_E_U2*igConf_U1U2[ijk]*kd[ijk==lmn] - 
dLnOf_alpha_E_U1*igConf_U0U1[ijk]*(ChrisConf_U1D0D1[ijk]*kd[ijk==lmn] + 
JB0_D0(j_JB0_D0,ijk,lmn)) - dLnOf_alpha_E_U1*igConf_U1U1[ijk]*
(ChrisConf_U1D1D1[ijk]*kd[ijk==lmn] + JB0_D1(j_JB0_D1,ijk,lmn)) - 
dLnOf_alpha_E_U1*igConf_U1U2[ijk]*(ChrisConf_U1D1D2[ijk]*kd[ijk==lmn] + 
JB0_D2(j_JB0_D2,ijk,lmn));


  double t10_E_U1 = 
t1p10_E_U1 + t2p10_E_U1 + t3p10_E_U1;





  double Epart_U1 = 
t10_E_U1 + t1_E_U1 + t2_E_U1 + t3_E_U1 + t4_E_U1 + t5_E_U1 + t6_E_U1 + 
t7_E_U1 + t8_E_U1 + t9_E_U1;

  E_Trans[j][i] = Epart_U1;

  DDM_SCHUR_JACOBIAN_EQ_Epart_CLOSE

  return 0;
}