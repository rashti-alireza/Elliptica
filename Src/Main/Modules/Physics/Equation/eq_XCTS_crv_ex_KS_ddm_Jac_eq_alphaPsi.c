/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2021 Alireza Rashti.
*/


#include "eq_header.h"
#include "maths_equation_solvings_lib.h"

void *eq_XCTS_curve_exc_KS_ddm_jacobian_eq_alphaPsi(void *vp1,void *vp2);
void *eq_XCTS_curve_exc_KS_ddm_jacobian_eq_alphaPsi(void *vp1,void *vp2)
{
  DDM_SCHUR_JACOBIAN_EQ_DECLARE
  Uint ijk,lmn;/* for Jacobian entries J[ijk][lmn] */
  const double kd[2] = {0.,1.};/* Kronecker delta */

  /* declaring: */
  JACOBIAN_DERIVATIVE(JalphaPsi_D1)
  JACOBIAN_DERIVATIVE(JalphaPsi_D0)
  JACOBIAN_DERIVATIVE(JalphaPsi_D2)
  JACOBIAN_DERIVATIVE(JJalphaPsi_D0D0)
  JACOBIAN_DERIVATIVE(JJalphaPsi_D0D1)
  JACOBIAN_DERIVATIVE(JJalphaPsi_D0D2)
  JACOBIAN_DERIVATIVE(JJalphaPsi_D1D2)
  JACOBIAN_DERIVATIVE(JJalphaPsi_D2D2)
  JACOBIAN_DERIVATIVE(JJalphaPsi_D1D1)
  READ_v(psi)
  READ_v(trRicciConf)
  READ_v(trK)
  READ_v(EConf)
  READ_v(SConf)
  READ_v(AConfIJ2)
  READ_v(igConf_U2U2)
  READ_v(igConf_U1U2)
  READ_v(igConf_U1U1)
  READ_v(igConf_U0U2)
  READ_v(igConf_U0U0)
  READ_v(igConf_U0U1)
  READ_v(ChrisConf_U1D0D0)
  READ_v(ChrisConf_U0D2D2)
  READ_v(ChrisConf_U1D0D2)
  READ_v(ChrisConf_U2D2D2)
  READ_v(ChrisConf_U2D1D1)
  READ_v(ChrisConf_U2D0D0)
  READ_v(ChrisConf_U1D1D1)
  READ_v(ChrisConf_U0D1D1)
  READ_v(ChrisConf_U0D1D2)
  READ_v(ChrisConf_U1D1D2)
  READ_v(ChrisConf_U0D0D1)
  READ_v(ChrisConf_U0D0D0)
  READ_v(ChrisConf_U2D0D1)
  READ_v(ChrisConf_U0D0D2)
  READ_v(ChrisConf_U2D1D2)
  READ_v(ChrisConf_U1D2D2)
  READ_v(ChrisConf_U2D0D2)
  READ_v(ChrisConf_U1D0D1)


  DDM_SCHUR_JACOBIAN_EQ_Bpart_OPEN

  double psi4_B = 
pow(psi[ijk], 4);

  double aij2_B = 
AConfIJ2[ijk]/psi4_B;

  double Bpart = 
-igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D0D0[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D0D0[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D0D0(j_JJalphaPsi_D0D0,ijk,lmn)) - 2.0*igConf_U0U1[ijk]*
(ChrisConf_U0D0D1[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D0D1[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D0D1[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D0D1(j_JJalphaPsi_D0D1,ijk,lmn)) - 2.0*igConf_U0U2[ijk]*
(ChrisConf_U0D0D2[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D0D2[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D0D2[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D0D2(j_JJalphaPsi_D0D2,ijk,lmn)) - igConf_U1U1[ijk]*
(ChrisConf_U0D1D1[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D1D1[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D1D1[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D1D1(j_JJalphaPsi_D1D1,ijk,lmn)) - 2.0*igConf_U1U2[ijk]*
(ChrisConf_U0D1D2[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D1D2[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D1D2[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D1D2(j_JJalphaPsi_D1D2,ijk,lmn)) - igConf_U2U2[ijk]*
(ChrisConf_U0D2D2[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D2D2[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D2D2[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D2D2(j_JJalphaPsi_D2D2,ijk,lmn)) - kd[ijk==lmn]*(0.875*
aij2_B/psi4_B + 0.41666666666666669*psi4_B*pow(trK[ijk], 2) + (1.0/8.0)*
trRicciConf[ijk]) - 2*M_PI*kd[ijk==lmn]*(EConf[ijk] + 2*SConf[ijk])/
pow(psi[ijk], 2);

  B[i][j] = Bpart;

  DDM_SCHUR_JACOBIAN_EQ_Bpart_CLOSE

  DDM_SCHUR_JACOBIAN_EQ_Epart_OPEN

  double psi4_E = 
pow(psi[ijk], 4);

  double aij2_E = 
AConfIJ2[ijk]/psi4_E;

  double Epart = 
-igConf_U0U0[ijk]*(ChrisConf_U0D0D0[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D0D0[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D0D0[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D0D0(j_JJalphaPsi_D0D0,ijk,lmn)) - 2.0*igConf_U0U1[ijk]*
(ChrisConf_U0D0D1[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D0D1[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D0D1[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D0D1(j_JJalphaPsi_D0D1,ijk,lmn)) - 2.0*igConf_U0U2[ijk]*
(ChrisConf_U0D0D2[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D0D2[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D0D2[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D0D2(j_JJalphaPsi_D0D2,ijk,lmn)) - igConf_U1U1[ijk]*
(ChrisConf_U0D1D1[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D1D1[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D1D1[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D1D1(j_JJalphaPsi_D1D1,ijk,lmn)) - 2.0*igConf_U1U2[ijk]*
(ChrisConf_U0D1D2[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D1D2[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D1D2[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D1D2(j_JJalphaPsi_D1D2,ijk,lmn)) - igConf_U2U2[ijk]*
(ChrisConf_U0D2D2[ijk]*JalphaPsi_D0(j_JalphaPsi_D0,ijk,lmn) +
ChrisConf_U1D2D2[ijk]*JalphaPsi_D1(j_JalphaPsi_D1,ijk,lmn) +
ChrisConf_U2D2D2[ijk]*JalphaPsi_D2(j_JalphaPsi_D2,ijk,lmn) -
JJalphaPsi_D2D2(j_JJalphaPsi_D2D2,ijk,lmn)) - kd[ijk==lmn]*(0.875*
aij2_E/psi4_E + 0.41666666666666669*psi4_E*pow(trK[ijk], 2) + (1.0/8.0)*
trRicciConf[ijk]) - 2*M_PI*kd[ijk==lmn]*(EConf[ijk] + 2*SConf[ijk])/
pow(psi[ijk], 2);

  E_Trans[j][i] = Epart;

  DDM_SCHUR_JACOBIAN_EQ_Epart_CLOSE

  return 0;
}