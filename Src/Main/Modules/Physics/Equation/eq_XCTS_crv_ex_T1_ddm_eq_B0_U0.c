/*
  These C codes generated by Cpi version 3.0
  Copyright (C) 2019-2022 Alireza Rashti.
*/


#include "eq_header.h"
#include "maths_equation_solvings_lib.h"

void *eq_XCTS_curve_exc_T1_ddm_eq_B0_U0(void *vp1,void *vp2);
void *eq_XCTS_curve_exc_T1_ddm_eq_B0_U0(void *vp1,void *vp2)
{
  DDM_SCHUR_EQ_DECLARE
  Uint ijk;/* node index */
  READ_v_UNUSED(AConfIJ_U0U0)
  READ_v_UNUSED(AConfIJ_U0U1)
  READ_v_UNUSED(AConfIJ_U0U2)
  READ_v_UNUSED(AConfIJ_U1U1)
  READ_v_UNUSED(AConfIJ_U1U2)
  READ_v_UNUSED(AConfIJ_U2U2)
  READ_v_UNUSED(dAConfIJ_U0U0D0)
  READ_v_UNUSED(dAConfIJ_U0U0D1)
  READ_v_UNUSED(dAConfIJ_U0U0D2)
  READ_v_UNUSED(dAConfIJ_U0U1D0)
  READ_v_UNUSED(dAConfIJ_U0U1D1)
  READ_v_UNUSED(dAConfIJ_U0U1D2)
  READ_v_UNUSED(dAConfIJ_U0U2D0)
  READ_v_UNUSED(dAConfIJ_U0U2D1)
  READ_v_UNUSED(dAConfIJ_U0U2D2)
  READ_v_UNUSED(dAConfIJ_U1U1D0)
  READ_v_UNUSED(dAConfIJ_U1U1D1)
  READ_v_UNUSED(dAConfIJ_U1U1D2)
  READ_v_UNUSED(dAConfIJ_U1U2D0)
  READ_v_UNUSED(dAConfIJ_U1U2D1)
  READ_v_UNUSED(dAConfIJ_U1U2D2)
  READ_v_UNUSED(dAConfIJ_U2U2D0)
  READ_v_UNUSED(dAConfIJ_U2U2D1)
  READ_v_UNUSED(dAConfIJ_U2U2D2)
  READ_v(psi)
  READ_v(alphaPsi)
  READ_v(dtrK_D0)
  READ_v(dtrK_D1)
  READ_v(dtrK_D2)
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
  READ_v_UNUSED(JConfP_U0)
  READ_v_UNUSED(JConfP_U1)
  READ_v_UNUSED(JConfP_U2)
  READ_v(JConfC)
  DDM_SCHUR_EQ_OPEN

  double F00_U0 =
(2.0*AConfIJ_U0U0[ijk]*ChrisConf_U0D0D0[ijk] + AConfIJ_U0U0[ijk]*
ChrisConf_U1D0D1[ijk] + AConfIJ_U0U0[ijk]*ChrisConf_U2D0D2[ijk] + 3.0*
AConfIJ_U0U1[ijk]*ChrisConf_U0D0D1[ijk] + AConfIJ_U0U1[ijk]*
ChrisConf_U1D1D1[ijk] + AConfIJ_U0U1[ijk]*ChrisConf_U2D1D2[ijk] + 3.0*
AConfIJ_U0U2[ijk]*ChrisConf_U0D0D2[ijk] + AConfIJ_U0U2[ijk]*
ChrisConf_U1D1D2[ijk] + AConfIJ_U0U2[ijk]*ChrisConf_U2D2D2[ijk] +
AConfIJ_U1U1[ijk]*ChrisConf_U0D1D1[ijk] + 2.0*AConfIJ_U1U2[ijk]*
ChrisConf_U0D1D2[ijk] + AConfIJ_U2U2[ijk]*ChrisConf_U0D2D2[ijk] +
dAConfIJ_U0U0D0[ijk] + dAConfIJ_U0U1D1[ijk] + dAConfIJ_U0U2D2[ijk])/
pow(psi[ijk], 3);



  double F0_U0 =
F00_U0/pow(psi[ijk], 4);



  double F1_U0 =
0.66666666666666663*(-dtrK_D0[ijk]*igConf_U0U0[ijk] - dtrK_D1[ijk]*
igConf_U0U1[ijk] - dtrK_D2[ijk]*igConf_U0U2[ijk])/psi[ijk];



  double s_U0 =
-8*M_PI*JConfP_U0[ijk]/pow(psi[ijk], 3);



  double ell_U0 =
JConfC[ijk]*(F0_U0 + F1_U0);



  double F_eq_U0 =
2*alphaPsi[ijk]*(ell_U0 + s_U0);



  schur_F[schur_ijk] = F_eq_U0;

  DDM_SCHUR_EQ_CLOSE

  return 0;
}
