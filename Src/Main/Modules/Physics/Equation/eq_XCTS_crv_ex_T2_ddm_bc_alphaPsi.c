/*
  These C codes generated by Cpi version 3.0
  Copyright (C) 2019-2023 Alireza Rashti.
*/


#include "eq_header.h"
#include "maths_equation_solvings_lib.h"

void *eq_XCTS_curve_exc_T2_ddm_bc_alphaPsi(void *vp1,void *vp2);
void *eq_XCTS_curve_exc_T2_ddm_bc_alphaPsi(void *vp1,void *vp2)
{
  DDM_SCHUR_BC_DECLARE
  Uint ijk;/* node index */
  READ_v(alphaPsi)
  READ_v(dalphaPsi_D0)
  READ_v(dalphaPsi_D1)
  READ_v(dalphaPsi_D2)
  if (patch->outerB)/* at outer boundary */
  {
  DDM_SCHUR_BC_OPEN

  double outerB_F =
alphaPsi[ijk] - 1;

  schur_F[schur_ijk] = outerB_F;

  DDM_SCHUR_BC_CLOSE
  }/* end of if (patch->outerB) */
  else if (patch->innerB)/* at inner boundary */
  {
  READ_v(bh_sConf_U0)
  READ_v(bh_sConf_U1)
  READ_v(bh_sConf_U2)
  DDM_SCHUR_BC_OPEN


  double innerB_F =
bh_sConf_U0[ijk]*dalphaPsi_D0[ijk] + bh_sConf_U1[ijk]*
dalphaPsi_D1[ijk] + bh_sConf_U2[ijk]*dalphaPsi_D2[ijk];

  schur_F[schur_ijk] = innerB_F;

  DDM_SCHUR_BC_CLOSE
  }/* end of else if (patch->innerB) */
  return 0;
}
