/*
  These C codes generated by Cpi version 1.0
  Copyright (C) 2019 Alireza Rashti.
*/


#include "sbh_headers.h"


void sbh_update_psi10A_UiUj(Patch_T *const patch)
{

  /* declaring: */
  PREP_FIELD(_A_UiUj_U2U2)
  PREP_FIELD(_A_UiUj_U1U2)
  PREP_FIELD(_A_UiUj_U1U1)
  PREP_FIELD(_A_UiUj_U0U2)
  PREP_FIELD(_A_UiUj_U0U1)
  PREP_FIELD(_A_UiUj_U0U0)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U2U2D2)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U2U2D0)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U2U2D1)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U1U1D2)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U1U1D0)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U1U1D1)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U0U0D2)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U0U0D0)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U0U0D1)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U0U1D2)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U0U1D1)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U0U1D0)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U1U2D1)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U1U2D0)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U1U2D2)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U0U2D0)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U0U2D1)
  DECLARE_AND_EMPTY_FIELD(_dA_UiUj_U0U2D2)
  PREP_FIELD(_Aij2)
  GET_FIELD(Beta_U1)
  GET_FIELD(Beta_U0)
  GET_FIELD(Beta_U2)
  GET_FIELD(dBeta_U2D2)
  GET_FIELD(dBeta_U2D0)
  GET_FIELD(dBeta_U2D1)
  GET_FIELD(dBeta_U0D0)
  GET_FIELD(dBeta_U1D1)
  GET_FIELD(dBeta_U0D1)
  GET_FIELD(dBeta_U0D2)
  GET_FIELD(dBeta_U1D2)
  GET_FIELD(dBeta_U1D0)
  GET_FIELD(psi)
  GET_FIELD(eta)
  GET_FIELD(_gamma_D2D2)
  GET_FIELD(_gamma_D0D2)
  GET_FIELD(_gamma_D0D0)
  GET_FIELD(_gamma_D0D1)
  GET_FIELD(_gamma_D1D2)
  GET_FIELD(_gamma_D1D1)
  GET_FIELD(_gammaI_U0U2)
  GET_FIELD(_gammaI_U0U0)
  GET_FIELD(_gammaI_U0U1)
  GET_FIELD(_gammaI_U1U2)
  GET_FIELD(_gammaI_U1U1)
  GET_FIELD(_gammaI_U2U2)
  GET_FIELD(_Gamma_U2D1D1)
  GET_FIELD(_Gamma_U2D1D2)
  GET_FIELD(_Gamma_U0D1D1)
  GET_FIELD(_Gamma_U2D0D2)
  GET_FIELD(_Gamma_U2D2D2)
  GET_FIELD(_Gamma_U0D1D2)
  GET_FIELD(_Gamma_U0D0D2)
  GET_FIELD(_Gamma_U0D0D1)
  GET_FIELD(_Gamma_U0D0D0)
  GET_FIELD(_Gamma_U1D2D2)
  GET_FIELD(_Gamma_U2D0D1)
  GET_FIELD(_Gamma_U0D2D2)
  GET_FIELD(_Gamma_U2D0D0)
  GET_FIELD(_Gamma_U1D0D2)
  GET_FIELD(_Gamma_U1D1D2)
  GET_FIELD(_Gamma_U1D0D0)
  GET_FIELD(_Gamma_U1D0D1)
  GET_FIELD(_Gamma_U1D1D1)


    unsigned nn = patch->nn;
    unsigned ijk;
    for(ijk = 0; ijk < nn; ++ijk)
    {
    double alphabar = 
eta[ijk]/pow(psi[ijk], 7);

    double _DB_UU_U1U2 = 
_gammaI_U0U1[ijk]*(Beta_U0[ijk]*_Gamma_U2D0D0[ijk] + Beta_U1[ijk]*
_Gamma_U2D0D1[ijk] + Beta_U2[ijk]*_Gamma_U2D0D2[ijk] + 
dBeta_U2D0[ijk]) + _gammaI_U1U1[ijk]*(Beta_U0[ijk]*_Gamma_U2D0D1[ijk] + 
Beta_U1[ijk]*_Gamma_U2D1D1[ijk] + Beta_U2[ijk]*_Gamma_U2D1D2[ijk] + 
dBeta_U2D1[ijk]) + _gammaI_U1U2[ijk]*(Beta_U0[ijk]*_Gamma_U2D0D2[ijk] + 
Beta_U1[ijk]*_Gamma_U2D1D2[ijk] + Beta_U2[ijk]*_Gamma_U2D2D2[ijk] + 
dBeta_U2D2[ijk]);

    double _DB_UU_U1U1 = 
_gammaI_U0U1[ijk]*(Beta_U0[ijk]*_Gamma_U1D0D0[ijk] + Beta_U1[ijk]*
_Gamma_U1D0D1[ijk] + Beta_U2[ijk]*_Gamma_U1D0D2[ijk] + 
dBeta_U1D0[ijk]) + _gammaI_U1U1[ijk]*(Beta_U0[ijk]*_Gamma_U1D0D1[ijk] + 
Beta_U1[ijk]*_Gamma_U1D1D1[ijk] + Beta_U2[ijk]*_Gamma_U1D1D2[ijk] + 
dBeta_U1D1[ijk]) + _gammaI_U1U2[ijk]*(Beta_U0[ijk]*_Gamma_U1D0D2[ijk] + 
Beta_U1[ijk]*_Gamma_U1D1D2[ijk] + Beta_U2[ijk]*_Gamma_U1D2D2[ijk] + 
dBeta_U1D2[ijk]);

    double _DB_UU_U1U0 = 
_gammaI_U0U1[ijk]*(Beta_U0[ijk]*_Gamma_U0D0D0[ijk] + Beta_U1[ijk]*
_Gamma_U0D0D1[ijk] + Beta_U2[ijk]*_Gamma_U0D0D2[ijk] + 
dBeta_U0D0[ijk]) + _gammaI_U1U1[ijk]*(Beta_U0[ijk]*_Gamma_U0D0D1[ijk] + 
Beta_U1[ijk]*_Gamma_U0D1D1[ijk] + Beta_U2[ijk]*_Gamma_U0D1D2[ijk] + 
dBeta_U0D1[ijk]) + _gammaI_U1U2[ijk]*(Beta_U0[ijk]*_Gamma_U0D0D2[ijk] + 
Beta_U1[ijk]*_Gamma_U0D1D2[ijk] + Beta_U2[ijk]*_Gamma_U0D2D2[ijk] + 
dBeta_U0D2[ijk]);

    double _DB_UU_U0U2 = 
_gammaI_U0U0[ijk]*(Beta_U0[ijk]*_Gamma_U2D0D0[ijk] + Beta_U1[ijk]*
_Gamma_U2D0D1[ijk] + Beta_U2[ijk]*_Gamma_U2D0D2[ijk] + 
dBeta_U2D0[ijk]) + _gammaI_U0U1[ijk]*(Beta_U0[ijk]*_Gamma_U2D0D1[ijk] + 
Beta_U1[ijk]*_Gamma_U2D1D1[ijk] + Beta_U2[ijk]*_Gamma_U2D1D2[ijk] + 
dBeta_U2D1[ijk]) + _gammaI_U0U2[ijk]*(Beta_U0[ijk]*_Gamma_U2D0D2[ijk] + 
Beta_U1[ijk]*_Gamma_U2D1D2[ijk] + Beta_U2[ijk]*_Gamma_U2D2D2[ijk] + 
dBeta_U2D2[ijk]);

    double _DB_UU_U0U0 = 
_gammaI_U0U0[ijk]*(Beta_U0[ijk]*_Gamma_U0D0D0[ijk] + Beta_U1[ijk]*
_Gamma_U0D0D1[ijk] + Beta_U2[ijk]*_Gamma_U0D0D2[ijk] + 
dBeta_U0D0[ijk]) + _gammaI_U0U1[ijk]*(Beta_U0[ijk]*_Gamma_U0D0D1[ijk] + 
Beta_U1[ijk]*_Gamma_U0D1D1[ijk] + Beta_U2[ijk]*_Gamma_U0D1D2[ijk] + 
dBeta_U0D1[ijk]) + _gammaI_U0U2[ijk]*(Beta_U0[ijk]*_Gamma_U0D0D2[ijk] + 
Beta_U1[ijk]*_Gamma_U0D1D2[ijk] + Beta_U2[ijk]*_Gamma_U0D2D2[ijk] + 
dBeta_U0D2[ijk]);

    double _DB_UU_U0U1 = 
_gammaI_U0U0[ijk]*(Beta_U0[ijk]*_Gamma_U1D0D0[ijk] + Beta_U1[ijk]*
_Gamma_U1D0D1[ijk] + Beta_U2[ijk]*_Gamma_U1D0D2[ijk] + 
dBeta_U1D0[ijk]) + _gammaI_U0U1[ijk]*(Beta_U0[ijk]*_Gamma_U1D0D1[ijk] + 
Beta_U1[ijk]*_Gamma_U1D1D1[ijk] + Beta_U2[ijk]*_Gamma_U1D1D2[ijk] + 
dBeta_U1D1[ijk]) + _gammaI_U0U2[ijk]*(Beta_U0[ijk]*_Gamma_U1D0D2[ijk] + 
Beta_U1[ijk]*_Gamma_U1D1D2[ijk] + Beta_U2[ijk]*_Gamma_U1D2D2[ijk] + 
dBeta_U1D2[ijk]);

    double _DB_UU_U2U0 = 
_gammaI_U0U2[ijk]*(Beta_U0[ijk]*_Gamma_U0D0D0[ijk] + Beta_U1[ijk]*
_Gamma_U0D0D1[ijk] + Beta_U2[ijk]*_Gamma_U0D0D2[ijk] + 
dBeta_U0D0[ijk]) + _gammaI_U1U2[ijk]*(Beta_U0[ijk]*_Gamma_U0D0D1[ijk] + 
Beta_U1[ijk]*_Gamma_U0D1D1[ijk] + Beta_U2[ijk]*_Gamma_U0D1D2[ijk] + 
dBeta_U0D1[ijk]) + _gammaI_U2U2[ijk]*(Beta_U0[ijk]*_Gamma_U0D0D2[ijk] + 
Beta_U1[ijk]*_Gamma_U0D1D2[ijk] + Beta_U2[ijk]*_Gamma_U0D2D2[ijk] + 
dBeta_U0D2[ijk]);

    double _DB_UU_U2U1 = 
_gammaI_U0U2[ijk]*(Beta_U0[ijk]*_Gamma_U1D0D0[ijk] + Beta_U1[ijk]*
_Gamma_U1D0D1[ijk] + Beta_U2[ijk]*_Gamma_U1D0D2[ijk] + 
dBeta_U1D0[ijk]) + _gammaI_U1U2[ijk]*(Beta_U0[ijk]*_Gamma_U1D0D1[ijk] + 
Beta_U1[ijk]*_Gamma_U1D1D1[ijk] + Beta_U2[ijk]*_Gamma_U1D1D2[ijk] + 
dBeta_U1D1[ijk]) + _gammaI_U2U2[ijk]*(Beta_U0[ijk]*_Gamma_U1D0D2[ijk] + 
Beta_U1[ijk]*_Gamma_U1D1D2[ijk] + Beta_U2[ijk]*_Gamma_U1D2D2[ijk] + 
dBeta_U1D2[ijk]);

    double _DB_UU_U2U2 = 
_gammaI_U0U2[ijk]*(Beta_U0[ijk]*_Gamma_U2D0D0[ijk] + Beta_U1[ijk]*
_Gamma_U2D0D1[ijk] + Beta_U2[ijk]*_Gamma_U2D0D2[ijk] + 
dBeta_U2D0[ijk]) + _gammaI_U1U2[ijk]*(Beta_U0[ijk]*_Gamma_U2D0D1[ijk] + 
Beta_U1[ijk]*_Gamma_U2D1D1[ijk] + Beta_U2[ijk]*_Gamma_U2D1D2[ijk] + 
dBeta_U2D1[ijk]) + _gammaI_U2U2[ijk]*(Beta_U0[ijk]*_Gamma_U2D0D2[ijk] + 
Beta_U1[ijk]*_Gamma_U2D1D2[ijk] + Beta_U2[ijk]*_Gamma_U2D2D2[ijk] + 
dBeta_U2D2[ijk]);

    double _LV_UU_U0U2 = 
-0.66666666666666663*_DB_UU_U0U0*_gammaI_U0U2[ijk]*_gamma_D0D0[ijk] - 
0.66666666666666663*_DB_UU_U0U1*_gammaI_U0U2[ijk]*_gamma_D0D1[ijk] - 
0.66666666666666663*_DB_UU_U0U2*_gammaI_U0U2[ijk]*_gamma_D0D2[ijk] + 
_DB_UU_U0U2 - 0.66666666666666663*_DB_UU_U1U0*_gammaI_U0U2[ijk]*
_gamma_D0D1[ijk] - 0.66666666666666663*_DB_UU_U1U1*_gammaI_U0U2[ijk]*
_gamma_D1D1[ijk] - 0.66666666666666663*_DB_UU_U1U2*_gammaI_U0U2[ijk]*
_gamma_D1D2[ijk] - 0.66666666666666663*_DB_UU_U2U0*_gammaI_U0U2[ijk]*
_gamma_D0D2[ijk] + _DB_UU_U2U0 - 0.66666666666666663*_DB_UU_U2U1*
_gammaI_U0U2[ijk]*_gamma_D1D2[ijk] - 0.66666666666666663*_DB_UU_U2U2*
_gammaI_U0U2[ijk]*_gamma_D2D2[ijk];

    double _LV_UU_U0U0 = 
-0.66666666666666663*_DB_UU_U0U0*_gammaI_U0U0[ijk]*_gamma_D0D0[ijk] + 
2*_DB_UU_U0U0 - 0.66666666666666663*_DB_UU_U0U1*_gammaI_U0U0[ijk]*
_gamma_D0D1[ijk] - 0.66666666666666663*_DB_UU_U0U2*_gammaI_U0U0[ijk]*
_gamma_D0D2[ijk] - 0.66666666666666663*_DB_UU_U1U0*_gammaI_U0U0[ijk]*
_gamma_D0D1[ijk] - 0.66666666666666663*_DB_UU_U1U1*_gammaI_U0U0[ijk]*
_gamma_D1D1[ijk] - 0.66666666666666663*_DB_UU_U1U2*_gammaI_U0U0[ijk]*
_gamma_D1D2[ijk] - 0.66666666666666663*_DB_UU_U2U0*_gammaI_U0U0[ijk]*
_gamma_D0D2[ijk] - 0.66666666666666663*_DB_UU_U2U1*_gammaI_U0U0[ijk]*
_gamma_D1D2[ijk] - 0.66666666666666663*_DB_UU_U2U2*_gammaI_U0U0[ijk]*
_gamma_D2D2[ijk];

    double _LV_UU_U0U1 = 
-0.66666666666666663*_DB_UU_U0U0*_gammaI_U0U1[ijk]*_gamma_D0D0[ijk] - 
0.66666666666666663*_DB_UU_U0U1*_gammaI_U0U1[ijk]*_gamma_D0D1[ijk] + 
_DB_UU_U0U1 - 0.66666666666666663*_DB_UU_U0U2*_gammaI_U0U1[ijk]*
_gamma_D0D2[ijk] - 0.66666666666666663*_DB_UU_U1U0*_gammaI_U0U1[ijk]*
_gamma_D0D1[ijk] + _DB_UU_U1U0 - 0.66666666666666663*_DB_UU_U1U1*
_gammaI_U0U1[ijk]*_gamma_D1D1[ijk] - 0.66666666666666663*_DB_UU_U1U2*
_gammaI_U0U1[ijk]*_gamma_D1D2[ijk] - 0.66666666666666663*_DB_UU_U2U0*
_gammaI_U0U1[ijk]*_gamma_D0D2[ijk] - 0.66666666666666663*_DB_UU_U2U1*
_gammaI_U0U1[ijk]*_gamma_D1D2[ijk] - 0.66666666666666663*_DB_UU_U2U2*
_gammaI_U0U1[ijk]*_gamma_D2D2[ijk];

    double _LV_UU_U2U2 = 
-0.66666666666666663*_DB_UU_U0U0*_gammaI_U2U2[ijk]*_gamma_D0D0[ijk] - 
0.66666666666666663*_DB_UU_U0U1*_gammaI_U2U2[ijk]*_gamma_D0D1[ijk] - 
0.66666666666666663*_DB_UU_U0U2*_gammaI_U2U2[ijk]*_gamma_D0D2[ijk] - 
0.66666666666666663*_DB_UU_U1U0*_gammaI_U2U2[ijk]*_gamma_D0D1[ijk] - 
0.66666666666666663*_DB_UU_U1U1*_gammaI_U2U2[ijk]*_gamma_D1D1[ijk] - 
0.66666666666666663*_DB_UU_U1U2*_gammaI_U2U2[ijk]*_gamma_D1D2[ijk] - 
0.66666666666666663*_DB_UU_U2U0*_gammaI_U2U2[ijk]*_gamma_D0D2[ijk] - 
0.66666666666666663*_DB_UU_U2U1*_gammaI_U2U2[ijk]*_gamma_D1D2[ijk] - 
0.66666666666666663*_DB_UU_U2U2*_gammaI_U2U2[ijk]*_gamma_D2D2[ijk] + 2*
_DB_UU_U2U2;

    double _LV_UU_U1U2 = 
-0.66666666666666663*_DB_UU_U0U0*_gammaI_U1U2[ijk]*_gamma_D0D0[ijk] - 
0.66666666666666663*_DB_UU_U0U1*_gammaI_U1U2[ijk]*_gamma_D0D1[ijk] - 
0.66666666666666663*_DB_UU_U0U2*_gammaI_U1U2[ijk]*_gamma_D0D2[ijk] - 
0.66666666666666663*_DB_UU_U1U0*_gammaI_U1U2[ijk]*_gamma_D0D1[ijk] - 
0.66666666666666663*_DB_UU_U1U1*_gammaI_U1U2[ijk]*_gamma_D1D1[ijk] - 
0.66666666666666663*_DB_UU_U1U2*_gammaI_U1U2[ijk]*_gamma_D1D2[ijk] + 
_DB_UU_U1U2 - 0.66666666666666663*_DB_UU_U2U0*_gammaI_U1U2[ijk]*
_gamma_D0D2[ijk] - 0.66666666666666663*_DB_UU_U2U1*_gammaI_U1U2[ijk]*
_gamma_D1D2[ijk] + _DB_UU_U2U1 - 0.66666666666666663*_DB_UU_U2U2*
_gammaI_U1U2[ijk]*_gamma_D2D2[ijk];

    double _LV_UU_U1U1 = 
-0.66666666666666663*_DB_UU_U0U0*_gammaI_U1U1[ijk]*_gamma_D0D0[ijk] - 
0.66666666666666663*_DB_UU_U0U1*_gammaI_U1U1[ijk]*_gamma_D0D1[ijk] - 
0.66666666666666663*_DB_UU_U0U2*_gammaI_U1U1[ijk]*_gamma_D0D2[ijk] - 
0.66666666666666663*_DB_UU_U1U0*_gammaI_U1U1[ijk]*_gamma_D0D1[ijk] - 
0.66666666666666663*_DB_UU_U1U1*_gammaI_U1U1[ijk]*_gamma_D1D1[ijk] + 2*
_DB_UU_U1U1 - 0.66666666666666663*_DB_UU_U1U2*_gammaI_U1U1[ijk]*
_gamma_D1D2[ijk] - 0.66666666666666663*_DB_UU_U2U0*_gammaI_U1U1[ijk]*
_gamma_D0D2[ijk] - 0.66666666666666663*_DB_UU_U2U1*_gammaI_U1U1[ijk]*
_gamma_D1D2[ijk] - 0.66666666666666663*_DB_UU_U2U2*_gammaI_U1U1[ijk]*
_gamma_D2D2[ijk];

    double _A_UUij_U2U2 = 
(1.0/2.0)*_LV_UU_U2U2/alphabar;

    double _A_UUij_U1U2 = 
(1.0/2.0)*_LV_UU_U1U2/alphabar;

    double _A_UUij_U1U1 = 
(1.0/2.0)*_LV_UU_U1U1/alphabar;

    double _A_UUij_U0U2 = 
(1.0/2.0)*_LV_UU_U0U2/alphabar;

    double _A_UUij_U0U1 = 
(1.0/2.0)*_LV_UU_U0U1/alphabar;

    double _A_UUij_U0U0 = 
(1.0/2.0)*_LV_UU_U0U0/alphabar;

    double _AijAij = 
pow(_A_UUij_U0U0, 2)*pow(_gamma_D0D0[ijk], 2) + 4.0*_A_UUij_U0U0*
_A_UUij_U0U1*_gamma_D0D0[ijk]*_gamma_D0D1[ijk] + 4.0*_A_UUij_U0U0*
_A_UUij_U0U2*_gamma_D0D0[ijk]*_gamma_D0D2[ijk] + 2.0*_A_UUij_U0U0*
_A_UUij_U1U1*pow(_gamma_D0D1[ijk], 2) + 4.0*_A_UUij_U0U0*_A_UUij_U1U2*
_gamma_D0D1[ijk]*_gamma_D0D2[ijk] + 2.0*_A_UUij_U0U0*_A_UUij_U2U2*
pow(_gamma_D0D2[ijk], 2) + 2.0*pow(_A_UUij_U0U1, 2)*_gamma_D0D0[ijk]*
_gamma_D1D1[ijk] + 2.0*pow(_A_UUij_U0U1, 2)*pow(_gamma_D0D1[ijk], 2) +
4.0*_A_UUij_U0U1*_A_UUij_U0U2*_gamma_D0D0[ijk]*_gamma_D1D2[ijk] + 4.0*
_A_UUij_U0U1*_A_UUij_U0U2*_gamma_D0D1[ijk]*_gamma_D0D2[ijk] + 4.0*
_A_UUij_U0U1*_A_UUij_U1U1*_gamma_D0D1[ijk]*_gamma_D1D1[ijk] + 4.0*
_A_UUij_U0U1*_A_UUij_U1U2*_gamma_D0D1[ijk]*_gamma_D1D2[ijk] + 4.0*
_A_UUij_U0U1*_A_UUij_U1U2*_gamma_D0D2[ijk]*_gamma_D1D1[ijk] + 4.0*
_A_UUij_U0U1*_A_UUij_U2U2*_gamma_D0D2[ijk]*_gamma_D1D2[ijk] + 2.0*
pow(_A_UUij_U0U2, 2)*_gamma_D0D0[ijk]*_gamma_D2D2[ijk] + 2.0*
pow(_A_UUij_U0U2, 2)*pow(_gamma_D0D2[ijk], 2) + 4.0*_A_UUij_U0U2*
_A_UUij_U1U1*_gamma_D0D1[ijk]*_gamma_D1D2[ijk] + 4.0*_A_UUij_U0U2*
_A_UUij_U1U2*_gamma_D0D1[ijk]*_gamma_D2D2[ijk] + 4.0*_A_UUij_U0U2*
_A_UUij_U1U2*_gamma_D0D2[ijk]*_gamma_D1D2[ijk] + 4.0*_A_UUij_U0U2*
_A_UUij_U2U2*_gamma_D0D2[ijk]*_gamma_D2D2[ijk] + pow(_A_UUij_U1U1, 2)*
pow(_gamma_D1D1[ijk], 2) + 4.0*_A_UUij_U1U1*_A_UUij_U1U2*
_gamma_D1D1[ijk]*_gamma_D1D2[ijk] + 2.0*_A_UUij_U1U1*_A_UUij_U2U2*
pow(_gamma_D1D2[ijk], 2) + 2.0*pow(_A_UUij_U1U2, 2)*_gamma_D1D1[ijk]*
_gamma_D2D2[ijk] + 2.0*pow(_A_UUij_U1U2, 2)*pow(_gamma_D1D2[ijk], 2) +
4.0*_A_UUij_U1U2*_A_UUij_U2U2*_gamma_D1D2[ijk]*_gamma_D2D2[ijk] +
pow(_A_UUij_U2U2, 2)*pow(_gamma_D2D2[ijk], 2);


    /* populating: */
    _A_UiUj_U2U2[ijk] = _A_UUij_U2U2;
    _A_UiUj_U1U2[ijk] = _A_UUij_U1U2;
    _A_UiUj_U1U1[ijk] = _A_UUij_U1U1;
    _A_UiUj_U0U2[ijk] = _A_UUij_U0U2;
    _A_UiUj_U0U1[ijk] = _A_UUij_U0U1;
    _A_UiUj_U0U0[ijk] = _A_UUij_U0U0;
    _Aij2[ijk] = _AijAij;
    }/*end of for(ijk = 0; ijk < nn; ++ijk)*/
  Field_T *f_A_UiUj_U2U2 = patch->pool[Ind("_A_UiUj_U2U2")];
  Field_T *f_A_UiUj_U1U2 = patch->pool[Ind("_A_UiUj_U1U2")];
  Field_T *f_A_UiUj_U1U1 = patch->pool[Ind("_A_UiUj_U1U1")];
  Field_T *f_A_UiUj_U0U2 = patch->pool[Ind("_A_UiUj_U0U2")];
  Field_T *f_A_UiUj_U0U1 = patch->pool[Ind("_A_UiUj_U0U1")];
  Field_T *f_A_UiUj_U0U0 = patch->pool[Ind("_A_UiUj_U0U0")];
  _dA_UiUj_U2U2D2->v = Partial_Derivative(f_A_UiUj_U2U2,"z");
  _dA_UiUj_U2U2D0->v = Partial_Derivative(f_A_UiUj_U2U2,"x");
  _dA_UiUj_U2U2D1->v = Partial_Derivative(f_A_UiUj_U2U2,"y");
  _dA_UiUj_U1U1D2->v = Partial_Derivative(f_A_UiUj_U1U1,"z");
  _dA_UiUj_U1U1D0->v = Partial_Derivative(f_A_UiUj_U1U1,"x");
  _dA_UiUj_U1U1D1->v = Partial_Derivative(f_A_UiUj_U1U1,"y");
  _dA_UiUj_U0U0D2->v = Partial_Derivative(f_A_UiUj_U0U0,"z");
  _dA_UiUj_U0U0D0->v = Partial_Derivative(f_A_UiUj_U0U0,"x");
  _dA_UiUj_U0U0D1->v = Partial_Derivative(f_A_UiUj_U0U0,"y");
  _dA_UiUj_U0U1D2->v = Partial_Derivative(f_A_UiUj_U0U1,"z");
  _dA_UiUj_U0U1D1->v = Partial_Derivative(f_A_UiUj_U0U1,"y");
  _dA_UiUj_U0U1D0->v = Partial_Derivative(f_A_UiUj_U0U1,"x");
  _dA_UiUj_U1U2D1->v = Partial_Derivative(f_A_UiUj_U1U2,"y");
  _dA_UiUj_U1U2D0->v = Partial_Derivative(f_A_UiUj_U1U2,"x");
  _dA_UiUj_U1U2D2->v = Partial_Derivative(f_A_UiUj_U1U2,"z");
  _dA_UiUj_U0U2D0->v = Partial_Derivative(f_A_UiUj_U0U2,"x");
  _dA_UiUj_U0U2D1->v = Partial_Derivative(f_A_UiUj_U0U2,"y");
  _dA_UiUj_U0U2D2->v = Partial_Derivative(f_A_UiUj_U0U2,"z");
}