/*
  These C codes generated by Cpi version 1.0
  Copyright (C) 2019 Alireza Rashti.
*/


#include "bbn_headers.h"

#define add_and_get_field(name) \
  if (_Ind(#name) >= 0)\
  {DECLARE_FIELD(name);REMOVE_FIELD(name);}\
  ADD_FIELD(name);REALLOC_v_WRITE_v(name);


void bbn_Rc_NS(double Rc[3],Grid_T *const grid);
void bbn_Rc_NS(double Rc[3],Grid_T *const grid)
{
  const double Madm = Pgetd("NS_ADM_mass");
  const double x_CM = Pgetd("x_CM");
  const double y_CM = Pgetd("y_CM");
  unsigned p;

  Rc[0] = 0;
  Rc[1] = 0;
  Rc[2] = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if (!IsItNSPatch(patch))
      continue;


  /* declaring: */
  READ_v(dpsi_D0)
  READ_v(dpsi_D1)
  READ_v(dpsi_D2)
  READ_v(ddpsi_D2D2)
  READ_v(ddpsi_D0D0)
  READ_v(ddpsi_D0D1)
  READ_v(ddpsi_D1D1)
  READ_v(ddpsi_D1D2)
  READ_v(ddpsi_D0D2)
  READ_v(_gamma_D2D2)
  READ_v(_gamma_D0D2)
  READ_v(_gamma_D0D0)
  READ_v(_gamma_D0D1)
  READ_v(_gamma_D1D2)
  READ_v(_gamma_D1D1)
  READ_v(_gammaI_U0U2)
  READ_v(_gammaI_U0U0)
  READ_v(_gammaI_U0U1)
  READ_v(_gammaI_U1U2)
  READ_v(_gammaI_U1U1)
  READ_v(_gammaI_U2U2)
  READ_v(_Gamma_U2D1D1)
  READ_v(_Gamma_U2D1D2)
  READ_v(_Gamma_U0D1D1)
  READ_v(_Gamma_U2D0D2)
  READ_v(_Gamma_U2D2D2)
  READ_v(_Gamma_U0D1D2)
  READ_v(_Gamma_U0D0D2)
  READ_v(_Gamma_U0D0D1)
  READ_v(_Gamma_U0D0D0)
  READ_v(_Gamma_U1D2D2)
  READ_v(_Gamma_U2D0D1)
  READ_v(_Gamma_U0D2D2)
  READ_v(_Gamma_U2D0D0)
  READ_v(_Gamma_U1D0D2)
  READ_v(_Gamma_U1D1D2)
  READ_v(_Gamma_U1D0D0)
  READ_v(_Gamma_U1D0D1)
  READ_v(_Gamma_U1D1D1)


    unsigned nn  = patch->nn;
    unsigned ijk;

    ADD_FIELD(Rc_integrandx)
    ADD_FIELD(Rc_integrandy)
    ADD_FIELD(Rc_integrandz)
    {
    REALLOC_v_WRITE_v(Rc_integrandx)
    REALLOC_v_WRITE_v(Rc_integrandy)
    REALLOC_v_WRITE_v(Rc_integrandz)
    for (ijk = 0; ijk < nn; ++ijk)
    {
    double D2psi = 
-_Gamma_U0D0D0[ijk]*_gammaI_U0U0[ijk]*dpsi_D0[ijk] - 2.0*
_Gamma_U0D0D1[ijk]*_gammaI_U0U1[ijk]*dpsi_D0[ijk] - 2.0*
_Gamma_U0D0D2[ijk]*_gammaI_U0U2[ijk]*dpsi_D0[ijk] - _Gamma_U0D1D1[ijk]*
_gammaI_U1U1[ijk]*dpsi_D0[ijk] - 2.0*_Gamma_U0D1D2[ijk]*
_gammaI_U1U2[ijk]*dpsi_D0[ijk] - _Gamma_U0D2D2[ijk]*_gammaI_U2U2[ijk]*
dpsi_D0[ijk] - _Gamma_U1D0D0[ijk]*_gammaI_U0U0[ijk]*dpsi_D1[ijk] - 2.0*
_Gamma_U1D0D1[ijk]*_gammaI_U0U1[ijk]*dpsi_D1[ijk] - 2.0*
_Gamma_U1D0D2[ijk]*_gammaI_U0U2[ijk]*dpsi_D1[ijk] - _Gamma_U1D1D1[ijk]*
_gammaI_U1U1[ijk]*dpsi_D1[ijk] - 2.0*_Gamma_U1D1D2[ijk]*
_gammaI_U1U2[ijk]*dpsi_D1[ijk] - _Gamma_U1D2D2[ijk]*_gammaI_U2U2[ijk]*
dpsi_D1[ijk] - _Gamma_U2D0D0[ijk]*_gammaI_U0U0[ijk]*dpsi_D2[ijk] - 2.0*
_Gamma_U2D0D1[ijk]*_gammaI_U0U1[ijk]*dpsi_D2[ijk] - 2.0*
_Gamma_U2D0D2[ijk]*_gammaI_U0U2[ijk]*dpsi_D2[ijk] - _Gamma_U2D1D1[ijk]*
_gammaI_U1U1[ijk]*dpsi_D2[ijk] - 2.0*_Gamma_U2D1D2[ijk]*
_gammaI_U1U2[ijk]*dpsi_D2[ijk] - _Gamma_U2D2D2[ijk]*_gammaI_U2U2[ijk]*
dpsi_D2[ijk] + _gammaI_U0U0[ijk]*ddpsi_D0D0[ijk] + 2.0*
_gammaI_U0U1[ijk]*ddpsi_D0D1[ijk] + 2.0*_gammaI_U0U2[ijk]*
ddpsi_D0D2[ijk] + _gammaI_U1U1[ijk]*ddpsi_D1D1[ijk] + 2.0*
_gammaI_U1U2[ijk]*ddpsi_D1D2[ijk] + _gammaI_U2U2[ijk]*
ddpsi_D2D2[ijk];

    double x = patch->node[ijk]->x[0];
    double y = patch->node[ijk]->x[1];
    double z = patch->node[ijk]->x[2];
    Rc_integrandx[ijk] = D2psi*(x-x_CM);
    Rc_integrandy[ijk] = D2psi*(y-y_CM);
    Rc_integrandz[ijk] = D2psi*(z);
    }
    }
    DECLARE_FIELD(Rc_integrandx)
    DECLARE_FIELD(Rc_integrandy)
    DECLARE_FIELD(Rc_integrandz)
    Integration_T *I = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->g00 = _gamma_D0D0;
    I->g01 = _gamma_D0D1;
    I->g02 = _gamma_D0D2;
    I->g11 = _gamma_D1D1;
    I->g12 = _gamma_D1D2;
    I->g22 = _gamma_D2D2;
    I->Spectral->f = Rc_integrandx;
    plan_integration(I);
    Rc[0] += execute_integration(I);
    I->Spectral->f = Rc_integrandy;
    plan_integration(I);
    Rc[1] += execute_integration(I);
    I->Spectral->f = Rc_integrandz;
    plan_integration(I);
    Rc[2] += execute_integration(I);
    free_integration(I);
    REMOVE_FIELD(Rc_integrandx)
    REMOVE_FIELD(Rc_integrandy)
    REMOVE_FIELD(Rc_integrandz)
  }

  Rc[0] /= (-2*M_PI*Madm);
  Rc[1] /= (-2*M_PI*Madm);
  Rc[2] /= (-2*M_PI*Madm);
}