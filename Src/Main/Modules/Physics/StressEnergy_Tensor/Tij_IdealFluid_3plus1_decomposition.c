/*
// Alireza Rashti
// August 2019
*/

#include "Tij_IdealFluid_3plus1_decomposition.h"

/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid, CTS for conformal thin sandwich method;
// thus, Tij_IF_CTS means stress enerfy of ideal fluid in CTS method. 
// given all of the fields needed it builds 
// the first component of fluid four velocity, i.e. 
// u_mu = (u_U0,u_U1,u_U2,u_U3). */
void Tij_IF_CTS_u0(Patch_T *const patch)
{
  if (!IsItNSPatch(patch))
    return;
  const unsigned nn = patch->nn;
  unsigned ijk;

  /* declaring: */
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
  GET_FIELD(enthalpy)
  GET_FIELD(W_U1)
  GET_FIELD(W_U0)
  GET_FIELD(W_U2)
  GET_FIELD(dphi_D2)
  GET_FIELD(dphi_D1)
  GET_FIELD(dphi_D0)
  GET_FIELD(eta)
  GET_FIELD(psi)
  GET_FIELD(u0)


  for(ijk = 0; ijk < nn; ++ijk)
  {
  double alpha = 
eta[ijk]/psi[ijk];

  double psim4 = 
pow(psi[ijk], -4);

  double psi4 = 
pow(psi[ijk], 4);

  double P2 = 
2.0*W_U0[ijk]*dphi_D0[ijk] + 2.0*W_U1[ijk]*dphi_D1[ijk] + 2.0*
W_U2[ijk]*dphi_D2[ijk] + psi4*(pow(W_U0[ijk], 2)*_gamma_D0D0[ijk] +
2.0*W_U0[ijk]*W_U1[ijk]*_gamma_D0D1[ijk] + 2.0*W_U0[ijk]*W_U2[ijk]*
_gamma_D0D2[ijk] + pow(W_U1[ijk], 2)*_gamma_D1D1[ijk] + 2.0*W_U1[ijk]*
W_U2[ijk]*_gamma_D1D2[ijk] + pow(W_U2[ijk], 2)*_gamma_D2D2[ijk]) +
psim4*(_gammaI_U0U0[ijk]*pow(dphi_D0[ijk], 2) + 2.0*_gammaI_U0U1[ijk]*
dphi_D0[ijk]*dphi_D1[ijk] + 2.0*_gammaI_U0U2[ijk]*dphi_D0[ijk]*
dphi_D2[ijk] + _gammaI_U1U1[ijk]*pow(dphi_D1[ijk], 2) + 2.0*
_gammaI_U1U2[ijk]*dphi_D1[ijk]*dphi_D2[ijk] + _gammaI_U2U2[ijk]*
pow(dphi_D2[ijk], 2));

  double u_mu0 = 
sqrt(P2 + pow(enthalpy[ijk], 2))/(alpha*enthalpy[ijk]);

  u0[ijk] = u_mu0;
  }
}


/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid, CTS for conformal thin sandwich method;
// thus, Tij_IF_CTS means stress enerfy of ideal fluid in CTS method. 
// given all of the fields needed for J = -gamma(i,-mu)*T(mu,nu)*n(-nu) 
// it builds "momentum current * psi^6", where psi is conformal factor,
// and puts it to _J_U?.
// note: if patch does not contain fluid, it does nothing.
// note: it depends on u0, so better first to call Tij_IF_u0 function. */
void Tij_IF_CTS_psi6J_Ui(Patch_T *const patch)
{
  if (!IsItNSPatch(patch))
    return;
  const unsigned nn = patch->nn;
  unsigned ijk;

  /* declaring: */
  GET_FIELD(_gammaI_U0U2)
  GET_FIELD(_gammaI_U0U0)
  GET_FIELD(_gammaI_U0U1)
  GET_FIELD(_gammaI_U1U2)
  GET_FIELD(_gammaI_U1U1)
  GET_FIELD(_gammaI_U2U2)
  GET_FIELD(rho0)
  GET_FIELD(W_U1)
  GET_FIELD(W_U0)
  GET_FIELD(W_U2)
  GET_FIELD(dphi_D2)
  GET_FIELD(dphi_D1)
  GET_FIELD(dphi_D0)
  GET_FIELD(_J_U0)
  GET_FIELD(_J_U1)
  GET_FIELD(_J_U2)
  GET_FIELD(eta)
  GET_FIELD(psi)
  GET_FIELD(u0)


  for(ijk = 0; ijk < nn; ++ijk)
  {
  double alpha = 
eta[ijk]/psi[ijk];

  double psim4 = 
pow(psi[ijk], -4);

  double psi6 = 
pow(psi[ijk], 6);

  double j_u_U0 = 
alpha*psi6*rho0[ijk]*u0[ijk]*(W_U0[ijk] + psim4*(_gammaI_U0U0[ijk]*
dphi_D0[ijk] + _gammaI_U0U1[ijk]*dphi_D1[ijk] + _gammaI_U0U2[ijk]*
dphi_D2[ijk]));

  double j_u_U1 = 
alpha*psi6*rho0[ijk]*u0[ijk]*(W_U1[ijk] + psim4*(_gammaI_U0U1[ijk]*
dphi_D0[ijk] + _gammaI_U1U1[ijk]*dphi_D1[ijk] + _gammaI_U1U2[ijk]*
dphi_D2[ijk]));

  double j_u_U2 = 
alpha*psi6*rho0[ijk]*u0[ijk]*(W_U2[ijk] + psim4*(_gammaI_U0U2[ijk]*
dphi_D0[ijk] + _gammaI_U1U2[ijk]*dphi_D1[ijk] + _gammaI_U2U2[ijk]*
dphi_D2[ijk]));


  /* populating: */
  _J_U0[ijk] = j_u_U0;
  _J_U1[ijk] = j_u_U1;
  _J_U2[ijk] = j_u_U2;
  }
}


/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid, CTS for conformal thin sandwich method;
// thus, Tij_IF_CTS means stress enerfy of ideal fluid in CTS method. 
// given all of the fields needed for E = T(mu,nu)*n(-mu)*n(-nu)
// it builds "total energy density * psi^6", where psi is conformal factor,
// and puts it to _E.
// note: if patch does not contain fluid, it does nothing.
// note: it depends on u0, so better first to call Tij_IF_u0 function. */
void Tij_IF_CTS_psi6E(Patch_T *const patch)
{
  if (!IsItNSPatch(patch))
    return;
  const unsigned nn = patch->nn;
  unsigned ijk;

  /* declaring: */
  GET_FIELD(enthalpy)
  GET_FIELD(eta)
  GET_FIELD(u0)
  GET_FIELD(psi)
  GET_FIELD(rho0)
  GET_FIELD(_E)


  EoS_T *eos = initialize_EoS();
  for(ijk = 0; ijk < nn; ++ijk)
  {
    eos->h   = enthalpy[ijk];
    double p = eos->pressure(eos);
    double alpha = 
eta[ijk]/psi[ijk];

    double psi6 = 
pow(psi[ijk], 6);

    double Ebar = 
pow(alpha, 2)*enthalpy[ijk]*psi6*rho0[ijk]*pow(u0[ijk], 2) -
p;

  _E[ijk] = Ebar;
  }
  free_EoS(eos);
}

/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid, CTS for conformal thin sandwich method;
// thus, Tij_IF_CTS means stress enerfy of ideal fluid in CTS method. 
// given all of the fields needed for 
// S = T(mu,nu)*gamma(i,j)*gamma(-i,-mu)*gamma(-j,-nu)
// it builds "total energy density * psi^6", 
// where psi is conformal factor, and puts it to _S. 
// note: if patch does not contain fluid, it does nothing.
// note: it depends on u0, so better first to call Tij_IF_u0 function. */
void Tij_IF_CTS_psi6S(Patch_T *const patch)
{
  if (!IsItNSPatch(patch))
    return;
  const unsigned nn = patch->nn;
  unsigned ijk;

  /* declaring: */
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
  GET_FIELD(enthalpy)
  GET_FIELD(W_U1)
  GET_FIELD(W_U0)
  GET_FIELD(W_U2)
  GET_FIELD(dphi_D2)
  GET_FIELD(dphi_D1)
  GET_FIELD(dphi_D0)
  GET_FIELD(psi)
  GET_FIELD(rho0)
  GET_FIELD(_S)


  EoS_T *eos = initialize_EoS();
  for(ijk = 0; ijk < nn; ++ijk)
  {
    eos->h   = enthalpy[ijk];
    double p = eos->pressure(eos);
    double psim4 = 
pow(psi[ijk], -4);

    double psi4 = 
pow(psi[ijk], 4);

    double psi6 = 
pow(psi[ijk], 6);

    double P2 = 
2.0*W_U0[ijk]*dphi_D0[ijk] + 2.0*W_U1[ijk]*dphi_D1[ijk] + 2.0*
W_U2[ijk]*dphi_D2[ijk] + psi4*(pow(W_U0[ijk], 2)*_gamma_D0D0[ijk] +
2.0*W_U0[ijk]*W_U1[ijk]*_gamma_D0D1[ijk] + 2.0*W_U0[ijk]*W_U2[ijk]*
_gamma_D0D2[ijk] + pow(W_U1[ijk], 2)*_gamma_D1D1[ijk] + 2.0*W_U1[ijk]*
W_U2[ijk]*_gamma_D1D2[ijk] + pow(W_U2[ijk], 2)*_gamma_D2D2[ijk]) +
psim4*(_gammaI_U0U0[ijk]*pow(dphi_D0[ijk], 2) + 2.0*_gammaI_U0U1[ijk]*
dphi_D0[ijk]*dphi_D1[ijk] + 2.0*_gammaI_U0U2[ijk]*dphi_D0[ijk]*
dphi_D2[ijk] + _gammaI_U1U1[ijk]*pow(dphi_D1[ijk], 2) + 2.0*
_gammaI_U1U2[ijk]*dphi_D1[ijk]*dphi_D2[ijk] + _gammaI_U2U2[ijk]*
pow(dphi_D2[ijk], 2));

    double Sbar = 
P2*psi6*rho0[ijk]/enthalpy[ijk] + 3*p;

  _S[ijk] = Sbar;
  }
  free_EoS(eos);
}

/* building u0, _J^i, _E and _S */
void Tij_IF_CTS_psi6Sources(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Building _J^i, _E and _S of CTS sources ...\n");
  
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    Tij_IF_CTS_u0(patch);
    Tij_IF_CTS_psi6J_Ui(patch);
    Tij_IF_CTS_psi6E(patch);
    Tij_IF_CTS_psi6S(patch);
  }
  
  printf("Building _J^i, _E and _S of CTS sources ==> Done.\n");
  pr_clock();
  pr_line_custom('=');

}
