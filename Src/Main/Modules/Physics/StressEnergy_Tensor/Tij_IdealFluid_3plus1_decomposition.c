/*
// Alireza Rashti
// August 2019
*/

#include "Tij_IdealFluid_3plus1_decomposition.h"

/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid. thus, 
// Tij_IF means stress enerfy of ideal fluid. 
// given all of the fields needed it builds 
// the first component of fluid four velocity, i.e. 
// u_mu = (u_U0,u_U1,u_U2,u_U3). */
void Tij_IF_u0(Patch_T *const patch)
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
// IF stands for ideal fluid. thus, 
// Tij_IF means stress enerfy of ideal fluid.
// given all of the fields needed for J = -gamma(i,-mu)*T(mu,nu)*n(-nu) 
// it builds "momentum current * psi^6", where psi is conformal factor,
// and puts it to _J_U?.
// note: if patch does not contain fluid, it does nothing. */
void Tij_IF_build_psi6J_Ui(Patch_T *const patch)
{
}




/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid. thus, 
// Tij_IF means stress enerfy of ideal fluid.
// given all of the fields needed for E = T(mu,nu)*n(-mu)*n(-nu)
// it builds "total energy density * psi^6", where psi is conformal factor,
// and puts it to _E.
// note: if patch does not contain fluid, it does nothing. */
void Tij_IF_build_psi6E(Patch_T *const patch)
{
}


/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid. thus, 
// Tij_IF means stress enerfy of ideal fluid.
// given all of the fields needed for 
// S = T(mu,nu)*gamma(i,j)*gamma(-i,-mu)*gamma(-j,-nu)
// it builds "total energy density * psi^6", 
// where psi is conformal factor, and puts it to _S. 
// note: if patch does not contain fluid, it does nothing. */
void Tij_IF_build_psi6S(Patch_T *const patch)
{
}
