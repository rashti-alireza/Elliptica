/*
// Alireza Rashti
// August 2019
*/


#include "bbn_fields.h"

/* allocating all of the fields needed for construction of Initial Data */
void bbn_allocate_fields(Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    /* only if the patch covers a part of the NS add the following fields */
    if (IsItNSPatch(patch))
    {
      /* scalar for the irrotational part of fluid i.e h*u = dphi+W in NS and 
      // its partial derivatives*/
      add_field("phi",0,patch,YES);
      add_field("phi_residual",0,patch,YES);
      
      add_field("OLD_phi",0,patch,YES);
      ADD_FIELD_NoMem(dphi_D2)
      ADD_FIELD_NoMem(dphi_D1)
      ADD_FIELD_NoMem(dphi_D0)
      
      ADD_FIELD_NoMem(ddphi_D2D2)
      ADD_FIELD_NoMem(ddphi_D1D2)
      ADD_FIELD_NoMem(ddphi_D1D1)
      ADD_FIELD_NoMem(ddphi_D0D2)
      ADD_FIELD_NoMem(ddphi_D0D0)
      ADD_FIELD_NoMem(ddphi_D0D1)

      /* enthalpy in NS and its partial derivatives */
      add_field("enthalpy",0,patch,YES);
      add_field("OLD_enthalpy",0,patch,YES);
      ADD_FIELD_NoMem(denthalpy_D2)
      ADD_FIELD_NoMem(denthalpy_D1)
      ADD_FIELD_NoMem(denthalpy_D0)
  
      /* rest mass density in NS and its partial derivatives */
      add_field("rho0",0,patch,YES);
      ADD_FIELD_NoMem(drho0_D2)
      ADD_FIELD_NoMem(drho0_D1)
      ADD_FIELD_NoMem(drho0_D0)
      
      /* the first component of fluid four velocity, i.e. 
      // u_mu = (u_U0,u_U1,u_U2,u_U3) and its partial derivatives */
      add_field("u0",0,patch,YES);
      ADD_FIELD_NoMem(du0_D2)
      ADD_FIELD_NoMem(du0_D1)
      ADD_FIELD_NoMem(du0_D0)
      
      /* spin part of fluid W^i */
      add_field("W_U0",0,patch,YES);
      add_field("W_U1",0,patch,YES);
      add_field("W_U2",0,patch,YES);
      
    }
    
    /* only if patch covers horzion */
    if (IsItHorizonPatch(patch))
    {
      /* normal vector on horizon */
      ADD_FIELD(_HS_U0);
      ADD_FIELD(_HS_U1);
      ADD_FIELD(_HS_U2);
      ADD_FIELD_NoMem(_dHS_U1D0)
      ADD_FIELD_NoMem(_dHS_U1D1)
      ADD_FIELD_NoMem(_dHS_U1D2)
      ADD_FIELD_NoMem(_dHS_U0D1)
      ADD_FIELD_NoMem(_dHS_U0D0)
      ADD_FIELD_NoMem(_dHS_U0D2)
      ADD_FIELD_NoMem(_dHS_U2D2)
      ADD_FIELD_NoMem(_dHS_U2D1)
      ADD_FIELD_NoMem(_dHS_U2D0)
    }
    
    /* Hamiltonian and Momentum constraints */
    ADD_FIELD(ham_constraint_1st);
    ADD_FIELD(mom_constraint_1st_U0);
    ADD_FIELD(mom_constraint_1st_U1);
    ADD_FIELD(mom_constraint_1st_U2);
  
    ADD_FIELD(ham_constraint_2nd);
    ADD_FIELD(mom_constraint_2nd_U0);
    ADD_FIELD(mom_constraint_2nd_U1);
    ADD_FIELD(mom_constraint_2nd_U2);
    
    /* K^ij and dK^ij (extrinsic curvature) */
    ADD_FIELD(K_UiUj_U2U2)
    ADD_FIELD(K_UiUj_U1U2)
    ADD_FIELD(K_UiUj_U1U1)
    ADD_FIELD(K_UiUj_U0U2)
    ADD_FIELD(K_UiUj_U0U1)
    ADD_FIELD(K_UiUj_U0U0)
    ADD_FIELD(dK_UiUj_U0U0D2)
    ADD_FIELD(dK_UiUj_U0U0D0)
    ADD_FIELD(dK_UiUj_U0U0D1)
    ADD_FIELD(dK_UiUj_U1U1D2)
    ADD_FIELD(dK_UiUj_U1U1D0)
    ADD_FIELD(dK_UiUj_U1U1D1)
    ADD_FIELD(dK_UiUj_U0U2D0)
    ADD_FIELD(dK_UiUj_U0U2D1)
    ADD_FIELD(dK_UiUj_U0U2D2)
    ADD_FIELD(dK_UiUj_U1U2D1)
    ADD_FIELD(dK_UiUj_U1U2D0)
    ADD_FIELD(dK_UiUj_U2U2D0)
    ADD_FIELD(dK_UiUj_U1U2D2)
    ADD_FIELD(dK_UiUj_U2U2D1)
    ADD_FIELD(dK_UiUj_U2U2D2)
    ADD_FIELD(dK_UiUj_U0U1D2)
    ADD_FIELD(dK_UiUj_U0U1D1)
    ADD_FIELD(dK_UiUj_U0U1D0)

    /* conformal total energy density */
    add_field("_E",0,patch,YES);
    
    /* conformal trace of stress tensor */
    add_field("_S",0,patch,YES);
    
    /* conformal momentum current */
    add_field("_J_U0",0,patch,YES);
    add_field("_J_U1",0,patch,YES);
    add_field("_J_U2",0,patch,YES);
    
    /* conformal factor and its derivative */
    add_field("psi",0,patch,YES);
    add_field("psi_residual",0,patch,YES);
    add_field("OLD_psi",0,patch,YES);
    ADD_FIELD_NoMem(dpsi_D2)
    ADD_FIELD_NoMem(dpsi_D1)
    ADD_FIELD_NoMem(dpsi_D0)
    
    ADD_FIELD_NoMem(ddpsi_D2D2)
    ADD_FIELD_NoMem(ddpsi_D1D2)
    ADD_FIELD_NoMem(ddpsi_D1D1)
    ADD_FIELD_NoMem(ddpsi_D0D2)
    ADD_FIELD_NoMem(ddpsi_D0D0)
    ADD_FIELD_NoMem(ddpsi_D0D1)

    /* eta = lapse * psi and its partial derivatives */
    add_field("eta",0,patch,YES);
    add_field("eta_residual",0,patch,YES);
    
    add_field("OLD_eta",0,patch,YES);
    
    ADD_FIELD_NoMem(deta_D2)
    ADD_FIELD_NoMem(deta_D1)
    ADD_FIELD_NoMem(deta_D0)

    ADD_FIELD_NoMem(ddeta_D2D2)
    ADD_FIELD_NoMem(ddeta_D1D2)
    ADD_FIELD_NoMem(ddeta_D1D1)
    ADD_FIELD_NoMem(ddeta_D0D2)
    ADD_FIELD_NoMem(ddeta_D0D0)
    ADD_FIELD_NoMem(ddeta_D0D1)    
    
    /* shift, Beta^i = B0^i+B1^i, 
    // B1^i = Omega_BHNS*(-y+y_CM,x-x_CM,0)+v_r/D*(x-x_CM,y-y_CM) 
    // and its partial derivative */
    add_field("B0_U0",0,patch,YES);
    add_field("B0_U1",0,patch,YES);
    add_field("B0_U2",0,patch,YES);
    add_field("B0_U0_residual",0,patch,YES);
    add_field("B0_U1_residual",0,patch,YES);
    add_field("B0_U2_residual",0,patch,YES);
    
    add_field("B1_U0",0,patch,YES);
    add_field("B1_U1",0,patch,YES);
    add_field("B1_U2",0,patch,YES);
    add_field("Beta_U0",0,patch,YES);
    add_field("Beta_U1",0,patch,YES);
    add_field("Beta_U2",0,patch,YES);
    add_field("OLD_B0_U0",0,patch,YES);
    add_field("OLD_B0_U1",0,patch,YES);
    add_field("OLD_B0_U2",0,patch,YES);
    
    ADD_FIELD_NoMem(dB0_U1D0)
    ADD_FIELD_NoMem(dB0_U1D1)
    ADD_FIELD_NoMem(dB0_U1D2)
    ADD_FIELD_NoMem(dB0_U0D1)
    ADD_FIELD_NoMem(dB0_U0D0)
    ADD_FIELD_NoMem(dB0_U0D2)
    ADD_FIELD_NoMem(dB0_U2D2)
    ADD_FIELD_NoMem(dB0_U2D1)
    ADD_FIELD_NoMem(dB0_U2D0)

    ADD_FIELD_NoMem(ddB0_U1D2D2)
    ADD_FIELD_NoMem(ddB0_U1D0D0)
    ADD_FIELD_NoMem(ddB0_U1D0D1)
    ADD_FIELD_NoMem(ddB0_U1D0D2)
    ADD_FIELD_NoMem(ddB0_U1D1D2)
    ADD_FIELD_NoMem(ddB0_U2D0D0)
    ADD_FIELD_NoMem(ddB0_U2D1D2)
    ADD_FIELD_NoMem(ddB0_U0D2D2)
    ADD_FIELD_NoMem(ddB0_U2D1D1)
    ADD_FIELD_NoMem(ddB0_U0D0D1)
    ADD_FIELD_NoMem(ddB0_U0D1D1)
    ADD_FIELD_NoMem(ddB0_U0D1D2)
    ADD_FIELD_NoMem(ddB0_U0D0D2)
    ADD_FIELD_NoMem(ddB0_U1D1D1)
    ADD_FIELD_NoMem(ddB0_U2D0D1)
    ADD_FIELD_NoMem(ddB0_U2D2D2)
    ADD_FIELD_NoMem(ddB0_U2D0D2)
    ADD_FIELD_NoMem(ddB0_U0D0D0)

    ADD_FIELD(dBeta_U1D0)
    ADD_FIELD(dBeta_U1D1)
    ADD_FIELD(dBeta_U1D2)
    ADD_FIELD(dBeta_U0D1)
    ADD_FIELD(dBeta_U0D0)
    ADD_FIELD(dBeta_U0D2)
    ADD_FIELD(dBeta_U2D2)
    ADD_FIELD(dBeta_U2D1)
    ADD_FIELD(dBeta_U2D0)

    ADD_FIELD(ddBeta_U1D2D2)
    ADD_FIELD(ddBeta_U1D0D0)
    ADD_FIELD(ddBeta_U1D0D1)
    ADD_FIELD(ddBeta_U1D0D2)
    ADD_FIELD(ddBeta_U1D1D2)
    ADD_FIELD(ddBeta_U2D0D0)
    ADD_FIELD(ddBeta_U2D1D2)
    ADD_FIELD(ddBeta_U0D2D2)
    ADD_FIELD(ddBeta_U2D1D1)
    ADD_FIELD(ddBeta_U0D0D1)
    ADD_FIELD(ddBeta_U0D1D1)
    ADD_FIELD(ddBeta_U0D1D2)
    ADD_FIELD(ddBeta_U0D0D2)
    ADD_FIELD(ddBeta_U1D1D1)
    ADD_FIELD(ddBeta_U2D0D1)
    ADD_FIELD(ddBeta_U2D2D2)
    ADD_FIELD(ddBeta_U2D0D2)
    ADD_FIELD(ddBeta_U0D0D0)
  
    ADD_FIELD(dB1_U1D0)
    ADD_FIELD(dB1_U1D1)
    ADD_FIELD(dB1_U1D2)
    ADD_FIELD(dB1_U0D1)
    ADD_FIELD(dB1_U0D0)
    ADD_FIELD(dB1_U0D2)
    ADD_FIELD(dB1_U2D2)
    ADD_FIELD(dB1_U2D1)
    ADD_FIELD(dB1_U2D0)

    ADD_FIELD(ddB1_U1D2D2)
    ADD_FIELD(ddB1_U1D0D0)
    ADD_FIELD(ddB1_U1D0D1)
    ADD_FIELD(ddB1_U1D0D2)
    ADD_FIELD(ddB1_U1D1D2)
    ADD_FIELD(ddB1_U2D0D0)
    ADD_FIELD(ddB1_U2D1D2)
    ADD_FIELD(ddB1_U0D2D2)
    ADD_FIELD(ddB1_U2D1D1)
    ADD_FIELD(ddB1_U0D0D1)
    ADD_FIELD(ddB1_U0D1D1)
    ADD_FIELD(ddB1_U0D1D2)
    ADD_FIELD(ddB1_U0D0D2)
    ADD_FIELD(ddB1_U1D1D1)
    ADD_FIELD(ddB1_U2D0D1)
    ADD_FIELD(ddB1_U2D2D2)
    ADD_FIELD(ddB1_U2D0D2)
    ADD_FIELD(ddB1_U0D0D0)
  
    /* conformal metric: _gamma_DiDj */
    ADD_FIELD(_gamma_D2D2)
    ADD_FIELD(_gamma_D0D2)
    ADD_FIELD(_gamma_D0D0)
    ADD_FIELD(_gamma_D0D1)
    ADD_FIELD(_gamma_D1D2)
    ADD_FIELD(_gamma_D1D1)

    /* conformal metric inverse _gammaI_UiUj I stands for inverse */
    ADD_FIELD(_gammaI_U0U2)
    ADD_FIELD(_gammaI_U0U0)
    ADD_FIELD(_gammaI_U0U1)
    ADD_FIELD(_gammaI_U1U2)
    ADD_FIELD(_gammaI_U1U1)
    ADD_FIELD(_gammaI_U2U2)

    /* Christoffer symbols made up of conformal metric */
    ADD_FIELD(_Gamma_U2D1D1)
    ADD_FIELD(_Gamma_U2D1D2)
    ADD_FIELD(_Gamma_U0D1D1)
    ADD_FIELD(_Gamma_U2D0D2)
    ADD_FIELD(_Gamma_U2D2D2)
    ADD_FIELD(_Gamma_U0D1D2)
    ADD_FIELD(_Gamma_U0D0D2)
    ADD_FIELD(_Gamma_U0D0D1)
    ADD_FIELD(_Gamma_U0D0D0)
    ADD_FIELD(_Gamma_U1D2D2)
    ADD_FIELD(_Gamma_U2D0D1)
    ADD_FIELD(_Gamma_U0D2D2)
    ADD_FIELD(_Gamma_U2D0D0)
    ADD_FIELD(_Gamma_U1D0D2)
    ADD_FIELD(_Gamma_U1D1D2)
    ADD_FIELD(_Gamma_U1D0D0)
    ADD_FIELD(_Gamma_U1D0D1)
    ADD_FIELD(_Gamma_U1D1D1)
    
    /* partial derivative of Christoffer symbols */
    ADD_FIELD_NoMem(_dGamma_U2D2D2D2)
    ADD_FIELD_NoMem(_dGamma_U2D2D2D0)
    ADD_FIELD_NoMem(_dGamma_U2D2D2D1)
    ADD_FIELD_NoMem(_dGamma_U2D0D0D2)
    ADD_FIELD_NoMem(_dGamma_U1D1D2D2)
    ADD_FIELD_NoMem(_dGamma_U2D0D0D0)
    ADD_FIELD_NoMem(_dGamma_U1D1D2D0)
    ADD_FIELD_NoMem(_dGamma_U1D1D2D1)
    ADD_FIELD_NoMem(_dGamma_U2D1D1D0)
    ADD_FIELD_NoMem(_dGamma_U2D0D0D1)
    ADD_FIELD_NoMem(_dGamma_U2D0D2D1)
    ADD_FIELD_NoMem(_dGamma_U1D0D1D0)
    ADD_FIELD_NoMem(_dGamma_U1D0D1D1)
    ADD_FIELD_NoMem(_dGamma_U1D0D1D2)
    ADD_FIELD_NoMem(_dGamma_U1D2D2D1)
    ADD_FIELD_NoMem(_dGamma_U1D0D0D1)
    ADD_FIELD_NoMem(_dGamma_U1D0D0D0)
    ADD_FIELD_NoMem(_dGamma_U1D0D0D2)
    ADD_FIELD_NoMem(_dGamma_U0D1D2D2)
    ADD_FIELD_NoMem(_dGamma_U0D1D2D1)
    ADD_FIELD_NoMem(_dGamma_U0D1D2D0)
    ADD_FIELD_NoMem(_dGamma_U2D0D2D0)
    ADD_FIELD_NoMem(_dGamma_U1D0D2D2)
    ADD_FIELD_NoMem(_dGamma_U1D0D2D1)
    ADD_FIELD_NoMem(_dGamma_U1D0D2D0)
    ADD_FIELD_NoMem(_dGamma_U2D1D1D2)
    ADD_FIELD_NoMem(_dGamma_U2D0D2D2)
    ADD_FIELD_NoMem(_dGamma_U0D0D1D0)
    ADD_FIELD_NoMem(_dGamma_U1D2D2D0)
    ADD_FIELD_NoMem(_dGamma_U2D1D2D1)
    ADD_FIELD_NoMem(_dGamma_U2D0D1D2)
    ADD_FIELD_NoMem(_dGamma_U2D0D1D1)
    ADD_FIELD_NoMem(_dGamma_U2D0D1D0)
    ADD_FIELD_NoMem(_dGamma_U2D1D2D2)
    ADD_FIELD_NoMem(_dGamma_U0D1D1D0)
    ADD_FIELD_NoMem(_dGamma_U0D1D1D1)
    ADD_FIELD_NoMem(_dGamma_U0D1D1D2)
    ADD_FIELD_NoMem(_dGamma_U1D2D2D2)
    ADD_FIELD_NoMem(_dGamma_U1D1D1D1)
    ADD_FIELD_NoMem(_dGamma_U1D1D1D0)
    ADD_FIELD_NoMem(_dGamma_U1D1D1D2)
    ADD_FIELD_NoMem(_dGamma_U0D0D1D1)
    ADD_FIELD_NoMem(_dGamma_U0D0D2D2)
    ADD_FIELD_NoMem(_dGamma_U0D0D2D0)
    ADD_FIELD_NoMem(_dGamma_U0D0D2D1)
    ADD_FIELD_NoMem(_dGamma_U2D1D2D0)
    ADD_FIELD_NoMem(_dGamma_U0D0D0D0)
    ADD_FIELD_NoMem(_dGamma_U0D0D0D1)
    ADD_FIELD_NoMem(_dGamma_U0D0D0D2)
    ADD_FIELD_NoMem(_dGamma_U2D1D1D1)
    ADD_FIELD_NoMem(_dGamma_U0D2D2D0)
    ADD_FIELD_NoMem(_dGamma_U0D2D2D1)
    ADD_FIELD_NoMem(_dGamma_U0D2D2D2)
    ADD_FIELD_NoMem(_dGamma_U0D0D1D2)

    /* metric: gamma_DiDj */
    ADD_FIELD(gamma_D2D2)
    ADD_FIELD(gamma_D0D2)
    ADD_FIELD(gamma_D0D0)
    ADD_FIELD(gamma_D0D1)
    ADD_FIELD(gamma_D1D2)
    ADD_FIELD(gamma_D1D1)

    /* inverse gammaI_UiUj I stands for inverse */
    ADD_FIELD(gammaI_U0U2)
    ADD_FIELD(gammaI_U0U0)
    ADD_FIELD(gammaI_U0U1)
    ADD_FIELD(gammaI_U1U2)
    ADD_FIELD(gammaI_U1U1)
    ADD_FIELD(gammaI_U2U2)
    
    /* gamma derivatives */
    ADD_FIELD_NoMem(dgamma_D0D0D1)
    ADD_FIELD_NoMem(dgamma_D0D0D0)
    ADD_FIELD_NoMem(dgamma_D2D2D2)
    ADD_FIELD_NoMem(dgamma_D0D0D2)
    ADD_FIELD_NoMem(dgamma_D0D2D1)
    ADD_FIELD_NoMem(dgamma_D1D1D0)
    ADD_FIELD_NoMem(dgamma_D1D1D2)
    ADD_FIELD_NoMem(dgamma_D2D2D0)
    ADD_FIELD_NoMem(dgamma_D2D2D1)
    ADD_FIELD_NoMem(dgamma_D0D1D0)
    ADD_FIELD_NoMem(dgamma_D0D1D1)
    ADD_FIELD_NoMem(dgamma_D0D1D2)
    ADD_FIELD_NoMem(dgamma_D0D2D0)
    ADD_FIELD_NoMem(dgamma_D1D1D1)
    ADD_FIELD_NoMem(dgamma_D1D2D1)
    ADD_FIELD_NoMem(dgamma_D1D2D2)
    ADD_FIELD_NoMem(dgamma_D1D2D0)
    ADD_FIELD_NoMem(dgamma_D0D2D2)
    
    /* Christoffer symbols made up of metric */
    ADD_FIELD(Gamma_U2D1D1)
    ADD_FIELD(Gamma_U2D1D2)
    ADD_FIELD(Gamma_U0D1D1)
    ADD_FIELD(Gamma_U2D0D2)
    ADD_FIELD(Gamma_U2D2D2)
    ADD_FIELD(Gamma_U0D1D2)
    ADD_FIELD(Gamma_U0D0D2)
    ADD_FIELD(Gamma_U0D0D1)
    ADD_FIELD(Gamma_U0D0D0)
    ADD_FIELD(Gamma_U1D2D2)
    ADD_FIELD(Gamma_U2D0D1)
    ADD_FIELD(Gamma_U0D2D2)
    ADD_FIELD(Gamma_U2D0D0)
    ADD_FIELD(Gamma_U1D0D2)
    ADD_FIELD(Gamma_U1D1D2)
    ADD_FIELD(Gamma_U1D0D0)
    ADD_FIELD(Gamma_U1D0D1)
    ADD_FIELD(Gamma_U1D1D1)
    
    /* partial derivative of Christoffer symbols */
    ADD_FIELD_NoMem(dGamma_U2D2D2D2)
    ADD_FIELD_NoMem(dGamma_U2D2D2D0)
    ADD_FIELD_NoMem(dGamma_U2D2D2D1)
    ADD_FIELD_NoMem(dGamma_U2D0D0D2)
    ADD_FIELD_NoMem(dGamma_U1D1D2D2)
    ADD_FIELD_NoMem(dGamma_U2D0D0D0)
    ADD_FIELD_NoMem(dGamma_U1D1D2D0)
    ADD_FIELD_NoMem(dGamma_U1D1D2D1)
    ADD_FIELD_NoMem(dGamma_U2D1D1D0)
    ADD_FIELD_NoMem(dGamma_U2D0D0D1)
    ADD_FIELD_NoMem(dGamma_U2D0D2D1)
    ADD_FIELD_NoMem(dGamma_U1D0D1D0)
    ADD_FIELD_NoMem(dGamma_U1D0D1D1)
    ADD_FIELD_NoMem(dGamma_U1D0D1D2)
    ADD_FIELD_NoMem(dGamma_U1D2D2D1)
    ADD_FIELD_NoMem(dGamma_U1D0D0D1)
    ADD_FIELD_NoMem(dGamma_U1D0D0D0)
    ADD_FIELD_NoMem(dGamma_U1D0D0D2)
    ADD_FIELD_NoMem(dGamma_U0D1D2D2)
    ADD_FIELD_NoMem(dGamma_U0D1D2D1)
    ADD_FIELD_NoMem(dGamma_U0D1D2D0)
    ADD_FIELD_NoMem(dGamma_U2D0D2D0)
    ADD_FIELD_NoMem(dGamma_U1D0D2D2)
    ADD_FIELD_NoMem(dGamma_U1D0D2D1)
    ADD_FIELD_NoMem(dGamma_U1D0D2D0)
    ADD_FIELD_NoMem(dGamma_U2D1D1D2)
    ADD_FIELD_NoMem(dGamma_U2D0D2D2)
    ADD_FIELD_NoMem(dGamma_U0D0D1D0)
    ADD_FIELD_NoMem(dGamma_U1D2D2D0)
    ADD_FIELD_NoMem(dGamma_U2D1D2D1)
    ADD_FIELD_NoMem(dGamma_U2D0D1D2)
    ADD_FIELD_NoMem(dGamma_U2D0D1D1)
    ADD_FIELD_NoMem(dGamma_U2D0D1D0)
    ADD_FIELD_NoMem(dGamma_U2D1D2D2)
    ADD_FIELD_NoMem(dGamma_U0D1D1D0)
    ADD_FIELD_NoMem(dGamma_U0D1D1D1)
    ADD_FIELD_NoMem(dGamma_U0D1D1D2)
    ADD_FIELD_NoMem(dGamma_U1D2D2D2)
    ADD_FIELD_NoMem(dGamma_U1D1D1D1)
    ADD_FIELD_NoMem(dGamma_U1D1D1D0)
    ADD_FIELD_NoMem(dGamma_U1D1D1D2)
    ADD_FIELD_NoMem(dGamma_U0D0D1D1)
    ADD_FIELD_NoMem(dGamma_U0D0D2D2)
    ADD_FIELD_NoMem(dGamma_U0D0D2D0)
    ADD_FIELD_NoMem(dGamma_U0D0D2D1)
    ADD_FIELD_NoMem(dGamma_U2D1D2D0)
    ADD_FIELD_NoMem(dGamma_U0D0D0D0)
    ADD_FIELD_NoMem(dGamma_U0D0D0D1)
    ADD_FIELD_NoMem(dGamma_U0D0D0D2)
    ADD_FIELD_NoMem(dGamma_U2D1D1D1)
    ADD_FIELD_NoMem(dGamma_U0D2D2D0)
    ADD_FIELD_NoMem(dGamma_U0D2D2D1)
    ADD_FIELD_NoMem(dGamma_U0D2D2D2)
    ADD_FIELD_NoMem(dGamma_U0D0D1D2)

    /* Ricci made up of conformal metric _gamma */
    add_field("_R",0,patch,YES);
    ADD_FIELD(_Ric_D1D2)
    ADD_FIELD(_Ric_D2D2)
    ADD_FIELD(_Ric_D0D2)
    ADD_FIELD(_Ric_D0D1)
    ADD_FIELD(_Ric_D0D0)
    ADD_FIELD(_Ric_D1D1)

    /* extrinsic curvature */
    add_field("K",0,patch,YES);
    ADD_FIELD_NoMem(dK_D2)
    ADD_FIELD_NoMem(dK_D1)
    ADD_FIELD_NoMem(dK_D0)

    /*_A_UiUj psi^10*A^{ij}*/
    ADD_FIELD(_A_UiUj_U2U2)
    ADD_FIELD(_A_UiUj_U1U2)
    ADD_FIELD(_A_UiUj_U1U1)
    ADD_FIELD(_A_UiUj_U0U2)
    ADD_FIELD(_A_UiUj_U0U1)
    ADD_FIELD(_A_UiUj_U0U0)
    
    /* contraction _A_{ij}*_A^{ij} */
    ADD_FIELD(_Aij2)
    
    /* derivatives of _A^{ij} */
    ADD_FIELD(_dA_UiUj_U2U2D2)
    ADD_FIELD(_dA_UiUj_U2U2D0)
    ADD_FIELD(_dA_UiUj_U2U2D1)
    ADD_FIELD(_dA_UiUj_U1U1D2)
    ADD_FIELD(_dA_UiUj_U1U1D0)
    ADD_FIELD(_dA_UiUj_U1U1D1)
    ADD_FIELD(_dA_UiUj_U0U0D2)
    ADD_FIELD(_dA_UiUj_U0U0D0)
    ADD_FIELD(_dA_UiUj_U0U0D1)
    ADD_FIELD(_dA_UiUj_U0U1D2)
    ADD_FIELD(_dA_UiUj_U0U1D1)
    ADD_FIELD(_dA_UiUj_U0U1D0)
    ADD_FIELD(_dA_UiUj_U1U2D1)
    ADD_FIELD(_dA_UiUj_U1U2D0)
    ADD_FIELD(_dA_UiUj_U1U2D2)
    ADD_FIELD(_dA_UiUj_U0U2D0)
    ADD_FIELD(_dA_UiUj_U0U2D1)
    ADD_FIELD(_dA_UiUj_U0U2D2)
    
  }
}

/* taking partial derivatives of fields needed in equations */
void bbn_partial_derivatives_fields(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Taking partial derivatives of fields ...\n");
  
  const unsigned np = grid->np;
  unsigned p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    /* phi derivatives */
    bbn_update_derivative_phi(patch);
    
    /* enthalpy derivatives */
    bbn_update_derivative_enthalpy(patch);
    
    /* rho0 derivatives */
    bbn_update_derivative_rho0(patch);
    
    /* u0 derivatives */
    bbn_update_derivative_u0(patch);
    
    /* normal vector on Horizon derivatives */
    bbn_update_derivative_HS(patch);
    
    /* K derivatives */
    bbn_update_derivative_K(patch);
    
    /* psi derivatives */
    bbn_update_derivative_psi(patch);
   
    /* eta derivatives */
    bbn_update_derivative_eta(patch);
    
    /* B0^i*/
    bbn_update_derivative_B0_U0(patch);
    bbn_update_derivative_B0_U1(patch);
    bbn_update_derivative_B0_U2(patch);
    
    /* B1^i*/
    bbn_update_derivative_B1_U0(patch);
    bbn_update_derivative_B1_U1(patch);
    bbn_update_derivative_B1_U2(patch);
        
    /* Beta^i*/
    bbn_update_derivative_Beta_U0(patch);
    bbn_update_derivative_Beta_U1(patch);
    bbn_update_derivative_Beta_U2(patch);

  }
  
  printf("Taking partial derivatives of fields ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* updating Beta_U0 */
void bbn_update_Beta_U0(Patch_T *const patch)
{
  PREP_FIELD(Beta_U0)
  GET_FIELD(B0_U0)
  GET_FIELD(B1_U0)
  
  unsigned ijk;
  FOR_ALL_POINTS(ijk,patch)
  {
    Beta_U0[ijk] = B0_U0[ijk]+B1_U0[ijk];
  }
  
}

/* updating Beta_U1 */
void bbn_update_Beta_U1(Patch_T *const patch)
{
  PREP_FIELD(Beta_U1)
  GET_FIELD(B0_U1)
  GET_FIELD(B1_U1)
  
  unsigned ijk;
  FOR_ALL_POINTS(ijk,patch)
  {
    Beta_U1[ijk] = B0_U1[ijk]+B1_U1[ijk];
  }
  
}

/* updating Beta_U2 */
void bbn_update_Beta_U2(Patch_T *const patch)
{
  PREP_FIELD(Beta_U2)
  GET_FIELD(B0_U2)
  GET_FIELD(B1_U2)
  
  unsigned ijk;
  FOR_ALL_POINTS(ijk,patch)
  {
    Beta_U2[ijk] = B0_U2[ijk]+B1_U2[ijk];
  }
  
}

/* updating derivative */
void bbn_update_derivative_phi(Patch_T *const patch)
{
  /* only if the patch covers a part of the NS add the following fields */
  if (IsItNSPatch(patch))
  {
    /* phi's derivatives */
    DECLARE_FIELD(phi)

    DECLARE_AND_EMPTY_FIELD(dphi_D2)
    DECLARE_AND_EMPTY_FIELD(dphi_D1)
    DECLARE_AND_EMPTY_FIELD(dphi_D0)
    
    DECLARE_AND_EMPTY_FIELD(ddphi_D2D2)
    DECLARE_AND_EMPTY_FIELD(ddphi_D1D2)
    DECLARE_AND_EMPTY_FIELD(ddphi_D1D1)
    DECLARE_AND_EMPTY_FIELD(ddphi_D0D2)
    DECLARE_AND_EMPTY_FIELD(ddphi_D0D0)
    DECLARE_AND_EMPTY_FIELD(ddphi_D0D1)

    dphi_D2->v = Partial_Derivative(phi,"z");
    dphi_D1->v = Partial_Derivative(phi,"y");
    dphi_D0->v = Partial_Derivative(phi,"x");
    
    ddphi_D2D2->v = Partial_Derivative(dphi_D2,"z");
    ddphi_D1D2->v = Partial_Derivative(dphi_D1,"z");
    ddphi_D1D1->v = Partial_Derivative(dphi_D1,"y");
    ddphi_D0D2->v = Partial_Derivative(dphi_D0,"z");
    ddphi_D0D0->v = Partial_Derivative(dphi_D0,"x");
    ddphi_D0D1->v = Partial_Derivative(dphi_D0,"y");
  }
}


/* updating derivative */
void bbn_update_derivative_enthalpy(Patch_T *const patch)
{
  /* only if the patch covers a part of the NS add the following fields */
  if (IsItNSPatch(patch))
  {
    /* enthalpy derivatives */
    DECLARE_FIELD(enthalpy)

    DECLARE_AND_EMPTY_FIELD(denthalpy_D2)
    DECLARE_AND_EMPTY_FIELD(denthalpy_D1)
    DECLARE_AND_EMPTY_FIELD(denthalpy_D0)
    
    denthalpy_D2->v = Partial_Derivative(enthalpy,"z");
    denthalpy_D1->v = Partial_Derivative(enthalpy,"y");
    denthalpy_D0->v = Partial_Derivative(enthalpy,"x");
  }
}

/* update rho0 */
void bbn_update_rho0(Patch_T *const patch)
{
  /* only if the patch covers a part of the NS add the following fields */
  if (IsItNSPatch(patch))
  {
    EoS_T *eos = initialize_EoS();
    GET_FIELD(enthalpy)
    PREP_FIELD(rho0)
    const unsigned nn = patch->nn;
    unsigned ijk;

    for (ijk = 0; ijk < nn; ++ijk)
    {
      eos->h    = enthalpy[ijk];
      rho0[ijk] = eos->rest_mass_density(eos);
      
      /* make sure rho0 won't get nan due to Gibbs phenomena or 
      // when finding NS surface. */
      if (!isfinite(rho0[ijk]))
      {
        //printf("Put rho0(h = %g) = 0.\n",enthalpy[ijk]);
        rho0[ijk] = 0;
      }
    }
    free_EoS(eos);
  }
}

/* updating derivative */
void bbn_update_derivative_rho0(Patch_T *const patch)
{
  /* only if the patch covers a part of the NS add the following fields */
  if (IsItNSPatch(patch))
  {
    /* rho0 derivatives */
    DECLARE_FIELD(rho0)

    DECLARE_AND_EMPTY_FIELD(drho0_D2)
    DECLARE_AND_EMPTY_FIELD(drho0_D1)
    DECLARE_AND_EMPTY_FIELD(drho0_D0)
    
    drho0_D2->v = Partial_Derivative(rho0,"z");
    drho0_D1->v = Partial_Derivative(rho0,"y");
    drho0_D0->v = Partial_Derivative(rho0,"x");
  }    

}

/* updating derivative */
void bbn_update_derivative_u0(Patch_T *const patch)
{
  /* only if the patch covers a part of the NS add the following fields */
  if (IsItNSPatch(patch))
  {
    /* u0 derivatives */
    DECLARE_FIELD(u0)

    DECLARE_AND_EMPTY_FIELD(du0_D2)
    DECLARE_AND_EMPTY_FIELD(du0_D1)
    DECLARE_AND_EMPTY_FIELD(du0_D0)
    
    du0_D2->v = Partial_Derivative(u0,"z");
    du0_D1->v = Partial_Derivative(u0,"y");
    du0_D0->v = Partial_Derivative(u0,"x");
  }
  
}

/* updating derivative */
void bbn_update_derivative_HS(Patch_T *const patch)
{
  /* only if patch covers horzion */
  if (IsItHorizonPatch(patch))
  {
    /* normal vector on horizon */
    DECLARE_FIELD(_HS_U0);
    DECLARE_FIELD(_HS_U1);
    DECLARE_FIELD(_HS_U2);
    
    DECLARE_AND_EMPTY_FIELD(_dHS_U1D0)
    DECLARE_AND_EMPTY_FIELD(_dHS_U1D1)
    DECLARE_AND_EMPTY_FIELD(_dHS_U1D2)
    DECLARE_AND_EMPTY_FIELD(_dHS_U0D1)
    DECLARE_AND_EMPTY_FIELD(_dHS_U0D0)
    DECLARE_AND_EMPTY_FIELD(_dHS_U0D2)
    DECLARE_AND_EMPTY_FIELD(_dHS_U2D2)
    DECLARE_AND_EMPTY_FIELD(_dHS_U2D1)
    DECLARE_AND_EMPTY_FIELD(_dHS_U2D0)
    
    _dHS_U1D0->v = Partial_Derivative(_HS_U1,"x");
    _dHS_U1D1->v = Partial_Derivative(_HS_U1,"y");
    _dHS_U1D2->v = Partial_Derivative(_HS_U1,"z");
    _dHS_U0D1->v = Partial_Derivative(_HS_U0,"y");
    _dHS_U0D0->v = Partial_Derivative(_HS_U0,"x");
    _dHS_U0D2->v = Partial_Derivative(_HS_U0,"z");
    _dHS_U2D2->v = Partial_Derivative(_HS_U2,"z");
    _dHS_U2D1->v = Partial_Derivative(_HS_U2,"y");
    _dHS_U2D0->v = Partial_Derivative(_HS_U2,"x");

  }
  
}

/* updating derivative */
void bbn_update_derivative_K(Patch_T *const patch)
{
  DECLARE_FIELD(K)

  DECLARE_AND_EMPTY_FIELD(dK_D2)
  DECLARE_AND_EMPTY_FIELD(dK_D1)
  DECLARE_AND_EMPTY_FIELD(dK_D0)

  dK_D2->v = Partial_Derivative(K,"z");
  dK_D1->v = Partial_Derivative(K,"y");
  dK_D0->v = Partial_Derivative(K,"x");
}

/* updating derivative */
void bbn_update_derivative_psi(Patch_T *const patch)
{
  DECLARE_FIELD(psi)

  DECLARE_AND_EMPTY_FIELD(dpsi_D2)
  DECLARE_AND_EMPTY_FIELD(dpsi_D1)
  DECLARE_AND_EMPTY_FIELD(dpsi_D0)
  
  DECLARE_AND_EMPTY_FIELD(ddpsi_D2D2)
  DECLARE_AND_EMPTY_FIELD(ddpsi_D1D2)
  DECLARE_AND_EMPTY_FIELD(ddpsi_D1D1)
  DECLARE_AND_EMPTY_FIELD(ddpsi_D0D2)
  DECLARE_AND_EMPTY_FIELD(ddpsi_D0D0)
  DECLARE_AND_EMPTY_FIELD(ddpsi_D0D1)

  dpsi_D2->v = Partial_Derivative(psi,"z");
  dpsi_D1->v = Partial_Derivative(psi,"y");
  dpsi_D0->v = Partial_Derivative(psi,"x");
  
  ddpsi_D2D2->v = Partial_Derivative(dpsi_D2,"z");
  ddpsi_D1D2->v = Partial_Derivative(dpsi_D1,"z");
  ddpsi_D1D1->v = Partial_Derivative(dpsi_D1,"y");
  ddpsi_D0D2->v = Partial_Derivative(dpsi_D0,"z");
  ddpsi_D0D0->v = Partial_Derivative(dpsi_D0,"x");
  ddpsi_D0D1->v = Partial_Derivative(dpsi_D0,"y");
}

/* updating derivative */
void bbn_update_derivative_eta(Patch_T *const patch)
{
  DECLARE_FIELD(eta)

  DECLARE_AND_EMPTY_FIELD(deta_D2)
  DECLARE_AND_EMPTY_FIELD(deta_D1)
  DECLARE_AND_EMPTY_FIELD(deta_D0)
  
  DECLARE_AND_EMPTY_FIELD(ddeta_D2D2)
  DECLARE_AND_EMPTY_FIELD(ddeta_D1D2)
  DECLARE_AND_EMPTY_FIELD(ddeta_D1D1)
  DECLARE_AND_EMPTY_FIELD(ddeta_D0D2)
  DECLARE_AND_EMPTY_FIELD(ddeta_D0D0)
  DECLARE_AND_EMPTY_FIELD(ddeta_D0D1)    

  deta_D2->v = Partial_Derivative(eta,"z");
  deta_D1->v = Partial_Derivative(eta,"y");
  deta_D0->v = Partial_Derivative(eta,"x");
  
  ddeta_D2D2->v = Partial_Derivative(deta_D2,"z");
  ddeta_D1D2->v = Partial_Derivative(deta_D1,"z");
  ddeta_D1D1->v = Partial_Derivative(deta_D1,"y");
  ddeta_D0D2->v = Partial_Derivative(deta_D0,"z");
  ddeta_D0D0->v = Partial_Derivative(deta_D0,"x");
  ddeta_D0D1->v = Partial_Derivative(deta_D0,"y");

}

/* updating derivative */
void bbn_update_derivative_B0_U0(Patch_T *const patch)
{
  DECLARE_FIELD(B0_U0)
  
  DECLARE_AND_EMPTY_FIELD(dB0_U0D1)
  DECLARE_AND_EMPTY_FIELD(dB0_U0D0)
  DECLARE_AND_EMPTY_FIELD(dB0_U0D2)
  DECLARE_AND_EMPTY_FIELD(ddB0_U0D2D2)
  DECLARE_AND_EMPTY_FIELD(ddB0_U0D0D1)
  DECLARE_AND_EMPTY_FIELD(ddB0_U0D1D1)
  DECLARE_AND_EMPTY_FIELD(ddB0_U0D1D2)
  DECLARE_AND_EMPTY_FIELD(ddB0_U0D0D2)
  DECLARE_AND_EMPTY_FIELD(ddB0_U0D0D0)

  dB0_U0D1->v = Partial_Derivative(B0_U0,"y");
  dB0_U0D0->v = Partial_Derivative(B0_U0,"x");
  dB0_U0D2->v = Partial_Derivative(B0_U0,"z");
  ddB0_U0D2D2->v = Partial_Derivative(dB0_U0D2,"z");
  ddB0_U0D0D1->v = Partial_Derivative(dB0_U0D0,"y");
  ddB0_U0D1D1->v = Partial_Derivative(dB0_U0D1,"y");
  ddB0_U0D1D2->v = Partial_Derivative(dB0_U0D1,"z");
  ddB0_U0D0D2->v = Partial_Derivative(dB0_U0D0,"z");
  ddB0_U0D0D0->v = Partial_Derivative(dB0_U0D0,"x");

}

/* updating derivative */
void bbn_update_derivative_B0_U1(Patch_T *const patch)
{
  
  DECLARE_FIELD(B0_U1)
  
  DECLARE_AND_EMPTY_FIELD(dB0_U1D0)
  DECLARE_AND_EMPTY_FIELD(dB0_U1D1)
  DECLARE_AND_EMPTY_FIELD(dB0_U1D2)
 
  DECLARE_AND_EMPTY_FIELD(ddB0_U1D2D2)
  DECLARE_AND_EMPTY_FIELD(ddB0_U1D0D0)
  DECLARE_AND_EMPTY_FIELD(ddB0_U1D0D1)
  DECLARE_AND_EMPTY_FIELD(ddB0_U1D0D2)
  DECLARE_AND_EMPTY_FIELD(ddB0_U1D1D2)
  DECLARE_AND_EMPTY_FIELD(ddB0_U1D1D1)

  dB0_U1D0->v = Partial_Derivative(B0_U1,"x");
  dB0_U1D1->v = Partial_Derivative(B0_U1,"y");
  dB0_U1D2->v = Partial_Derivative(B0_U1,"z");

  ddB0_U1D2D2->v = Partial_Derivative(dB0_U1D2,"z");
  ddB0_U1D0D0->v = Partial_Derivative(dB0_U1D0,"x");
  ddB0_U1D0D1->v = Partial_Derivative(dB0_U1D0,"y");
  ddB0_U1D0D2->v = Partial_Derivative(dB0_U1D0,"z");
  ddB0_U1D1D2->v = Partial_Derivative(dB0_U1D1,"z");
  ddB0_U1D1D1->v = Partial_Derivative(dB0_U1D1,"y");

}

/* updating derivative */
void bbn_update_derivative_B0_U2(Patch_T *const patch)
{
  DECLARE_FIELD(B0_U2)
  
  DECLARE_AND_EMPTY_FIELD(dB0_U2D2)
  DECLARE_AND_EMPTY_FIELD(dB0_U2D1)
  DECLARE_AND_EMPTY_FIELD(dB0_U2D0)
 
  DECLARE_AND_EMPTY_FIELD(ddB0_U2D0D0)
  DECLARE_AND_EMPTY_FIELD(ddB0_U2D1D2)
  DECLARE_AND_EMPTY_FIELD(ddB0_U2D1D1)
  DECLARE_AND_EMPTY_FIELD(ddB0_U2D0D1)
  DECLARE_AND_EMPTY_FIELD(ddB0_U2D2D2)
  DECLARE_AND_EMPTY_FIELD(ddB0_U2D0D2)

  dB0_U2D2->v = Partial_Derivative(B0_U2,"z");
  dB0_U2D1->v = Partial_Derivative(B0_U2,"y");
  dB0_U2D0->v = Partial_Derivative(B0_U2,"x");

  ddB0_U2D0D0->v = Partial_Derivative(dB0_U2D0,"x");
  ddB0_U2D1D2->v = Partial_Derivative(dB0_U2D1,"z");
  ddB0_U2D1D1->v = Partial_Derivative(dB0_U2D1,"y");
  ddB0_U2D0D1->v = Partial_Derivative(dB0_U2D0,"y");
  ddB0_U2D2D2->v = Partial_Derivative(dB0_U2D2,"z");
  ddB0_U2D0D2->v = Partial_Derivative(dB0_U2D0,"z");
  
}

/* updating derivative */
void bbn_update_derivative_B1_U2(Patch_T *const patch)
{
  //const double D          = GetParameterD_E("BH_NS_separation");
  //const double Vr         = GetParameterD_E("BH_NS_infall_velocity");
  //const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity"); 
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  PREP_FIELD(dB1_U2D2)
  PREP_FIELD(dB1_U2D1)
  PREP_FIELD(dB1_U2D0)
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    dB1_U2D2[ijk] = 0;
    dB1_U2D1[ijk] = 0;
    dB1_U2D0[ijk] = 0;
  }
  
  DECLARE_AND_EMPTY_FIELD(ddB1_U2D0D0)
  DECLARE_AND_EMPTY_FIELD(ddB1_U2D1D2)
  DECLARE_AND_EMPTY_FIELD(ddB1_U2D1D1)
  DECLARE_AND_EMPTY_FIELD(ddB1_U2D0D1)
  DECLARE_AND_EMPTY_FIELD(ddB1_U2D2D2)
  DECLARE_AND_EMPTY_FIELD(ddB1_U2D0D2)
  
  /* since the second derivative is zero */
  ddB1_U2D0D0->v = alloc_double(nn);
  ddB1_U2D1D2->v = alloc_double(nn);
  ddB1_U2D1D1->v = alloc_double(nn);
  ddB1_U2D0D1->v = alloc_double(nn);
  ddB1_U2D2D2->v = alloc_double(nn);
  ddB1_U2D0D2->v = alloc_double(nn);
}

/* updating derivative */
void bbn_update_derivative_B1_U1(Patch_T *const patch)
{
  const double D          = GetParameterD_E("BH_NS_separation");
  const double Vr         = GetParameterD_E("BH_NS_infall_velocity");
  const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity"); 
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  PREP_FIELD(dB1_U1D2)
  PREP_FIELD(dB1_U1D1)
  PREP_FIELD(dB1_U1D0)
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    dB1_U1D2[ijk] = 0;
    dB1_U1D1[ijk] = Vr/D;
    dB1_U1D0[ijk] = Omega_BHNS;
  }
  
  DECLARE_AND_EMPTY_FIELD(ddB1_U1D0D0)
  DECLARE_AND_EMPTY_FIELD(ddB1_U1D1D2)
  DECLARE_AND_EMPTY_FIELD(ddB1_U1D1D1)
  DECLARE_AND_EMPTY_FIELD(ddB1_U1D0D1)
  DECLARE_AND_EMPTY_FIELD(ddB1_U1D2D2)
  DECLARE_AND_EMPTY_FIELD(ddB1_U1D0D2)
  
  /* since the second derivative is zero */
  ddB1_U1D0D0->v = alloc_double(nn);
  ddB1_U1D1D2->v = alloc_double(nn);
  ddB1_U1D1D1->v = alloc_double(nn);
  ddB1_U1D0D1->v = alloc_double(nn);
  ddB1_U1D2D2->v = alloc_double(nn);
  ddB1_U1D0D2->v = alloc_double(nn);
}

/* updating derivative */
void bbn_update_derivative_B1_U0(Patch_T *const patch)
{
  const double D          = GetParameterD_E("BH_NS_separation");
  const double Vr         = GetParameterD_E("BH_NS_infall_velocity");
  const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity"); 
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  PREP_FIELD(dB1_U0D2)
  PREP_FIELD(dB1_U0D1)
  PREP_FIELD(dB1_U0D0)
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    dB1_U0D2[ijk] = 0;
    dB1_U0D1[ijk] = -Omega_BHNS;
    dB1_U0D0[ijk] = Vr/D;
  }
  
  DECLARE_AND_EMPTY_FIELD(ddB1_U0D0D0)
  DECLARE_AND_EMPTY_FIELD(ddB1_U0D1D2)
  DECLARE_AND_EMPTY_FIELD(ddB1_U0D1D1)
  DECLARE_AND_EMPTY_FIELD(ddB1_U0D0D1)
  DECLARE_AND_EMPTY_FIELD(ddB1_U0D2D2)
  DECLARE_AND_EMPTY_FIELD(ddB1_U0D0D2)
  
  /* since the second derivative is zero */
  ddB1_U0D0D0->v = alloc_double(nn);
  ddB1_U0D1D2->v = alloc_double(nn);
  ddB1_U0D1D1->v = alloc_double(nn);
  ddB1_U0D0D1->v = alloc_double(nn);
  ddB1_U0D2D2->v = alloc_double(nn);
  ddB1_U0D0D2->v = alloc_double(nn);
}


/* updating derivative. 
// NOTE: make sure dB0,ddB0,dB1,ddB1 are already updated.*/
void bbn_update_derivative_Beta_U0(Patch_T *const patch)
{
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  bbn_update_derivative_B0_U0(patch);
  bbn_update_derivative_B1_U0(patch);  
  
  GET_FIELD(dB0_U0D1)
  GET_FIELD(dB0_U0D0)
  GET_FIELD(dB0_U0D2)
  GET_FIELD(ddB0_U0D2D2)
  GET_FIELD(ddB0_U0D0D1)
  GET_FIELD(ddB0_U0D1D1)
  GET_FIELD(ddB0_U0D1D2)
  GET_FIELD(ddB0_U0D0D2)
  GET_FIELD(ddB0_U0D0D0)

  GET_FIELD(dB1_U0D1)
  GET_FIELD(dB1_U0D0)
  GET_FIELD(dB1_U0D2)
  GET_FIELD(ddB1_U0D2D2)
  GET_FIELD(ddB1_U0D0D1)
  GET_FIELD(ddB1_U0D1D1)
  GET_FIELD(ddB1_U0D1D2)
  GET_FIELD(ddB1_U0D0D2)
  GET_FIELD(ddB1_U0D0D0)
    
  PREP_FIELD(dBeta_U0D1)
  PREP_FIELD(dBeta_U0D0)
  PREP_FIELD(dBeta_U0D2)
  PREP_FIELD(ddBeta_U0D2D2)
  PREP_FIELD(ddBeta_U0D0D1)
  PREP_FIELD(ddBeta_U0D1D1)
  PREP_FIELD(ddBeta_U0D1D2)
  PREP_FIELD(ddBeta_U0D0D2)
  PREP_FIELD(ddBeta_U0D0D0)

  for (ijk = 0; ijk < nn; ++ijk)
  {
    dBeta_U0D1[ijk]    = dB0_U0D1[ijk]+dB1_U0D1[ijk];
    dBeta_U0D0[ijk]    = dB0_U0D0[ijk]+dB1_U0D0[ijk];
    dBeta_U0D2[ijk]    = dB0_U0D2[ijk]+dB1_U0D2[ijk];
    ddBeta_U0D2D2[ijk] = ddB0_U0D2D2[ijk]+ddB1_U0D2D2[ijk];
    ddBeta_U0D0D1[ijk] = ddB0_U0D0D1[ijk]+ddB1_U0D0D1[ijk];
    ddBeta_U0D1D1[ijk] = ddB0_U0D1D1[ijk]+ddB1_U0D1D1[ijk];
    ddBeta_U0D1D2[ijk] = ddB0_U0D1D2[ijk]+ddB1_U0D1D2[ijk];
    ddBeta_U0D0D2[ijk] = ddB0_U0D0D2[ijk]+ddB1_U0D0D2[ijk];
    ddBeta_U0D0D0[ijk] = ddB0_U0D0D0[ijk]+ddB1_U0D0D0[ijk];
  }

}

/* updating derivative. 
// NOTE: make sure dB0,ddB0,dB1,ddB1 are already updated.*/
void bbn_update_derivative_Beta_U1(Patch_T *const patch)
{
  const unsigned nn = patch->nn;
  unsigned ijk;
  bbn_update_derivative_B0_U1(patch);
  bbn_update_derivative_B1_U1(patch);
      
  GET_FIELD(dB0_U1D1)
  GET_FIELD(dB0_U1D0)
  GET_FIELD(dB0_U1D2)
  GET_FIELD(ddB0_U1D2D2)
  GET_FIELD(ddB0_U1D0D1)
  GET_FIELD(ddB0_U1D1D1)
  GET_FIELD(ddB0_U1D1D2)
  GET_FIELD(ddB0_U1D0D2)
  GET_FIELD(ddB0_U1D0D0)

  GET_FIELD(dB1_U1D1)
  GET_FIELD(dB1_U1D0)
  GET_FIELD(dB1_U1D2)
  GET_FIELD(ddB1_U1D2D2)
  GET_FIELD(ddB1_U1D0D1)
  GET_FIELD(ddB1_U1D1D1)
  GET_FIELD(ddB1_U1D1D2)
  GET_FIELD(ddB1_U1D0D2)
  GET_FIELD(ddB1_U1D0D0)
    
  PREP_FIELD(dBeta_U1D1)
  PREP_FIELD(dBeta_U1D0)
  PREP_FIELD(dBeta_U1D2)
  PREP_FIELD(ddBeta_U1D2D2)
  PREP_FIELD(ddBeta_U1D0D1)
  PREP_FIELD(ddBeta_U1D1D1)
  PREP_FIELD(ddBeta_U1D1D2)
  PREP_FIELD(ddBeta_U1D0D2)
  PREP_FIELD(ddBeta_U1D0D0)

  for (ijk = 0; ijk < nn; ++ijk)
  {
    dBeta_U1D1[ijk]    = dB0_U1D1[ijk]+dB1_U1D1[ijk];
    dBeta_U1D0[ijk]    = dB0_U1D0[ijk]+dB1_U1D0[ijk];
    dBeta_U1D2[ijk]    = dB0_U1D2[ijk]+dB1_U1D2[ijk];
    ddBeta_U1D2D2[ijk] = ddB0_U1D2D2[ijk]+ddB1_U1D2D2[ijk];
    ddBeta_U1D0D1[ijk] = ddB0_U1D0D1[ijk]+ddB1_U1D0D1[ijk];
    ddBeta_U1D1D1[ijk] = ddB0_U1D1D1[ijk]+ddB1_U1D1D1[ijk];
    ddBeta_U1D1D2[ijk] = ddB0_U1D1D2[ijk]+ddB1_U1D1D2[ijk];
    ddBeta_U1D0D2[ijk] = ddB0_U1D0D2[ijk]+ddB1_U1D0D2[ijk];
    ddBeta_U1D0D0[ijk] = ddB0_U1D0D0[ijk]+ddB1_U1D0D0[ijk];
  }

}


/* updating derivative. 
// NOTE: make sure dB0,ddB0,dB1,ddB1 are already updated.*/
void bbn_update_derivative_Beta_U2(Patch_T *const patch)
{
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  bbn_update_derivative_B0_U2(patch);
  bbn_update_derivative_B1_U2(patch);
  
  GET_FIELD(dB0_U2D1)
  GET_FIELD(dB0_U2D0)
  GET_FIELD(dB0_U2D2)
  GET_FIELD(ddB0_U2D2D2)
  GET_FIELD(ddB0_U2D0D1)
  GET_FIELD(ddB0_U2D1D1)
  GET_FIELD(ddB0_U2D1D2)
  GET_FIELD(ddB0_U2D0D2)
  GET_FIELD(ddB0_U2D0D0)

  GET_FIELD(dB1_U2D1)
  GET_FIELD(dB1_U2D0)
  GET_FIELD(dB1_U2D2)
  GET_FIELD(ddB1_U2D2D2)
  GET_FIELD(ddB1_U2D0D1)
  GET_FIELD(ddB1_U2D1D1)
  GET_FIELD(ddB1_U2D1D2)
  GET_FIELD(ddB1_U2D0D2)
  GET_FIELD(ddB1_U2D0D0)
    
  PREP_FIELD(dBeta_U2D1)
  PREP_FIELD(dBeta_U2D0)
  PREP_FIELD(dBeta_U2D2)
  PREP_FIELD(ddBeta_U2D2D2)
  PREP_FIELD(ddBeta_U2D0D1)
  PREP_FIELD(ddBeta_U2D1D1)
  PREP_FIELD(ddBeta_U2D1D2)
  PREP_FIELD(ddBeta_U2D0D2)
  PREP_FIELD(ddBeta_U2D0D0)

  for (ijk = 0; ijk < nn; ++ijk)
  {
    dBeta_U2D1[ijk]    = dB0_U2D1[ijk]+dB1_U2D1[ijk];
    dBeta_U2D0[ijk]    = dB0_U2D0[ijk]+dB1_U2D0[ijk];
    dBeta_U2D2[ijk]    = dB0_U2D2[ijk]+dB1_U2D2[ijk];
    ddBeta_U2D2D2[ijk] = ddB0_U2D2D2[ijk]+ddB1_U2D2D2[ijk];
    ddBeta_U2D0D1[ijk] = ddB0_U2D0D1[ijk]+ddB1_U2D0D1[ijk];
    ddBeta_U2D1D1[ijk] = ddB0_U2D1D1[ijk]+ddB1_U2D1D1[ijk];
    ddBeta_U2D1D2[ijk] = ddB0_U2D1D2[ijk]+ddB1_U2D1D2[ijk];
    ddBeta_U2D0D2[ijk] = ddB0_U2D0D2[ijk]+ddB1_U2D0D2[ijk];
    ddBeta_U2D0D0[ijk] = ddB0_U2D0D0[ijk]+ddB1_U2D0D0[ijk];
  }

}

/* after finding new NS surface, root finder might find h ~ 1
// so put it h = 1, to prevent nan in matter fields  */
static void cleaning_enthalpy(Patch_T *const patch)
{
  if(!IsItNSSurface(patch))
    return;
    
  GET_FIELD(enthalpy)
  
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    if (LSSEQL(enthalpy[ijk],1))
      enthalpy[ijk] = 1;
  }

}

/* updating enthalpy and its derivative */
void bbn_update_enthalpy_and_denthalpy(Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if(!IsItNSPatch(patch))
      continue;
      
    Tij_IF_CTS_enthalpy(patch);
    bbn_update_derivative_enthalpy(patch);  
  }

}

/* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S
// which used in stress energy tensor in relaxed fashion.
// furthermore, it updates the central density parameter. 
// flag = 1 puts enthalpy = 1 if it is less than 1.
// note: dphi^i and W^i are assumed ready. */
void bbn_update_stress_energy_tensor(Grid_T *const grid,const int flag)
{
  assert(grid);
  
  pr_line_custom('=');
  printf("Updating enthalpy, rest-mass density and their derivatives ...\n");
  
  const double W1  = GetParameterD_E("Solving_Field_Update_Weight");
  const double W2  = 1-W1;
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    if(!IsItNSPatch(patch))
      continue;
      
    DECLARE_FIELD(enthalpy)
    DECLARE_FIELD(OLD_enthalpy)
    free_coeffs(enthalpy);
    free_coeffs(OLD_enthalpy);
    
    /* swapping pointers */
    void *tmp       = OLD_enthalpy->v;
    OLD_enthalpy->v = enthalpy->v;
    enthalpy->v     = tmp;
    
    /* update enthalpy */
    Tij_IF_CTS_enthalpy(patch);
    for (ijk = 0; ijk < nn; ++ijk)
      enthalpy->v[ijk] = W1*enthalpy->v[ijk]+W2*OLD_enthalpy->v[ijk];

    if (flag == 1)
      cleaning_enthalpy(patch);
      
    bbn_update_derivative_enthalpy(patch);
    bbn_update_rho0(patch);
    bbn_update_derivative_rho0(patch);
  
  }
  
  /* update central density parameter */
  Patch_T *patch = GetPatch("left_centeral_box",grid);
  DECLARE_FIELD(rho0);
  double rho_center;
  const double *x = patch->c;/* NS center */
  double X[3] = {0};
  
  X_of_x(X,x,patch);
  Interpolation_T *interp_s = init_interpolation();
  interp_s->field = rho0;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  rho_center = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  update_parameter_double_format("rho_center",rho_center);
  
  printf("Updating enthalpy, rest-mass density and their derivatives ==> Done.\n");
  pr_clock();
  pr_line_custom('=');

  Tij_IF_CTS_psi6Sources(grid);
}

/* updaing B1 */
void bbn_update_B1_U012(Patch_T *const patch)
{
  const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity");
  const double Vr   = GetParameterD_E("BH_NS_infall_velocity");
  const double D    = GetParameterD_E("BH_NS_separation");
  const double y_CM = GetParameterD_E("y_CM");
  const double x_CM = GetParameterD_E("x_CM");
  const unsigned nn = patch->nn;
  unsigned ijk;
    
  /* B^1 */
  PREP_FIELD(B1_U0)
  PREP_FIELD(B1_U1)
  PREP_FIELD(B1_U2)
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double x     = patch->node[ijk]->x[0];
    double y     = patch->node[ijk]->x[1];
    
    B1_U0[ijk] = Omega_BHNS*(-y+y_CM)+Vr*(x-x_CM)/D;
    B1_U1[ijk] = Omega_BHNS*(x-x_CM)+Vr*(y-y_CM)/D;
    B1_U2[ijk] = 0;
  }

}
