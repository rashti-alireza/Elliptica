/*
// Alireza Rashti
// August 2019
*/


#include "bbn_fields.h"

/* declaring all of the fields needed for construction of Initial Data
// in the given grid */
void bbn_add_fields(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Adding fields ...\n");
  
  /* if it is ready */
  if (Pgeti("use_previous_data"))
  {
    printf("~> Using the fields of the previous grid.\n");
    printf("} Adding fields ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
    return;
  }
  
  unsigned p;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    bbn_add_fields_in_patch(patch);
  }
  
  printf("} Adding fields ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* declaring all of the fields needed for construction of Initial Data
// in the given patch */
void bbn_add_fields_in_patch(Patch_T *const patch)
{
  /* scalar for the irrotational part of fluid i.e h*u = dphi+W 
  // in NS and its partial derivatives*/
  ADD_FIELD(phi);
  ADD_FIELD(phi_residual);
  
  ADD_FIELD(OLD_phi);
  ADD_FIELD(dphi_D2)
  ADD_FIELD(dphi_D1)
  ADD_FIELD(dphi_D0)
  
  ADD_FIELD(ddphi_D2D2)
  ADD_FIELD(ddphi_D1D2)
  ADD_FIELD(ddphi_D1D1)
  ADD_FIELD(ddphi_D0D2)
  ADD_FIELD(ddphi_D0D0)
  ADD_FIELD(ddphi_D0D1)

  /* enthalpy in NS and its partial derivatives */
  ADD_FIELD(enthalpy);
  ADD_FIELD(denthalpy_D2)
  ADD_FIELD(denthalpy_D1)
  ADD_FIELD(denthalpy_D0)

  /* rest mass density in NS and its partial derivatives */
  ADD_FIELD(rho0);
  ADD_FIELD(drho0_D2)
  ADD_FIELD(drho0_D1)
  ADD_FIELD(drho0_D0)
  
  /* the first component of fluid four velocity, i.e. 
  // u_mu = (u_U0,u_U1,u_U2,u_U3) and its partial derivatives */
  ADD_FIELD(u0);
  ADD_FIELD(du0_D2)
  ADD_FIELD(du0_D1)
  ADD_FIELD(du0_D0)
  
  /* spin part of fluid W^i */
  ADD_FIELD(W_U0);
  ADD_FIELD(W_U1);
  ADD_FIELD(W_U2);
    
  /* normal vector on horizon */
  ADD_FIELD(_HS_U0);
  ADD_FIELD(_HS_U1);
  ADD_FIELD(_HS_U2);
  ADD_FIELD(_dHS_U1D0)
  ADD_FIELD(_dHS_U1D1)
  ADD_FIELD(_dHS_U1D2)
  ADD_FIELD(_dHS_U0D1)
  ADD_FIELD(_dHS_U0D0)
  ADD_FIELD(_dHS_U0D2)
  ADD_FIELD(_dHS_U2D2)
  ADD_FIELD(_dHS_U2D1)
  ADD_FIELD(_dHS_U2D0)
  
  /* Hamiltonian and Momentum constraints */
  ADD_FIELD(ham_constraint);
  ADD_FIELD(mom_constraint_U0);
  ADD_FIELD(mom_constraint_U1);
  ADD_FIELD(mom_constraint_U2);

  /* Hamiltonian and Momentum constraints for 
  // second method of calculations */ 
  ADD_FIELD(ham_constraint_2nd);
  ADD_FIELD(mom_constraint_2nd_U0);
  ADD_FIELD(mom_constraint_2nd_U1);
  ADD_FIELD(mom_constraint_2nd_U2);
  
  /* conformal total energy density */
  ADD_FIELD(_E);
  
  /* conformal trace of stress tensor */
  ADD_FIELD(_S);
  
  /* conformal momentum current */
  ADD_FIELD(_J_U0);
  ADD_FIELD(_J_U1);
  ADD_FIELD(_J_U2);
  
  /* conformal factor and its derivative */
  ADD_FIELD(psi);
  ADD_FIELD(psi_residual);
  ADD_FIELD(OLD_psi);
  ADD_FIELD(dpsi_D2)
  ADD_FIELD(dpsi_D1)
  ADD_FIELD(dpsi_D0)
  
  ADD_FIELD(ddpsi_D2D2)
  ADD_FIELD(ddpsi_D1D2)
  ADD_FIELD(ddpsi_D1D1)
  ADD_FIELD(ddpsi_D0D2)
  ADD_FIELD(ddpsi_D0D0)
  ADD_FIELD(ddpsi_D0D1)

  /* eta = lapse * psi and its partial derivatives */
  ADD_FIELD(eta);
  ADD_FIELD(eta_residual);
  
  ADD_FIELD(OLD_eta)
  
  ADD_FIELD(deta_D2)
  ADD_FIELD(deta_D1)
  ADD_FIELD(deta_D0)

  ADD_FIELD(ddeta_D2D2)
  ADD_FIELD(ddeta_D1D2)
  ADD_FIELD(ddeta_D1D1)
  ADD_FIELD(ddeta_D0D2)
  ADD_FIELD(ddeta_D0D0)
  ADD_FIELD(ddeta_D0D1)    
  
  /* shift, Beta^i = B0^i+B1^i, 
  // B1^i = Omega_BHNS*(-y+y_CM,x-x_CM,0)+v_r/D*(x-x_CM,y-y_CM) 
  // and its partial derivative */
  ADD_FIELD(B0_U0);
  ADD_FIELD(B0_U1);
  ADD_FIELD(B0_U2);
  ADD_FIELD(B0_U0_residual);
  ADD_FIELD(B0_U1_residual);
  ADD_FIELD(B0_U2_residual);
  
  ADD_FIELD(B1_U0);
  ADD_FIELD(B1_U1);
  ADD_FIELD(B1_U2);
  ADD_FIELD(Beta_U0);
  ADD_FIELD(Beta_U1);
  ADD_FIELD(Beta_U2);
  ADD_FIELD(OLD_B0_U0);
  ADD_FIELD(OLD_B0_U1);
  ADD_FIELD(OLD_B0_U2);
  
  ADD_FIELD(dB0_U1D0)
  ADD_FIELD(dB0_U1D1)
  ADD_FIELD(dB0_U1D2)
  ADD_FIELD(dB0_U0D1)
  ADD_FIELD(dB0_U0D0)
  ADD_FIELD(dB0_U0D2)
  ADD_FIELD(dB0_U2D2)
  ADD_FIELD(dB0_U2D1)
  ADD_FIELD(dB0_U2D0)

  ADD_FIELD(ddB0_U1D2D2)
  ADD_FIELD(ddB0_U1D0D0)
  ADD_FIELD(ddB0_U1D0D1)
  ADD_FIELD(ddB0_U1D0D2)
  ADD_FIELD(ddB0_U1D1D2)
  ADD_FIELD(ddB0_U2D0D0)
  ADD_FIELD(ddB0_U2D1D2)
  ADD_FIELD(ddB0_U0D2D2)
  ADD_FIELD(ddB0_U2D1D1)
  ADD_FIELD(ddB0_U0D0D1)
  ADD_FIELD(ddB0_U0D1D1)
  ADD_FIELD(ddB0_U0D1D2)
  ADD_FIELD(ddB0_U0D0D2)
  ADD_FIELD(ddB0_U1D1D1)
  ADD_FIELD(ddB0_U2D0D1)
  ADD_FIELD(ddB0_U2D2D2)
  ADD_FIELD(ddB0_U2D0D2)
  ADD_FIELD(ddB0_U0D0D0)

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
  ADD_FIELD(_dGamma_U2D2D2D2)
  ADD_FIELD(_dGamma_U2D2D2D0)
  ADD_FIELD(_dGamma_U2D2D2D1)
  ADD_FIELD(_dGamma_U2D0D0D2)
  ADD_FIELD(_dGamma_U1D1D2D2)
  ADD_FIELD(_dGamma_U2D0D0D0)
  ADD_FIELD(_dGamma_U1D1D2D0)
  ADD_FIELD(_dGamma_U1D1D2D1)
  ADD_FIELD(_dGamma_U2D1D1D0)
  ADD_FIELD(_dGamma_U2D0D0D1)
  ADD_FIELD(_dGamma_U2D0D2D1)
  ADD_FIELD(_dGamma_U1D0D1D0)
  ADD_FIELD(_dGamma_U1D0D1D1)
  ADD_FIELD(_dGamma_U1D0D1D2)
  ADD_FIELD(_dGamma_U1D2D2D1)
  ADD_FIELD(_dGamma_U1D0D0D1)
  ADD_FIELD(_dGamma_U1D0D0D0)
  ADD_FIELD(_dGamma_U1D0D0D2)
  ADD_FIELD(_dGamma_U0D1D2D2)
  ADD_FIELD(_dGamma_U0D1D2D1)
  ADD_FIELD(_dGamma_U0D1D2D0)
  ADD_FIELD(_dGamma_U2D0D2D0)
  ADD_FIELD(_dGamma_U1D0D2D2)
  ADD_FIELD(_dGamma_U1D0D2D1)
  ADD_FIELD(_dGamma_U1D0D2D0)
  ADD_FIELD(_dGamma_U2D1D1D2)
  ADD_FIELD(_dGamma_U2D0D2D2)
  ADD_FIELD(_dGamma_U0D0D1D0)
  ADD_FIELD(_dGamma_U1D2D2D0)
  ADD_FIELD(_dGamma_U2D1D2D1)
  ADD_FIELD(_dGamma_U2D0D1D2)
  ADD_FIELD(_dGamma_U2D0D1D1)
  ADD_FIELD(_dGamma_U2D0D1D0)
  ADD_FIELD(_dGamma_U2D1D2D2)
  ADD_FIELD(_dGamma_U0D1D1D0)
  ADD_FIELD(_dGamma_U0D1D1D1)
  ADD_FIELD(_dGamma_U0D1D1D2)
  ADD_FIELD(_dGamma_U1D2D2D2)
  ADD_FIELD(_dGamma_U1D1D1D1)
  ADD_FIELD(_dGamma_U1D1D1D0)
  ADD_FIELD(_dGamma_U1D1D1D2)
  ADD_FIELD(_dGamma_U0D0D1D1)
  ADD_FIELD(_dGamma_U0D0D2D2)
  ADD_FIELD(_dGamma_U0D0D2D0)
  ADD_FIELD(_dGamma_U0D0D2D1)
  ADD_FIELD(_dGamma_U2D1D2D0)
  ADD_FIELD(_dGamma_U0D0D0D0)
  ADD_FIELD(_dGamma_U0D0D0D1)
  ADD_FIELD(_dGamma_U0D0D0D2)
  ADD_FIELD(_dGamma_U2D1D1D1)
  ADD_FIELD(_dGamma_U0D2D2D0)
  ADD_FIELD(_dGamma_U0D2D2D1)
  ADD_FIELD(_dGamma_U0D2D2D2)
  ADD_FIELD(_dGamma_U0D0D1D2)

  /* Ricci made up of conformal metric _gamma */
  ADD_FIELD(_R);
  ADD_FIELD(_Ric_D1D2)
  ADD_FIELD(_Ric_D2D2)
  ADD_FIELD(_Ric_D0D2)
  ADD_FIELD(_Ric_D0D1)
  ADD_FIELD(_Ric_D0D0)
  ADD_FIELD(_Ric_D1D1)

  /* extrinsic curvature */
  ADD_FIELD(K);
  ADD_FIELD(dK_D2)
  ADD_FIELD(dK_D1)
  ADD_FIELD(dK_D0)

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

/* taking partial derivatives of fields needed in equations */
void bbn_partial_derivatives_fields(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Taking partial derivatives of fields ...\n");
  
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
  
  printf("} Taking partial derivatives of fields ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* updating Beta_U0 */
void bbn_update_Beta_U0(Patch_T *const patch)
{
  REALLOC_v_WRITE_v(Beta_U0)
  READ_v(B0_U0)
  READ_v(B1_U0)
  
  unsigned ijk;
  FOR_ALL_POINTS(ijk,patch)
  {
    Beta_U0[ijk] = B0_U0[ijk]+B1_U0[ijk];
  }
  
}

/* updating Beta_U1 */
void bbn_update_Beta_U1(Patch_T *const patch)
{
  REALLOC_v_WRITE_v(Beta_U1)
  READ_v(B0_U1)
  READ_v(B1_U1)
  
  unsigned ijk;
  FOR_ALL_POINTS(ijk,patch)
  {
    Beta_U1[ijk] = B0_U1[ijk]+B1_U1[ijk];
  }
  
}

/* updating Beta_U2 */
void bbn_update_Beta_U2(Patch_T *const patch)
{
  REALLOC_v_WRITE_v(Beta_U2)
  READ_v(B0_U2)
  READ_v(B1_U2)
  
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
    READ_v(enthalpy)
    REALLOC_v_WRITE_v(rho0)
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
    if (!rho0->v) return;/* if rho0 has not been set yet */
    
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
    if (!u0->v) return;/* if u0 has not been set yet */
    
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

/* add ddK fields and take 2nd order derivatives of K 
// mainly needed for extrapolation of initiala data at BH filler. */
void bbn_add_and_take_2nd_derivatives_K(Patch_T *const patch)
{
  ADD_FIELD(ddK_D2D2)
  ADD_FIELD(ddK_D1D2)
  ADD_FIELD(ddK_D1D1)
  ADD_FIELD(ddK_D0D2)
  ADD_FIELD(ddK_D0D0)
  ADD_FIELD(ddK_D0D1)
  
  DECLARE_FIELD(ddK_D2D2)
  DECLARE_FIELD(ddK_D1D2)
  DECLARE_FIELD(ddK_D1D1)
  DECLARE_FIELD(ddK_D0D2)
  DECLARE_FIELD(ddK_D0D0)
  DECLARE_FIELD(ddK_D0D1)
  
  DECLARE_FIELD(dK_D2)
  DECLARE_FIELD(dK_D1)
  DECLARE_FIELD(dK_D0)
  
  ddK_D2D2->v = Partial_Derivative(dK_D2,"z");
  ddK_D1D2->v = Partial_Derivative(dK_D1,"z");
  ddK_D1D1->v = Partial_Derivative(dK_D1,"y");
  ddK_D0D2->v = Partial_Derivative(dK_D0,"z");
  ddK_D0D0->v = Partial_Derivative(dK_D0,"x");
  ddK_D0D1->v = Partial_Derivative(dK_D0,"y");
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
  //const double D          = Pgetd("BH_NS_separation");
  //const double Vr         = Pgetd("BH_NS_infall_velocity");
  //const double Omega_BHNS = Pgetd("BH_NS_angular_velocity"); 
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  REALLOC_v_WRITE_v(dB1_U2D2)
  REALLOC_v_WRITE_v(dB1_U2D1)
  REALLOC_v_WRITE_v(dB1_U2D0)
  
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
  const double D          = Pgetd("BH_NS_separation");
  const double Vr         = Pgetd("BH_NS_infall_velocity");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity"); 
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  REALLOC_v_WRITE_v(dB1_U1D2)
  REALLOC_v_WRITE_v(dB1_U1D1)
  REALLOC_v_WRITE_v(dB1_U1D0)
  
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
  const double D          = Pgetd("BH_NS_separation");
  const double Vr         = Pgetd("BH_NS_infall_velocity");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity"); 
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  REALLOC_v_WRITE_v(dB1_U0D2)
  REALLOC_v_WRITE_v(dB1_U0D1)
  REALLOC_v_WRITE_v(dB1_U0D0)
  
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
  
  READ_v(dB0_U0D1)
  READ_v(dB0_U0D0)
  READ_v(dB0_U0D2)
  READ_v(ddB0_U0D2D2)
  READ_v(ddB0_U0D0D1)
  READ_v(ddB0_U0D1D1)
  READ_v(ddB0_U0D1D2)
  READ_v(ddB0_U0D0D2)
  READ_v(ddB0_U0D0D0)

  READ_v(dB1_U0D1)
  READ_v(dB1_U0D0)
  READ_v(dB1_U0D2)
  READ_v(ddB1_U0D2D2)
  READ_v(ddB1_U0D0D1)
  READ_v(ddB1_U0D1D1)
  READ_v(ddB1_U0D1D2)
  READ_v(ddB1_U0D0D2)
  READ_v(ddB1_U0D0D0)
    
  REALLOC_v_WRITE_v(dBeta_U0D1)
  REALLOC_v_WRITE_v(dBeta_U0D0)
  REALLOC_v_WRITE_v(dBeta_U0D2)
  REALLOC_v_WRITE_v(ddBeta_U0D2D2)
  REALLOC_v_WRITE_v(ddBeta_U0D0D1)
  REALLOC_v_WRITE_v(ddBeta_U0D1D1)
  REALLOC_v_WRITE_v(ddBeta_U0D1D2)
  REALLOC_v_WRITE_v(ddBeta_U0D0D2)
  REALLOC_v_WRITE_v(ddBeta_U0D0D0)

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
      
  READ_v(dB0_U1D1)
  READ_v(dB0_U1D0)
  READ_v(dB0_U1D2)
  READ_v(ddB0_U1D2D2)
  READ_v(ddB0_U1D0D1)
  READ_v(ddB0_U1D1D1)
  READ_v(ddB0_U1D1D2)
  READ_v(ddB0_U1D0D2)
  READ_v(ddB0_U1D0D0)

  READ_v(dB1_U1D1)
  READ_v(dB1_U1D0)
  READ_v(dB1_U1D2)
  READ_v(ddB1_U1D2D2)
  READ_v(ddB1_U1D0D1)
  READ_v(ddB1_U1D1D1)
  READ_v(ddB1_U1D1D2)
  READ_v(ddB1_U1D0D2)
  READ_v(ddB1_U1D0D0)
    
  REALLOC_v_WRITE_v(dBeta_U1D1)
  REALLOC_v_WRITE_v(dBeta_U1D0)
  REALLOC_v_WRITE_v(dBeta_U1D2)
  REALLOC_v_WRITE_v(ddBeta_U1D2D2)
  REALLOC_v_WRITE_v(ddBeta_U1D0D1)
  REALLOC_v_WRITE_v(ddBeta_U1D1D1)
  REALLOC_v_WRITE_v(ddBeta_U1D1D2)
  REALLOC_v_WRITE_v(ddBeta_U1D0D2)
  REALLOC_v_WRITE_v(ddBeta_U1D0D0)

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
  
  READ_v(dB0_U2D1)
  READ_v(dB0_U2D0)
  READ_v(dB0_U2D2)
  READ_v(ddB0_U2D2D2)
  READ_v(ddB0_U2D0D1)
  READ_v(ddB0_U2D1D1)
  READ_v(ddB0_U2D1D2)
  READ_v(ddB0_U2D0D2)
  READ_v(ddB0_U2D0D0)

  READ_v(dB1_U2D1)
  READ_v(dB1_U2D0)
  READ_v(dB1_U2D2)
  READ_v(ddB1_U2D2D2)
  READ_v(ddB1_U2D0D1)
  READ_v(ddB1_U2D1D1)
  READ_v(ddB1_U2D1D2)
  READ_v(ddB1_U2D0D2)
  READ_v(ddB1_U2D0D0)
    
  REALLOC_v_WRITE_v(dBeta_U2D1)
  REALLOC_v_WRITE_v(dBeta_U2D0)
  REALLOC_v_WRITE_v(dBeta_U2D2)
  REALLOC_v_WRITE_v(ddBeta_U2D2D2)
  REALLOC_v_WRITE_v(ddBeta_U2D0D1)
  REALLOC_v_WRITE_v(ddBeta_U2D1D1)
  REALLOC_v_WRITE_v(ddBeta_U2D1D2)
  REALLOC_v_WRITE_v(ddBeta_U2D0D2)
  REALLOC_v_WRITE_v(ddBeta_U2D0D0)

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
    
  WRITE_v(enthalpy)
  
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
  printf("{ Updating enthalpy, rest-mass density and their derivatives ...\n");
  
  const double W1  = Pgetd("NS_enthalpy_update_weight");
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
    free_coeffs(enthalpy);
    
    /* if there is any enthalpy it is legitimate for relaxed update */
    if (enthalpy->v)
    {
      double *old_value = enthalpy->v;
      enthalpy->v       = 0;
      
      /* update enthalpy */
      Tij_IF_CTS_enthalpy(patch);
      for (ijk = 0; ijk < nn; ++ijk)
        enthalpy->v[ijk] = W1*enthalpy->v[ijk]+W2*old_value[ijk];
      
      free(old_value);
    }
    else/* if there is no previous enthalpy */
    {
      /* update enthalpy */
      Tij_IF_CTS_enthalpy(patch);
    }

    if (flag == 1)
      cleaning_enthalpy(patch);
      
    bbn_update_derivative_enthalpy(patch);
    bbn_update_rho0(patch);
    bbn_update_derivative_rho0(patch);
  
  }
  
  /* update central quantites */
  Patch_T *patch = GetPatch("left_central_box",grid);
  DECLARE_FIELD(enthalpy)
  double h_center;
  const double *x = patch->c;/* NS center */
  double X[3] = {0};
  
  X_of_x(X,x,patch);
  EoS_T *eos = initialize_EoS();
  Interpolation_T *interp_s = init_interpolation();
  interp_s->field = enthalpy;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  h_center = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  eos->h = h_center;
  Psetd("rho_center",eos->rest_mass_density(eos));
  Psetd("pressure_center",eos->pressure(eos));
  Psetd("energy_density_center",eos->energy_density(eos));
  free_EoS(eos);
  
  printf("} Updating enthalpy, rest-mass density and their derivatives ==> Done.\n");
  pr_clock();
  pr_line_custom('=');

  Tij_IF_CTS_psi6Sources(grid);
  
  /* update u0 derivatives */
  FOR_ALL_PATCHES(p,grid)
  {
    patch = grid->patch[p];
    bbn_update_derivative_u0(patch);
  }
    
}

/* updaing B1 */
void bbn_update_B1_U012(Patch_T *const patch)
{
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double Vr   = Pgetd("BH_NS_infall_velocity");
  const double D    = Pgetd("BH_NS_separation");
  const double y_CM = Pgetd("y_CM");
  const double x_CM = Pgetd("x_CM");
  const unsigned nn = patch->nn;
  unsigned ijk;
    
  /* B^1 */
  REALLOC_v_WRITE_v(B1_U0)
  REALLOC_v_WRITE_v(B1_U1)
  REALLOC_v_WRITE_v(B1_U2)
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double x     = patch->node[ijk]->x[0];
    double y     = patch->node[ijk]->x[1];
    
    B1_U0[ijk] = Omega_BHNS*(-y+y_CM)+Vr*(x-x_CM)/D;
    B1_U1[ijk] = Omega_BHNS*(x-x_CM)+Vr*(y-y_CM)/D;
    B1_U2[ijk] = 0;
  }

}

/* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
// _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
void bbn_update_Aij(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Updating _A^{ij}, _dA^{ij} and _A^{ij}*A_{ij} ...\n");
  unsigned p;

  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (IsItInsideBHPatch(patch))
      continue;
      
    bbn_update_psi10A_UiUj(patch);
  }
  
  printf("} Updating _A^{ij}, _dA^{ij} and _A^{ij}*A_{ij} ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}
