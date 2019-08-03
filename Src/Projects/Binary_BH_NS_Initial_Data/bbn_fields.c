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
      /* scalar for the irrotational part of fluid i.e h*u = dphi+W in NS */
      add_field("phi",0,patch,YES);
      
      /* enthalpy in NS */
      add_field("enthalpy",0,patch,YES);
      
      /* rest mass density in NS */
      add_field("rho0",0,patch,YES);
      
      /* the first component of fluid four velocity, i.e. 
      // u_mu = (u_U0,u_U1,u_U2,u_U3) */
      add_field("u0",0,patch,YES);
      
      /* spin part of fluid W^i */
      add_field("W_U0",0,patch,YES);
      add_field("W_U1",0,patch,YES);
      add_field("W_U2",0,patch,YES);
      
      /* conformal total energy density */
      add_field("_E",0,patch,YES);
      
      /* conformal trace of stress tensor */
      add_field("_S",0,patch,YES);
      
      /* conformal momentum current */
      add_field("_J_U0",0,patch,YES);
      add_field("_J_U1",0,patch,YES);
      add_field("_J_U2",0,patch,YES);
    
    }
    
    /* conformal factor */
    add_field("psi",0,patch,YES);
    
    /* eta = lapse * psi */
    add_field("eta",0,patch,YES);
    
    /* shift, Beta^i = B0^i+B1^i, B1^i = omega*(-y+y_CM,x,0)+v_r/D*(x,y-y_CM) */
    add_field("B0_U0",0,patch,YES);
    add_field("B0_U1",0,patch,YES);
    add_field("B0_U2",0,patch,YES);
    add_field("B1_U0",0,patch,YES);
    add_field("B1_U1",0,patch,YES);
    add_field("B1_U2",0,patch,YES);
    add_field("Beta_U0",0,patch,YES);
    add_field("Beta_U1",0,patch,YES);
    add_field("Beta_U2",0,patch,YES);
    
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

    /* Ricci scalar made up of conformal metric _gamma */
    add_field("_R",0,patch,YES);
    
    /* extrinsic curvature */
    add_field("K",0,patch,YES);
    
  }
}

/* taking partial derivatives of fields needed in equations */
void bbn_partial_derivatives_fields(Grid_T *const grid)
{
}

