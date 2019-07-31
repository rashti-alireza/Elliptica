/*
// Alireza Rashti
// July 2019
*/

#include "free_data.h"

/* populating the free data part of initial data that we chose ourself */
void populate_free_data(Grid_T *const grid)
{
  /* populate conformal metric and its inverse */
  _gammas(grid);
  
  /* Christoffer symbols made up of conformal metric */
  //_Gamma(grid);
  
  /* Ricci scalar made up of conformal metric _gamma */
  //_R(grid);
  
  /* extrinsic curvature */
  //K(grid);
}

/* to make christoffer symbol we need derivative of the metric, 
// this function does that. */
static void preparing_conformal_metric_derivatives(Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    /* declaring conformal metric */
    DECLARE_FIELD(_gamma_D2D2)
    DECLARE_FIELD(_gamma_D0D2)
    DECLARE_FIELD(_gamma_D0D0)
    DECLARE_FIELD(_gamma_D0D1)
    DECLARE_FIELD(_gamma_D1D2)
    DECLARE_FIELD(_gamma_D1D1)

    /* add fields for derivative of the conformal metric */
    ADD_FIELD(_dgamma_D0D0D1)
    ADD_FIELD(_dgamma_D0D0D0)
    ADD_FIELD(_dgamma_D2D2D2)
    ADD_FIELD(_dgamma_D0D0D2)
    ADD_FIELD(_dgamma_D0D2D1)
    ADD_FIELD(_dgamma_D1D1D0)
    ADD_FIELD(_dgamma_D1D1D2)
    ADD_FIELD(_dgamma_D2D2D0)
    ADD_FIELD(_dgamma_D2D2D1)
    ADD_FIELD(_dgamma_D0D1D0)
    ADD_FIELD(_dgamma_D0D1D1)
    ADD_FIELD(_dgamma_D0D1D2)
    ADD_FIELD(_dgamma_D0D2D0)
    ADD_FIELD(_dgamma_D1D1D1)
    ADD_FIELD(_dgamma_D1D2D1)
    ADD_FIELD(_dgamma_D1D2D2)
    ADD_FIELD(_dgamma_D1D2D0)
    ADD_FIELD(_dgamma_D0D2D2)

    DECLARE_FIELD(_dgamma_D0D0D1)
    DECLARE_FIELD(_dgamma_D0D0D0)
    DECLARE_FIELD(_dgamma_D2D2D2)
    DECLARE_FIELD(_dgamma_D0D0D2)
    DECLARE_FIELD(_dgamma_D0D2D1)
    DECLARE_FIELD(_dgamma_D1D1D0)
    DECLARE_FIELD(_dgamma_D1D1D2)
    DECLARE_FIELD(_dgamma_D2D2D0)
    DECLARE_FIELD(_dgamma_D2D2D1)
    DECLARE_FIELD(_dgamma_D0D1D0)
    DECLARE_FIELD(_dgamma_D0D1D1)
    DECLARE_FIELD(_dgamma_D0D1D2)
    DECLARE_FIELD(_dgamma_D0D2D0)
    DECLARE_FIELD(_dgamma_D1D1D1)
    DECLARE_FIELD(_dgamma_D1D2D1)
    DECLARE_FIELD(_dgamma_D1D2D2)
    DECLARE_FIELD(_dgamma_D1D2D0)
    DECLARE_FIELD(_dgamma_D0D2D2)

    /* filling the values */
    _dgamma_D0D0D0->v = Partial_Derivative(_gamma_D0D0,"x");
    _dgamma_D0D0D1->v = Partial_Derivative(_gamma_D0D0,"y");
    _dgamma_D0D0D2->v = Partial_Derivative(_gamma_D0D0,"z");
    _dgamma_D1D2D2->v = Partial_Derivative(_gamma_D1D2,"z");
    _dgamma_D1D1D0->v = Partial_Derivative(_gamma_D1D1,"x");
    _dgamma_D0D1D1->v = Partial_Derivative(_gamma_D0D1,"y");
    _dgamma_D0D1D0->v = Partial_Derivative(_gamma_D0D1,"x");
    _dgamma_D1D1D2->v = Partial_Derivative(_gamma_D1D1,"z");
    _dgamma_D0D1D2->v = Partial_Derivative(_gamma_D0D1,"z");
    _dgamma_D0D2D2->v = Partial_Derivative(_gamma_D0D2,"z");
    _dgamma_D1D2D1->v = Partial_Derivative(_gamma_D1D2,"y");
    _dgamma_D0D2D0->v = Partial_Derivative(_gamma_D0D2,"x");
    _dgamma_D0D2D1->v = Partial_Derivative(_gamma_D0D2,"y");
    _dgamma_D2D2D2->v = Partial_Derivative(_gamma_D2D2,"z");
    _dgamma_D2D2D0->v = Partial_Derivative(_gamma_D2D2,"x");
    _dgamma_D1D2D0->v = Partial_Derivative(_gamma_D1D2,"x");
    _dgamma_D1D1D1->v = Partial_Derivative(_gamma_D1D1,"y");
    _dgamma_D2D2D1->v = Partial_Derivative(_gamma_D2D2,"y");

  }
  
}

/* freeing conformal metric derivatives */
static void free_conformal_metric_derivatives(Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    /* declare _dgamma */
    DECLARE_FIELD(_dgamma_D0D0D1)
    DECLARE_FIELD(_dgamma_D0D0D0)
    DECLARE_FIELD(_dgamma_D2D2D2)
    DECLARE_FIELD(_dgamma_D0D0D2)
    DECLARE_FIELD(_dgamma_D0D2D1)
    DECLARE_FIELD(_dgamma_D1D1D0)
    DECLARE_FIELD(_dgamma_D1D1D2)
    DECLARE_FIELD(_dgamma_D2D2D0)
    DECLARE_FIELD(_dgamma_D2D2D1)
    DECLARE_FIELD(_dgamma_D0D1D0)
    DECLARE_FIELD(_dgamma_D0D1D1)
    DECLARE_FIELD(_dgamma_D0D1D2)
    DECLARE_FIELD(_dgamma_D0D2D0)
    DECLARE_FIELD(_dgamma_D1D1D1)
    DECLARE_FIELD(_dgamma_D1D2D1)
    DECLARE_FIELD(_dgamma_D1D2D2)
    DECLARE_FIELD(_dgamma_D1D2D0)
    DECLARE_FIELD(_dgamma_D0D2D2)

    /* removing fields */
    REMOVE_FIELD(_dgamma_D0D0D1)
    REMOVE_FIELD(_dgamma_D0D0D0)
    REMOVE_FIELD(_dgamma_D2D2D2)
    REMOVE_FIELD(_dgamma_D0D0D2)
    REMOVE_FIELD(_dgamma_D0D2D1)
    REMOVE_FIELD(_dgamma_D1D1D0)
    REMOVE_FIELD(_dgamma_D1D1D2)
    REMOVE_FIELD(_dgamma_D2D2D0)
    REMOVE_FIELD(_dgamma_D2D2D1)
    REMOVE_FIELD(_dgamma_D0D1D0)
    REMOVE_FIELD(_dgamma_D0D1D1)
    REMOVE_FIELD(_dgamma_D0D1D2)
    REMOVE_FIELD(_dgamma_D0D2D0)
    REMOVE_FIELD(_dgamma_D1D1D1)
    REMOVE_FIELD(_dgamma_D1D2D1)
    REMOVE_FIELD(_dgamma_D1D2D2)
    REMOVE_FIELD(_dgamma_D1D2D0)
    REMOVE_FIELD(_dgamma_D0D2D2)

  }
}

/* populate conformal metric and its inverse */
static void _gammas(Grid_T *const grid)
{
  /* roll off rate at exp(-(r/r0)^4)  */
  const double r0   = 0.5*GetParameterD_E("BH_NS_separation");
  const double M_BH = GetParameterD_E("BH_mass");
  const double a    = GetParameterD_E("BH_dimensionless_spin")*M_BH;
  const double a2   = SQR(a);
  double H,k0,k1,k2;/* in ds^2 = (eta_ij+2*H*ki*kj)dx^i*dx^j */
  /* center of BH */
  const double C_BH = 0.5*GetParameterD_E("BH_NS_separation");
  unsigned p,ijk,nn;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    nn = patch->nn;
    
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
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x   = patch->node[ijk]->x[0];
      double y   = patch->node[ijk]->x[1]-C_BH;
      double z   = patch->node[ijk]->x[2];
      double r2 = SQR(x)+SQR(y)+SQR(z);
      double r  = sqrt(r2);
      double rbar2  = 0.5*(r2-a2+sqrt(SQR(r2-a2)+4*a2*SQR(z)));
      double rbar   = sqrt(rbar2);
      double e   = exp(-pow(r/r0,4));
      
      k0 = (rbar*x+a*y)/(rbar2+a2);
      k1 = (rbar*y-a*x)/(rbar2+a2);
      k2 = z/rbar;
      H  = M_BH*rbar/(rbar2+a2*SQR(k2));
      double C = 2.*H*e;
      double A = 1./(1+C*(SQR(k0)+SQR(k1)+SQR(k2)));
      
      _gamma_D0D0[ijk] = 1.+C*k0*k0;
      _gamma_D0D1[ijk] = C*k0*k1;
      _gamma_D0D2[ijk] = C*k0*k2;
      _gamma_D1D1[ijk] = 1+C*k1*k1;
      _gamma_D1D2[ijk] = C*k1*k2;
      _gamma_D2D2[ijk] = 1+C*k2*k2;
      
      _gammaI_U0U0[ijk] = A*(1+C*(SQR(k1)+SQR(k2)));
      _gammaI_U0U1[ijk] = -A*C*k0*k1;
      _gammaI_U0U2[ijk] = -A*C*k0*k2;
      _gammaI_U1U1[ijk] = A*(1+C*(SQR(k0)+SQR(k2)));
      _gammaI_U1U2[ijk] = -A*C*k1*k2;
      _gammaI_U2U2[ijk] = A*(1+C*(SQR(k0)+SQR(k1)));
      
      /* quick test check _gamma * _gammaI = delta */
      if (1)
      {
          double delta_U0D0 = 
        _gammaI_U0U0[ijk]*_gamma_D0D0[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D0D2[ijk];

          double delta_U0D1 = 
        _gammaI_U0U0[ijk]*_gamma_D0D1[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D1D2[ijk];

          double delta_U0D2 = 
        _gammaI_U0U0[ijk]*_gamma_D0D2[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U0U2[ijk]*_gamma_D2D2[ijk];

          double delta_U1D2 = 
        _gammaI_U0U1[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U1U2[ijk]*_gamma_D2D2[ijk];

          double delta_U1D0 = 
        _gammaI_U0U1[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D0D2[ijk];

         double delta_U1D1 = 
        _gammaI_U0U1[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D1D2[ijk];

          double delta_U2D2 = 
        _gammaI_U0U2[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U2U2[ijk]*_gamma_D2D2[ijk];

          double delta_U2D0 = 
        _gammaI_U0U2[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D0D2[ijk];

          double delta_U2D1 = 
        _gammaI_U0U2[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D1D2[ijk];

        if(!EQL(delta_U1D1,1))  abortEr("_gammaI is not correct!\n");
        if(!EQL(delta_U0D1,0))  abortEr("_gammaI is not correct!\n");
        if(!EQL(delta_U0D2,0))  abortEr("_gammaI is not correct!\n");
        if(!EQL(delta_U1D2,0))  abortEr("_gammaI is not correct!\n");
        if(!EQL(delta_U0D0,1))  abortEr("_gammaI is not correct!\n");
        if(!EQL(delta_U2D1,0))  abortEr("_gammaI is not correct!\n");
        if(!EQL(delta_U2D2,1))  abortEr("_gammaI is not correct!\n");
        if(!EQL(delta_U2D0,0))  abortEr("_gammaI is not correct!\n");
        if(!EQL(delta_U1D0,0))  abortEr("_gammaI is not correct!\n");

      }
      
    }
  }
}

