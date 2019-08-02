/*
// Alireza Rashti
// July 2019
*/

#include "free_data.h"

/* populating the free data part of initial data that we chose ourself */
void populate_free_data(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Populating free data and related ...\n");
  
  /* populate conformal metric and its inverse */
  _gammas(grid);
  printf("Conformal metric and its inverse ~> Done.\n");
  
  /* Christoffer symbols made up of conformal metric */
  _Gamma(grid);
  printf("Christoffer symbols (_Gamma)     ~> Done.\n");
  
  /* partial derivtive of _Gamma, used in covariant derivative and _R */
  _dGamma(grid);
  printf("Partial derivatives of _Gamma    ~> Done.\n");
  
  /* Ricci scalar made up of conformal metric _gamma */
  _Ricci(grid);
  printf("Ricci scalar (_R)                ~> Done.\n");
  
  /* trace of Kerr Schild extrinsic curvature */
  tr_KSKij(grid);
  printf("Trace of Kerr Schild BH:tr(K_ij) ~> Done.\n");
  
  printf("Populating free data and related ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* partial derivtive of _Gamma, used in covariant derivative and _R */
static void _dGamma(Grid_T *const grid)
{
  const unsigned np = grid->np;
  unsigned p;

  OpenMP_1d_Pragma(omp parallel for)
  for(p = 0; p < np; ++p)
  {
    Patch_T *patch = grid->patch[p];
  
    /* partial derivative of _Gamma */
    DECLARE_FIELD(_dGamma_U2D2D2D2)
    DECLARE_FIELD(_dGamma_U2D2D2D0)
    DECLARE_FIELD(_dGamma_U2D2D2D1)
    DECLARE_FIELD(_dGamma_U2D0D0D2)
    DECLARE_FIELD(_dGamma_U1D1D2D2)
    DECLARE_FIELD(_dGamma_U2D0D0D0)
    DECLARE_FIELD(_dGamma_U1D1D2D0)
    DECLARE_FIELD(_dGamma_U1D1D2D1)
    DECLARE_FIELD(_dGamma_U2D1D1D0)
    DECLARE_FIELD(_dGamma_U2D0D0D1)
    DECLARE_FIELD(_dGamma_U2D0D2D1)
    DECLARE_FIELD(_dGamma_U1D0D1D0)
    DECLARE_FIELD(_dGamma_U1D0D1D1)
    DECLARE_FIELD(_dGamma_U1D0D1D2)
    DECLARE_FIELD(_dGamma_U1D2D2D1)
    DECLARE_FIELD(_dGamma_U1D0D0D1)
    DECLARE_FIELD(_dGamma_U1D0D0D0)
    DECLARE_FIELD(_dGamma_U1D0D0D2)
    DECLARE_FIELD(_dGamma_U0D1D2D2)
    DECLARE_FIELD(_dGamma_U0D1D2D1)
    DECLARE_FIELD(_dGamma_U0D1D2D0)
    DECLARE_FIELD(_dGamma_U2D0D2D0)
    DECLARE_FIELD(_dGamma_U1D0D2D2)
    DECLARE_FIELD(_dGamma_U1D0D2D1)
    DECLARE_FIELD(_dGamma_U1D0D2D0)
    DECLARE_FIELD(_dGamma_U2D1D1D2)
    DECLARE_FIELD(_dGamma_U2D0D2D2)
    DECLARE_FIELD(_dGamma_U0D0D1D0)
    DECLARE_FIELD(_dGamma_U1D2D2D0)
    DECLARE_FIELD(_dGamma_U2D1D2D1)
    DECLARE_FIELD(_dGamma_U2D0D1D2)
    DECLARE_FIELD(_dGamma_U2D0D1D1)
    DECLARE_FIELD(_dGamma_U2D0D1D0)
    DECLARE_FIELD(_dGamma_U2D1D2D2)
    DECLARE_FIELD(_dGamma_U0D1D1D0)
    DECLARE_FIELD(_dGamma_U0D1D1D1)
    DECLARE_FIELD(_dGamma_U0D1D1D2)
    DECLARE_FIELD(_dGamma_U1D2D2D2)
    DECLARE_FIELD(_dGamma_U1D1D1D1)
    DECLARE_FIELD(_dGamma_U1D1D1D0)
    DECLARE_FIELD(_dGamma_U1D1D1D2)
    DECLARE_FIELD(_dGamma_U0D0D1D1)
    DECLARE_FIELD(_dGamma_U0D0D2D2)
    DECLARE_FIELD(_dGamma_U0D0D2D0)
    DECLARE_FIELD(_dGamma_U0D0D2D1)
    DECLARE_FIELD(_dGamma_U2D1D2D0)
    DECLARE_FIELD(_dGamma_U0D0D0D0)
    DECLARE_FIELD(_dGamma_U0D0D0D1)
    DECLARE_FIELD(_dGamma_U0D0D0D2)
    DECLARE_FIELD(_dGamma_U2D1D1D1)
    DECLARE_FIELD(_dGamma_U0D2D2D0)
    DECLARE_FIELD(_dGamma_U0D2D2D1)
    DECLARE_FIELD(_dGamma_U0D2D2D2)
    DECLARE_FIELD(_dGamma_U0D0D1D2)
    
    /* _Gamma */
    DECLARE_FIELD(_Gamma_U2D1D1)
    DECLARE_FIELD(_Gamma_U2D1D2)
    DECLARE_FIELD(_Gamma_U0D1D1)
    DECLARE_FIELD(_Gamma_U2D0D2)
    DECLARE_FIELD(_Gamma_U2D2D2)
    DECLARE_FIELD(_Gamma_U0D1D2)
    DECLARE_FIELD(_Gamma_U0D0D2)
    DECLARE_FIELD(_Gamma_U0D0D1)
    DECLARE_FIELD(_Gamma_U0D0D0)
    DECLARE_FIELD(_Gamma_U1D2D2)
    DECLARE_FIELD(_Gamma_U2D0D1)
    DECLARE_FIELD(_Gamma_U0D2D2)
    DECLARE_FIELD(_Gamma_U2D0D0)
    DECLARE_FIELD(_Gamma_U1D0D2)
    DECLARE_FIELD(_Gamma_U1D1D2)
    DECLARE_FIELD(_Gamma_U1D0D0)
    DECLARE_FIELD(_Gamma_U1D0D1)
    DECLARE_FIELD(_Gamma_U1D1D1)

    /* populating partial derivative of _Gamma */
    _dGamma_U2D2D2D2->v = Partial_Derivative(_Gamma_U2D2D2,"z");
    _dGamma_U2D2D2D0->v = Partial_Derivative(_Gamma_U2D2D2,"x");
    _dGamma_U2D2D2D1->v = Partial_Derivative(_Gamma_U2D2D2,"y");
    _dGamma_U2D0D0D2->v = Partial_Derivative(_Gamma_U2D0D0,"z");
    _dGamma_U1D1D2D2->v = Partial_Derivative(_Gamma_U1D1D2,"z");
    _dGamma_U2D0D0D0->v = Partial_Derivative(_Gamma_U2D0D0,"x");
    _dGamma_U1D1D2D0->v = Partial_Derivative(_Gamma_U1D1D2,"x");
    _dGamma_U1D1D2D1->v = Partial_Derivative(_Gamma_U1D1D2,"y");
    _dGamma_U2D1D1D0->v = Partial_Derivative(_Gamma_U2D1D1,"x");
    _dGamma_U2D0D0D1->v = Partial_Derivative(_Gamma_U2D0D0,"y");
    _dGamma_U2D0D2D1->v = Partial_Derivative(_Gamma_U2D0D2,"y");
    _dGamma_U1D0D1D0->v = Partial_Derivative(_Gamma_U1D0D1,"x");
    _dGamma_U1D0D1D1->v = Partial_Derivative(_Gamma_U1D0D1,"y");
    _dGamma_U1D0D1D2->v = Partial_Derivative(_Gamma_U1D0D1,"z");
    _dGamma_U1D2D2D1->v = Partial_Derivative(_Gamma_U1D2D2,"y");
    _dGamma_U1D0D0D1->v = Partial_Derivative(_Gamma_U1D0D0,"y");
    _dGamma_U1D0D0D0->v = Partial_Derivative(_Gamma_U1D0D0,"x");
    _dGamma_U1D0D0D2->v = Partial_Derivative(_Gamma_U1D0D0,"z");
    _dGamma_U0D1D2D2->v = Partial_Derivative(_Gamma_U0D1D2,"z");
    _dGamma_U0D1D2D1->v = Partial_Derivative(_Gamma_U0D1D2,"y");
    _dGamma_U0D1D2D0->v = Partial_Derivative(_Gamma_U0D1D2,"x");
    _dGamma_U2D0D2D0->v = Partial_Derivative(_Gamma_U2D0D2,"x");
    _dGamma_U1D0D2D2->v = Partial_Derivative(_Gamma_U1D0D2,"z");
    _dGamma_U1D0D2D1->v = Partial_Derivative(_Gamma_U1D0D2,"y");
    _dGamma_U1D0D2D0->v = Partial_Derivative(_Gamma_U1D0D2,"x");
    _dGamma_U2D1D1D2->v = Partial_Derivative(_Gamma_U2D1D1,"z");
    _dGamma_U2D0D2D2->v = Partial_Derivative(_Gamma_U2D0D2,"z");
    _dGamma_U0D0D1D0->v = Partial_Derivative(_Gamma_U0D0D1,"x");
    _dGamma_U1D2D2D0->v = Partial_Derivative(_Gamma_U1D2D2,"x");
    _dGamma_U2D1D2D1->v = Partial_Derivative(_Gamma_U2D1D2,"y");
    _dGamma_U2D0D1D2->v = Partial_Derivative(_Gamma_U2D0D1,"z");
    _dGamma_U2D0D1D1->v = Partial_Derivative(_Gamma_U2D0D1,"y");
    _dGamma_U2D0D1D0->v = Partial_Derivative(_Gamma_U2D0D1,"x");
    _dGamma_U2D1D2D2->v = Partial_Derivative(_Gamma_U2D1D2,"z");
    _dGamma_U0D1D1D0->v = Partial_Derivative(_Gamma_U0D1D1,"x");
    _dGamma_U0D1D1D1->v = Partial_Derivative(_Gamma_U0D1D1,"y");
    _dGamma_U0D1D1D2->v = Partial_Derivative(_Gamma_U0D1D1,"z");
    _dGamma_U1D2D2D2->v = Partial_Derivative(_Gamma_U1D2D2,"z");
    _dGamma_U1D1D1D1->v = Partial_Derivative(_Gamma_U1D1D1,"y");
    _dGamma_U1D1D1D0->v = Partial_Derivative(_Gamma_U1D1D1,"x");
    _dGamma_U1D1D1D2->v = Partial_Derivative(_Gamma_U1D1D1,"z");
    _dGamma_U0D0D1D1->v = Partial_Derivative(_Gamma_U0D0D1,"y");
    _dGamma_U0D0D2D2->v = Partial_Derivative(_Gamma_U0D0D2,"z");
    _dGamma_U0D0D2D0->v = Partial_Derivative(_Gamma_U0D0D2,"x");
    _dGamma_U0D0D2D1->v = Partial_Derivative(_Gamma_U0D0D2,"y");
    _dGamma_U2D1D2D0->v = Partial_Derivative(_Gamma_U2D1D2,"x");
    _dGamma_U0D0D0D0->v = Partial_Derivative(_Gamma_U0D0D0,"x");
    _dGamma_U0D0D0D1->v = Partial_Derivative(_Gamma_U0D0D0,"y");
    _dGamma_U0D0D0D2->v = Partial_Derivative(_Gamma_U0D0D0,"z");
    _dGamma_U2D1D1D1->v = Partial_Derivative(_Gamma_U2D1D1,"y");
    _dGamma_U0D2D2D0->v = Partial_Derivative(_Gamma_U0D2D2,"x");
    _dGamma_U0D2D2D1->v = Partial_Derivative(_Gamma_U0D2D2,"y");
    _dGamma_U0D2D2D2->v = Partial_Derivative(_Gamma_U0D2D2,"z");
    _dGamma_U0D0D1D2->v = Partial_Derivative(_Gamma_U0D0D1,"z");

  }
}

/* to make christoffer symbol we need derivative of the metric, 
// this function does that. */
static void preparing_conformal_metric_derivatives(Patch_T *const patch)
{
  /* declaring conformal metric */
  DECLARE_FIELD(_gamma_D2D2)
  DECLARE_FIELD(_gamma_D0D2)
  DECLARE_FIELD(_gamma_D0D0)
  DECLARE_FIELD(_gamma_D0D1)
  DECLARE_FIELD(_gamma_D1D2)
  DECLARE_FIELD(_gamma_D1D1)

  /* add fields for derivative of the conformal metric */
  ADD_FIELD_NoMem(_dgamma_D0D0D1)
  ADD_FIELD_NoMem(_dgamma_D0D0D0)
  ADD_FIELD_NoMem(_dgamma_D2D2D2)
  ADD_FIELD_NoMem(_dgamma_D0D0D2)
  ADD_FIELD_NoMem(_dgamma_D0D2D1)
  ADD_FIELD_NoMem(_dgamma_D1D1D0)
  ADD_FIELD_NoMem(_dgamma_D1D1D2)
  ADD_FIELD_NoMem(_dgamma_D2D2D0)
  ADD_FIELD_NoMem(_dgamma_D2D2D1)
  ADD_FIELD_NoMem(_dgamma_D0D1D0)
  ADD_FIELD_NoMem(_dgamma_D0D1D1)
  ADD_FIELD_NoMem(_dgamma_D0D1D2)
  ADD_FIELD_NoMem(_dgamma_D0D2D0)
  ADD_FIELD_NoMem(_dgamma_D1D1D1)
  ADD_FIELD_NoMem(_dgamma_D1D2D1)
  ADD_FIELD_NoMem(_dgamma_D1D2D2)
  ADD_FIELD_NoMem(_dgamma_D1D2D0)
  ADD_FIELD_NoMem(_dgamma_D0D2D2)

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

/* freeing conformal metric derivatives */
static void free_conformal_metric_derivatives(Patch_T *const patch)
{
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

/* populate conformal metric and its inverse */
static void _gammas(Grid_T *const grid)
{
  /* roll off distance at exp(-(r/r0)^4)  */
  const double r0   = GetParameterD_E("RollOff_distance");
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
      if (0)
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

/* Christoffer symbol */
static void _Gamma(Grid_T *const grid)
{
  const unsigned np = grid->np;
  unsigned p;

  OpenMP_1d_Pragma(omp parallel for)
  for(p = 0; p < np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    preparing_conformal_metric_derivatives(patch);

    /* declaring: */
    GET_FIELD(_gammaI_U0U2)
    GET_FIELD(_gammaI_U0U0)
    GET_FIELD(_gammaI_U0U1)
    GET_FIELD(_gammaI_U1U2)
    GET_FIELD(_gammaI_U1U1)
    GET_FIELD(_gammaI_U2U2)
    GET_FIELD(_dgamma_D0D0D1)
    GET_FIELD(_dgamma_D0D0D0)
    GET_FIELD(_dgamma_D2D2D2)
    GET_FIELD(_dgamma_D0D0D2)
    GET_FIELD(_dgamma_D0D2D1)
    GET_FIELD(_dgamma_D1D1D0)
    GET_FIELD(_dgamma_D1D1D2)
    GET_FIELD(_dgamma_D2D2D0)
    GET_FIELD(_dgamma_D2D2D1)
    GET_FIELD(_dgamma_D0D1D0)
    GET_FIELD(_dgamma_D0D1D1)
    GET_FIELD(_dgamma_D0D1D2)
    GET_FIELD(_dgamma_D0D2D0)
    GET_FIELD(_dgamma_D1D1D1)
    GET_FIELD(_dgamma_D1D2D1)
    GET_FIELD(_dgamma_D1D2D2)
    GET_FIELD(_dgamma_D1D2D0)
    GET_FIELD(_dgamma_D0D2D2)
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
    double GAMMA_U1D0D0 = 
0.5*_dgamma_D0D0D0[ijk]*_gammaI_U0U1[ijk] - 0.5*_gammaI_U1U1[ijk]*
(_dgamma_D0D0D1[ijk] - 2*_dgamma_D0D1D0[ijk]) - 0.5*_gammaI_U1U2[ijk]*
(_dgamma_D0D0D2[ijk] - 2*_dgamma_D0D2D0[ijk]);

    double GAMMA_U0D2D2 = 
0.5*_dgamma_D2D2D2[ijk]*_gammaI_U0U2[ijk] + 0.5*_gammaI_U0U0[ijk]*(2*
_dgamma_D0D2D2[ijk] - _dgamma_D2D2D0[ijk]) + 0.5*_gammaI_U0U1[ijk]*(2*
_dgamma_D1D2D2[ijk] - _dgamma_D2D2D1[ijk]);

    double GAMMA_U1D0D2 = 
0.5*_dgamma_D0D0D2[ijk]*_gammaI_U0U1[ijk] + 0.5*_dgamma_D2D2D0[ijk]*
_gammaI_U1U2[ijk] + 0.5*_gammaI_U1U1[ijk]*(_dgamma_D0D1D2[ijk] - 
_dgamma_D0D2D1[ijk] + _dgamma_D1D2D0[ijk]);

    double GAMMA_U2D2D2 = 
0.5*_dgamma_D2D2D2[ijk]*_gammaI_U2U2[ijk] + 0.5*_gammaI_U0U2[ijk]*(2*
_dgamma_D0D2D2[ijk] - _dgamma_D2D2D0[ijk]) + 0.5*_gammaI_U1U2[ijk]*(2*
_dgamma_D1D2D2[ijk] - _dgamma_D2D2D1[ijk]);

    double GAMMA_U2D1D1 = 
0.5*_dgamma_D1D1D1[ijk]*_gammaI_U1U2[ijk] + 0.5*_gammaI_U0U2[ijk]*(2*
_dgamma_D0D1D1[ijk] - _dgamma_D1D1D0[ijk]) - 0.5*_gammaI_U2U2[ijk]*
(_dgamma_D1D1D2[ijk] - 2*_dgamma_D1D2D1[ijk]);

    double GAMMA_U1D0D1 = 
0.5*_dgamma_D0D0D1[ijk]*_gammaI_U0U1[ijk] + 0.5*_dgamma_D1D1D0[ijk]*
_gammaI_U1U1[ijk] + 0.5*_gammaI_U1U2[ijk]*(-_dgamma_D0D1D2[ijk] + 
_dgamma_D0D2D1[ijk] + _dgamma_D1D2D0[ijk]);

    double GAMMA_U1D2D2 = 
0.5*_dgamma_D2D2D2[ijk]*_gammaI_U1U2[ijk] + 0.5*_gammaI_U0U1[ijk]*(2*
_dgamma_D0D2D2[ijk] - _dgamma_D2D2D0[ijk]) + 0.5*_gammaI_U1U1[ijk]*(2*
_dgamma_D1D2D2[ijk] - _dgamma_D2D2D1[ijk]);

    double GAMMA_U2D1D2 = 
0.5*_dgamma_D1D1D2[ijk]*_gammaI_U1U2[ijk] + 0.5*_dgamma_D2D2D1[ijk]*
_gammaI_U2U2[ijk] + 0.5*_gammaI_U0U2[ijk]*(_dgamma_D0D1D2[ijk] + 
_dgamma_D0D2D1[ijk] - _dgamma_D1D2D0[ijk]);

    double GAMMA_U2D0D2 = 
0.5*_dgamma_D0D0D2[ijk]*_gammaI_U0U2[ijk] + 0.5*_dgamma_D2D2D0[ijk]*
_gammaI_U2U2[ijk] + 0.5*_gammaI_U1U2[ijk]*(_dgamma_D0D1D2[ijk] - 
_dgamma_D0D2D1[ijk] + _dgamma_D1D2D0[ijk]);

    double GAMMA_U2D0D1 = 
0.5*_dgamma_D0D0D1[ijk]*_gammaI_U0U2[ijk] + 0.5*_dgamma_D1D1D0[ijk]*
_gammaI_U1U2[ijk] + 0.5*_gammaI_U2U2[ijk]*(-_dgamma_D0D1D2[ijk] + 
_dgamma_D0D2D1[ijk] + _dgamma_D1D2D0[ijk]);

    double GAMMA_U2D0D0 = 
0.5*_dgamma_D0D0D0[ijk]*_gammaI_U0U2[ijk] - 0.5*_gammaI_U1U2[ijk]*
(_dgamma_D0D0D1[ijk] - 2*_dgamma_D0D1D0[ijk]) - 0.5*_gammaI_U2U2[ijk]*
(_dgamma_D0D0D2[ijk] - 2*_dgamma_D0D2D0[ijk]);

    double GAMMA_U0D0D2 = 
0.5*_dgamma_D0D0D2[ijk]*_gammaI_U0U0[ijk] + 0.5*_dgamma_D2D2D0[ijk]*
_gammaI_U0U2[ijk] + 0.5*_gammaI_U0U1[ijk]*(_dgamma_D0D1D2[ijk] - 
_dgamma_D0D2D1[ijk] + _dgamma_D1D2D0[ijk]);

    double GAMMA_U1D1D1 = 
0.5*_dgamma_D1D1D1[ijk]*_gammaI_U1U1[ijk] + 0.5*_gammaI_U0U1[ijk]*(2*
_dgamma_D0D1D1[ijk] - _dgamma_D1D1D0[ijk]) - 0.5*_gammaI_U1U2[ijk]*
(_dgamma_D1D1D2[ijk] - 2*_dgamma_D1D2D1[ijk]);

    double GAMMA_U0D1D1 = 
0.5*_dgamma_D1D1D1[ijk]*_gammaI_U0U1[ijk] + 0.5*_gammaI_U0U0[ijk]*(2*
_dgamma_D0D1D1[ijk] - _dgamma_D1D1D0[ijk]) - 0.5*_gammaI_U0U2[ijk]*
(_dgamma_D1D1D2[ijk] - 2*_dgamma_D1D2D1[ijk]);

    double GAMMA_U0D0D0 = 
0.5*_dgamma_D0D0D0[ijk]*_gammaI_U0U0[ijk] - 0.5*_gammaI_U0U1[ijk]*
(_dgamma_D0D0D1[ijk] - 2*_dgamma_D0D1D0[ijk]) - 0.5*_gammaI_U0U2[ijk]*
(_dgamma_D0D0D2[ijk] - 2*_dgamma_D0D2D0[ijk]);

    double GAMMA_U1D1D2 = 
0.5*_dgamma_D1D1D2[ijk]*_gammaI_U1U1[ijk] + 0.5*_dgamma_D2D2D1[ijk]*
_gammaI_U1U2[ijk] + 0.5*_gammaI_U0U1[ijk]*(_dgamma_D0D1D2[ijk] + 
_dgamma_D0D2D1[ijk] - _dgamma_D1D2D0[ijk]);

    double GAMMA_U0D0D1 = 
0.5*_dgamma_D0D0D1[ijk]*_gammaI_U0U0[ijk] + 0.5*_dgamma_D1D1D0[ijk]*
_gammaI_U0U1[ijk] + 0.5*_gammaI_U0U2[ijk]*(-_dgamma_D0D1D2[ijk] + 
_dgamma_D0D2D1[ijk] + _dgamma_D1D2D0[ijk]);

    double GAMMA_U0D1D2 = 
0.5*_dgamma_D1D1D2[ijk]*_gammaI_U0U1[ijk] + 0.5*_dgamma_D2D2D1[ijk]*
_gammaI_U0U2[ijk] + 0.5*_gammaI_U0U0[ijk]*(_dgamma_D0D1D2[ijk] + 
_dgamma_D0D2D1[ijk] - _dgamma_D1D2D0[ijk]);


    /* populating: */
    _Gamma_U2D1D1[ijk] = GAMMA_U2D1D1;
    _Gamma_U2D1D2[ijk] = GAMMA_U2D1D2;
    _Gamma_U0D1D1[ijk] = GAMMA_U0D1D1;
    _Gamma_U2D0D2[ijk] = GAMMA_U2D0D2;
    _Gamma_U2D2D2[ijk] = GAMMA_U2D2D2;
    _Gamma_U0D1D2[ijk] = GAMMA_U0D1D2;
    _Gamma_U0D0D2[ijk] = GAMMA_U0D0D2;
    _Gamma_U0D0D1[ijk] = GAMMA_U0D0D1;
    _Gamma_U0D0D0[ijk] = GAMMA_U0D0D0;
    _Gamma_U1D2D2[ijk] = GAMMA_U1D2D2;
    _Gamma_U2D0D1[ijk] = GAMMA_U2D0D1;
    _Gamma_U0D2D2[ijk] = GAMMA_U0D2D2;
    _Gamma_U2D0D0[ijk] = GAMMA_U2D0D0;
    _Gamma_U1D0D2[ijk] = GAMMA_U1D0D2;
    _Gamma_U1D1D2[ijk] = GAMMA_U1D1D2;
    _Gamma_U1D0D0[ijk] = GAMMA_U1D0D0;
    _Gamma_U1D0D1[ijk] = GAMMA_U1D0D1;
    _Gamma_U1D1D1[ijk] = GAMMA_U1D1D1;
    }/*end of for(ijk = 0; ijk < nn; ++ijk)*/
    free_conformal_metric_derivatives(patch);
  }
}

/* Ricci scalar */
static void _Ricci(Grid_T *const grid)
{
  const unsigned np = grid->np;
  unsigned p;

  OpenMP_1d_Pragma(omp parallel for)
  for(p = 0; p < np; ++p)
  {
    Patch_T *patch = grid->patch[p];

    /* declaring: */
    GET_FIELD(_gammaI_U0U2)
    GET_FIELD(_gammaI_U0U0)
    GET_FIELD(_gammaI_U0U1)
    GET_FIELD(_gammaI_U1U2)
    GET_FIELD(_gammaI_U1U1)
    GET_FIELD(_gammaI_U2U2)
    GET_FIELD(_dGamma_U2D2D2D0)
    GET_FIELD(_dGamma_U2D2D2D1)
    GET_FIELD(_dGamma_U2D0D0D2)
    GET_FIELD(_dGamma_U1D1D2D2)
    GET_FIELD(_dGamma_U1D1D2D0)
    GET_FIELD(_dGamma_U1D1D2D1)
    GET_FIELD(_dGamma_U2D0D2D1)
    GET_FIELD(_dGamma_U1D0D1D0)
    GET_FIELD(_dGamma_U1D0D1D1)
    GET_FIELD(_dGamma_U1D0D1D2)
    GET_FIELD(_dGamma_U1D2D2D1)
    GET_FIELD(_dGamma_U1D0D0D1)
    GET_FIELD(_dGamma_U0D1D2D0)
    GET_FIELD(_dGamma_U2D0D2D0)
    GET_FIELD(_dGamma_U1D0D2D1)
    GET_FIELD(_dGamma_U2D1D1D2)
    GET_FIELD(_dGamma_U2D0D2D2)
    GET_FIELD(_dGamma_U0D0D1D0)
    GET_FIELD(_dGamma_U2D1D2D1)
    GET_FIELD(_dGamma_U2D0D1D2)
    GET_FIELD(_dGamma_U2D1D2D2)
    GET_FIELD(_dGamma_U0D1D1D0)
    GET_FIELD(_dGamma_U1D1D1D0)
    GET_FIELD(_dGamma_U1D1D1D2)
    GET_FIELD(_dGamma_U0D0D1D1)
    GET_FIELD(_dGamma_U0D0D2D2)
    GET_FIELD(_dGamma_U0D0D2D0)
    GET_FIELD(_dGamma_U0D0D2D1)
    GET_FIELD(_dGamma_U2D1D2D0)
    GET_FIELD(_dGamma_U0D0D0D1)
    GET_FIELD(_dGamma_U0D0D0D2)
    GET_FIELD(_dGamma_U0D2D2D0)
    GET_FIELD(_dGamma_U0D0D1D2)
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
    GET_FIELD(_R)
    for(ijk = 0; ijk < nn; ++ijk)
    {
    double R = 
_gammaI_U0U0[ijk]*(_Gamma_U0D0D0[ijk]*_Gamma_U1D0D1[ijk] +
_Gamma_U0D0D0[ijk]*_Gamma_U2D0D2[ijk] - _Gamma_U0D0D1[ijk]*
_Gamma_U1D0D0[ijk] - _Gamma_U0D0D2[ijk]*_Gamma_U2D0D0[ijk] +
_Gamma_U1D0D0[ijk]*_Gamma_U1D1D1[ijk] + _Gamma_U1D0D0[ijk]*
_Gamma_U2D1D2[ijk] - pow(_Gamma_U1D0D1[ijk], 2) - 2.0*
_Gamma_U1D0D2[ijk]*_Gamma_U2D0D1[ijk] + _Gamma_U1D1D2[ijk]*
_Gamma_U2D0D0[ijk] + _Gamma_U2D0D0[ijk]*_Gamma_U2D2D2[ijk] -
pow(_Gamma_U2D0D2[ijk], 2) + _dGamma_U1D0D0D1[ijk] -
_dGamma_U1D0D1D0[ijk] + _dGamma_U2D0D0D2[ijk] - _dGamma_U2D0D2D0[ijk]) +
_gammaI_U0U1[ijk]*(_Gamma_U0D0D1[ijk]*_Gamma_U1D0D1[ijk] +
_Gamma_U0D0D1[ijk]*_Gamma_U2D0D2[ijk] - _Gamma_U0D1D1[ijk]*
_Gamma_U1D0D0[ijk] - _Gamma_U0D1D2[ijk]*_Gamma_U2D0D0[ijk] +
_Gamma_U1D0D1[ijk]*_Gamma_U2D1D2[ijk] - _Gamma_U1D0D2[ijk]*
_Gamma_U2D1D1[ijk] + _Gamma_U2D0D1[ijk]*_Gamma_U2D2D2[ijk] -
_Gamma_U2D0D2[ijk]*_Gamma_U2D1D2[ijk] - _dGamma_U0D0D0D1[ijk] +
_dGamma_U0D0D1D0[ijk] + _dGamma_U2D0D1D2[ijk] - _dGamma_U2D0D2D1[ijk]) +
_gammaI_U0U1[ijk]*(_Gamma_U0D0D1[ijk]*_Gamma_U1D0D1[ijk] +
_Gamma_U0D0D1[ijk]*_Gamma_U2D0D2[ijk] - _Gamma_U0D1D1[ijk]*
_Gamma_U1D0D0[ijk] - _Gamma_U0D1D2[ijk]*_Gamma_U2D0D0[ijk] +
_Gamma_U1D0D1[ijk]*_Gamma_U2D1D2[ijk] - _Gamma_U1D0D2[ijk]*
_Gamma_U2D1D1[ijk] + _Gamma_U2D0D1[ijk]*_Gamma_U2D2D2[ijk] -
_Gamma_U2D0D2[ijk]*_Gamma_U2D1D2[ijk] + _dGamma_U1D0D1D1[ijk] -
_dGamma_U1D1D1D0[ijk] + _dGamma_U2D0D1D2[ijk] - _dGamma_U2D1D2D0[ijk]) +
_gammaI_U0U2[ijk]*(_Gamma_U0D0D2[ijk]*_Gamma_U1D0D1[ijk] +
_Gamma_U0D0D2[ijk]*_Gamma_U2D0D2[ijk] - _Gamma_U0D1D2[ijk]*
_Gamma_U1D0D0[ijk] - _Gamma_U0D2D2[ijk]*_Gamma_U2D0D0[ijk] -
_Gamma_U1D0D1[ijk]*_Gamma_U1D1D2[ijk] + _Gamma_U1D0D2[ijk]*
_Gamma_U1D1D1[ijk] + _Gamma_U1D1D2[ijk]*_Gamma_U2D0D2[ijk] -
_Gamma_U1D2D2[ijk]*_Gamma_U2D0D1[ijk] - _dGamma_U0D0D0D2[ijk] +
_dGamma_U0D0D2D0[ijk] - _dGamma_U1D0D1D2[ijk] + _dGamma_U1D0D2D1[ijk]) +
_gammaI_U0U2[ijk]*(_Gamma_U0D0D2[ijk]*_Gamma_U1D0D1[ijk] +
_Gamma_U0D0D2[ijk]*_Gamma_U2D0D2[ijk] - _Gamma_U0D1D2[ijk]*
_Gamma_U1D0D0[ijk] - _Gamma_U0D2D2[ijk]*_Gamma_U2D0D0[ijk] -
_Gamma_U1D0D1[ijk]*_Gamma_U1D1D2[ijk] + _Gamma_U1D0D2[ijk]*
_Gamma_U1D1D1[ijk] + _Gamma_U1D1D2[ijk]*_Gamma_U2D0D2[ijk] -
_Gamma_U1D2D2[ijk]*_Gamma_U2D0D1[ijk] + _dGamma_U1D0D2D1[ijk] -
_dGamma_U1D1D2D0[ijk] + _dGamma_U2D0D2D2[ijk] - _dGamma_U2D2D2D0[ijk]) +
_gammaI_U1U1[ijk]*(_Gamma_U0D0D0[ijk]*_Gamma_U0D1D1[ijk] -
pow(_Gamma_U0D0D1[ijk], 2) + _Gamma_U0D0D1[ijk]*_Gamma_U1D1D1[ijk] +
_Gamma_U0D0D2[ijk]*_Gamma_U2D1D1[ijk] - _Gamma_U0D1D1[ijk]*
_Gamma_U1D0D1[ijk] + _Gamma_U0D1D1[ijk]*_Gamma_U2D0D2[ijk] - 2.0*
_Gamma_U0D1D2[ijk]*_Gamma_U2D0D1[ijk] + _Gamma_U1D1D1[ijk]*
_Gamma_U2D1D2[ijk] - _Gamma_U1D1D2[ijk]*_Gamma_U2D1D1[ijk] +
_Gamma_U2D1D1[ijk]*_Gamma_U2D2D2[ijk] - pow(_Gamma_U2D1D2[ijk], 2) -
_dGamma_U0D0D1D1[ijk] + _dGamma_U0D1D1D0[ijk] + _dGamma_U2D1D1D2[ijk] -
_dGamma_U2D1D2D1[ijk]) + _gammaI_U1U2[ijk]*(_Gamma_U0D0D0[ijk]*
_Gamma_U0D1D2[ijk] - _Gamma_U0D0D1[ijk]*_Gamma_U0D0D2[ijk] +
_Gamma_U0D0D1[ijk]*_Gamma_U1D1D2[ijk] + _Gamma_U0D0D2[ijk]*
_Gamma_U2D1D2[ijk] - _Gamma_U0D1D1[ijk]*_Gamma_U1D0D2[ijk] -
_Gamma_U0D2D2[ijk]*_Gamma_U2D0D1[ijk] + _Gamma_U1D1D2[ijk]*
_Gamma_U2D1D2[ijk] - _Gamma_U1D2D2[ijk]*_Gamma_U2D1D1[ijk] -
_dGamma_U0D0D1D2[ijk] + _dGamma_U0D1D2D0[ijk] - _dGamma_U1D1D1D2[ijk] +
_dGamma_U1D1D2D1[ijk]) + _gammaI_U1U2[ijk]*(_Gamma_U0D0D0[ijk]*
_Gamma_U0D1D2[ijk] - _Gamma_U0D0D1[ijk]*_Gamma_U0D0D2[ijk] +
_Gamma_U0D0D1[ijk]*_Gamma_U1D1D2[ijk] + _Gamma_U0D0D2[ijk]*
_Gamma_U2D1D2[ijk] - _Gamma_U0D1D1[ijk]*_Gamma_U1D0D2[ijk] -
_Gamma_U0D2D2[ijk]*_Gamma_U2D0D1[ijk] + _Gamma_U1D1D2[ijk]*
_Gamma_U2D1D2[ijk] - _Gamma_U1D2D2[ijk]*_Gamma_U2D1D1[ijk] -
_dGamma_U0D0D2D1[ijk] + _dGamma_U0D1D2D0[ijk] + _dGamma_U2D1D2D2[ijk] -
_dGamma_U2D2D2D1[ijk]) + _gammaI_U2U2[ijk]*(_Gamma_U0D0D0[ijk]*
_Gamma_U0D2D2[ijk] + _Gamma_U0D0D1[ijk]*_Gamma_U1D2D2[ijk] -
pow(_Gamma_U0D0D2[ijk], 2) + _Gamma_U0D0D2[ijk]*_Gamma_U2D2D2[ijk] -
2.0*_Gamma_U0D1D2[ijk]*_Gamma_U1D0D2[ijk] + _Gamma_U0D2D2[ijk]*
_Gamma_U1D0D1[ijk] - _Gamma_U0D2D2[ijk]*_Gamma_U2D0D2[ijk] +
_Gamma_U1D1D1[ijk]*_Gamma_U1D2D2[ijk] - pow(_Gamma_U1D1D2[ijk], 2) +
_Gamma_U1D1D2[ijk]*_Gamma_U2D2D2[ijk] - _Gamma_U1D2D2[ijk]*
_Gamma_U2D1D2[ijk] - _dGamma_U0D0D2D2[ijk] + _dGamma_U0D2D2D0[ijk] -
_dGamma_U1D1D2D2[ijk] + _dGamma_U1D2D2D1[ijk]);

      _R[ijk] = R;
    }/*end of for(ijk = 0; ijk < nn; ++ijk)*/
  }
}

/* trace of Kerr Schild extrinsic curvature */
static void tr_KSKij(Grid_T *const grid)
{
  const unsigned np = grid->np;
  unsigned p;

  OpenMP_1d_Pragma(omp parallel for)
  for(p = 0; p < np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    /* populating Kerr Schild gammas, lapse and shift: */
    populate_KSgammas_KSalpha_KSBeta(patch);

    /* populating Chirstoffer symbols */
    populating_KSGamma(patch);
    
    /* d of KS shift vector */
    partial_derivative_KSBeta(patch);

    /* populate KS trKij */
    populate_KS_trKij(patch);
    
    /* now free */
    free_KSfields(patch);
  }/* end of for(p = 0; p < np; ++p) */
  
}

/* free unwanted Kerr Schild fields made for trKij */
static void free_KSfields(Patch_T *const patch)
{
  
  DECLARE_FIELD(KSgamma_D1D1)
  DECLARE_FIELD(KSgamma_D1D2)
  DECLARE_FIELD(KSgamma_D0D0)
  DECLARE_FIELD(KSgamma_D0D1)
  DECLARE_FIELD(KSgamma_D0D2)
  DECLARE_FIELD(KSgamma_D2D2)
  DECLARE_FIELD(KSgammaI_U2U2)
  DECLARE_FIELD(KSgammaI_U1U1)
  DECLARE_FIELD(KSgammaI_U1U2)
  DECLARE_FIELD(KSgammaI_U0U0)
  DECLARE_FIELD(KSgammaI_U0U1)
  DECLARE_FIELD(KSgammaI_U0U2)
  DECLARE_FIELD(dKSgamma_D2D2D0)
  DECLARE_FIELD(dKSgamma_D1D2D2)
  DECLARE_FIELD(dKSgamma_D0D2D1)
  DECLARE_FIELD(dKSgamma_D0D1D2)
  DECLARE_FIELD(dKSgamma_D0D2D0)
  DECLARE_FIELD(dKSgamma_D0D1D0)
  DECLARE_FIELD(dKSgamma_D0D1D1)
  DECLARE_FIELD(dKSgamma_D1D1D1)
  DECLARE_FIELD(dKSgamma_D1D1D0)
  DECLARE_FIELD(dKSgamma_D1D2D0)
  DECLARE_FIELD(dKSgamma_D0D0D2)
  DECLARE_FIELD(dKSgamma_D0D0D1)
  DECLARE_FIELD(dKSgamma_D0D0D0)
  DECLARE_FIELD(dKSgamma_D0D2D2)
  DECLARE_FIELD(dKSgamma_D1D1D2)
  DECLARE_FIELD(dKSgamma_D2D2D2)
  DECLARE_FIELD(dKSgamma_D2D2D1)
  DECLARE_FIELD(dKSgamma_D1D2D1)
  DECLARE_FIELD(KSGamma_U1D0D0)
  DECLARE_FIELD(KSGamma_U0D2D2)
  DECLARE_FIELD(KSGamma_U1D0D2)
  DECLARE_FIELD(KSGamma_U1D0D1)
  DECLARE_FIELD(KSGamma_U2D1D2)
  DECLARE_FIELD(KSGamma_U0D1D1)
  DECLARE_FIELD(KSGamma_U2D0D0)
  DECLARE_FIELD(KSGamma_U1D1D2)
  DECLARE_FIELD(KSGamma_U2D0D1)
  DECLARE_FIELD(KSGamma_U2D0D2)
  DECLARE_FIELD(KSGamma_U0D0D1)
  DECLARE_FIELD(KSGamma_U0D0D0)
  DECLARE_FIELD(KSGamma_U0D1D2)
  DECLARE_FIELD(KSGamma_U0D0D2)
  DECLARE_FIELD(KSGamma_U2D2D2)
  DECLARE_FIELD(KSGamma_U1D1D1)
  DECLARE_FIELD(KSGamma_U2D1D1)
  DECLARE_FIELD(KSGamma_U1D2D2)
  DECLARE_FIELD(KSalpha)
  DECLARE_FIELD(KSB_D1)
  DECLARE_FIELD(KSB_D0)
  DECLARE_FIELD(KSB_D2)
  DECLARE_FIELD(dKSB_D2D0)
  DECLARE_FIELD(dKSB_D0D0)
  DECLARE_FIELD(dKSB_D2D2)
  DECLARE_FIELD(dKSB_D0D1)
  DECLARE_FIELD(dKSB_D0D2)
  DECLARE_FIELD(dKSB_D1D2)
  DECLARE_FIELD(dKSB_D1D1)
  DECLARE_FIELD(dKSB_D1D0)
  DECLARE_FIELD(dKSB_D2D1)
  
  REMOVE_FIELD(KSgamma_D1D1)
  REMOVE_FIELD(KSgamma_D1D2)
  REMOVE_FIELD(KSgamma_D0D0)
  REMOVE_FIELD(KSgamma_D0D1)
  REMOVE_FIELD(KSgamma_D0D2)
  REMOVE_FIELD(KSgamma_D2D2)
  REMOVE_FIELD(KSgammaI_U2U2)
  REMOVE_FIELD(KSgammaI_U1U1)
  REMOVE_FIELD(KSgammaI_U1U2)
  REMOVE_FIELD(KSgammaI_U0U0)
  REMOVE_FIELD(KSgammaI_U0U1)
  REMOVE_FIELD(KSgammaI_U0U2)
  REMOVE_FIELD(dKSgamma_D2D2D0)
  REMOVE_FIELD(dKSgamma_D1D2D2)
  REMOVE_FIELD(dKSgamma_D0D2D1)
  REMOVE_FIELD(dKSgamma_D0D1D2)
  REMOVE_FIELD(dKSgamma_D0D2D0)
  REMOVE_FIELD(dKSgamma_D0D1D0)
  REMOVE_FIELD(dKSgamma_D0D1D1)
  REMOVE_FIELD(dKSgamma_D1D1D1)
  REMOVE_FIELD(dKSgamma_D1D1D0)
  REMOVE_FIELD(dKSgamma_D1D2D0)
  REMOVE_FIELD(dKSgamma_D0D0D2)
  REMOVE_FIELD(dKSgamma_D0D0D1)
  REMOVE_FIELD(dKSgamma_D0D0D0)
  REMOVE_FIELD(dKSgamma_D0D2D2)
  REMOVE_FIELD(dKSgamma_D1D1D2)
  REMOVE_FIELD(dKSgamma_D2D2D2)
  REMOVE_FIELD(dKSgamma_D2D2D1)
  REMOVE_FIELD(dKSgamma_D1D2D1)
  REMOVE_FIELD(KSGamma_U1D0D0)
  REMOVE_FIELD(KSGamma_U0D2D2)
  REMOVE_FIELD(KSGamma_U1D0D2)
  REMOVE_FIELD(KSGamma_U1D0D1)
  REMOVE_FIELD(KSGamma_U2D1D2)
  REMOVE_FIELD(KSGamma_U0D1D1)
  REMOVE_FIELD(KSGamma_U2D0D0)
  REMOVE_FIELD(KSGamma_U1D1D2)
  REMOVE_FIELD(KSGamma_U2D0D1)
  REMOVE_FIELD(KSGamma_U2D0D2)
  REMOVE_FIELD(KSGamma_U0D0D1)
  REMOVE_FIELD(KSGamma_U0D0D0)
  REMOVE_FIELD(KSGamma_U0D1D2)
  REMOVE_FIELD(KSGamma_U0D0D2)
  REMOVE_FIELD(KSGamma_U2D2D2)
  REMOVE_FIELD(KSGamma_U1D1D1)
  REMOVE_FIELD(KSGamma_U2D1D1)
  REMOVE_FIELD(KSGamma_U1D2D2)
  REMOVE_FIELD(KSalpha)
  REMOVE_FIELD(KSB_D1)
  REMOVE_FIELD(KSB_D0)
  REMOVE_FIELD(KSB_D2)
  REMOVE_FIELD(dKSB_D2D0)
  REMOVE_FIELD(dKSB_D0D0)
  REMOVE_FIELD(dKSB_D2D2)
  REMOVE_FIELD(dKSB_D0D1)
  REMOVE_FIELD(dKSB_D0D2)
  REMOVE_FIELD(dKSB_D1D2)
  REMOVE_FIELD(dKSB_D1D1)
  REMOVE_FIELD(dKSB_D1D0)
  REMOVE_FIELD(dKSB_D2D1)

}

/* taking partial derivative of KSB */
static void partial_derivative_KSBeta(Patch_T *const patch)
{
  ADD_FIELD_NoMem(dKSB_D2D0)
  ADD_FIELD_NoMem(dKSB_D0D0)
  ADD_FIELD_NoMem(dKSB_D2D2)
  ADD_FIELD_NoMem(dKSB_D0D1)
  ADD_FIELD_NoMem(dKSB_D0D2)
  ADD_FIELD_NoMem(dKSB_D1D2)
  ADD_FIELD_NoMem(dKSB_D1D1)
  ADD_FIELD_NoMem(dKSB_D1D0)
  ADD_FIELD_NoMem(dKSB_D2D1)
  
  DECLARE_FIELD(dKSB_D2D0)
  DECLARE_FIELD(dKSB_D0D0)
  DECLARE_FIELD(dKSB_D2D2)
  DECLARE_FIELD(dKSB_D0D1)
  DECLARE_FIELD(dKSB_D0D2)
  DECLARE_FIELD(dKSB_D1D2)
  DECLARE_FIELD(dKSB_D1D1)
  DECLARE_FIELD(dKSB_D1D0)
  DECLARE_FIELD(dKSB_D2D1)
  DECLARE_FIELD(KSB_D1)
  DECLARE_FIELD(KSB_D0)
  DECLARE_FIELD(KSB_D2)

  /* Partial Derivative of KSB field */
  dKSB_D2D0->v = Partial_Derivative(KSB_D2,"x");
  dKSB_D0D0->v = Partial_Derivative(KSB_D0,"x");
  dKSB_D2D2->v = Partial_Derivative(KSB_D2,"z");
  dKSB_D0D1->v = Partial_Derivative(KSB_D0,"y");
  dKSB_D0D2->v = Partial_Derivative(KSB_D0,"z");
  dKSB_D1D2->v = Partial_Derivative(KSB_D1,"z");
  dKSB_D1D1->v = Partial_Derivative(KSB_D1,"y");
  dKSB_D1D0->v = Partial_Derivative(KSB_D1,"x");
  dKSB_D2D1->v = Partial_Derivative(KSB_D2,"y");
}

/* populating Kerr Schild gammas , lapse and shift vector */
static void populate_KSgammas_KSalpha_KSBeta(Patch_T *const patch)
{
  const double M_BH = GetParameterD_E("BH_mass");
  const double a    = GetParameterD_E("BH_dimensionless_spin")*M_BH;
  const double a2   = SQR(a);
  const double C_BH = 0.5*GetParameterD_E("BH_NS_separation");
  const unsigned nn = patch->nn;
  unsigned ijk;
  double H,k0,k1,k2;/* in ds^2 = (eta_ij+2*H*ki*kj)dx^i*dx^j */
      
  /* add Kerr Schild gammas */
  ADD_FIELD(KSgamma_D2D2)
  ADD_FIELD(KSgamma_D0D2)
  ADD_FIELD(KSgamma_D0D0)
  ADD_FIELD(KSgamma_D0D1)
  ADD_FIELD(KSgamma_D1D2)
  ADD_FIELD(KSgamma_D1D1)
  ADD_FIELD(KSgammaI_U0U2)
  ADD_FIELD(KSgammaI_U0U0)
  ADD_FIELD(KSgammaI_U0U1)
  ADD_FIELD(KSgammaI_U1U2)
  ADD_FIELD(KSgammaI_U1U1)
  ADD_FIELD(KSgammaI_U2U2)
  
  /* get Kerr Schild gammas */
  GET_FIELD(KSgamma_D2D2)
  GET_FIELD(KSgamma_D0D2)
  GET_FIELD(KSgamma_D0D0)
  GET_FIELD(KSgamma_D0D1)
  GET_FIELD(KSgamma_D1D2)
  GET_FIELD(KSgamma_D1D1)
  GET_FIELD(KSgammaI_U0U2)
  GET_FIELD(KSgammaI_U0U0)
  GET_FIELD(KSgammaI_U0U1)
  GET_FIELD(KSgammaI_U1U2)
  GET_FIELD(KSgammaI_U1U1)
  GET_FIELD(KSgammaI_U2U2)
  
  /* Kerr Shild lapse and shift*/
  ADD_FIELD(KSalpha)
  GET_FIELD(KSalpha)
  
  ADD_FIELD(KSB_D0)
  ADD_FIELD(KSB_D1)
  ADD_FIELD(KSB_D2)
  GET_FIELD(KSB_D0)
  GET_FIELD(KSB_D1)
  GET_FIELD(KSB_D2)
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double x   = patch->node[ijk]->x[0];
    double y   = patch->node[ijk]->x[1]-C_BH;
    double z   = patch->node[ijk]->x[2];
    double r2 = SQR(x)+SQR(y)+SQR(z);
    double rbar2  = 0.5*(r2-a2+sqrt(SQR(r2-a2)+4*a2*SQR(z)));
    double rbar   = sqrt(rbar2);
    
    k0 = (rbar*x+a*y)/(rbar2+a2);
    k1 = (rbar*y-a*x)/(rbar2+a2);
    k2 = z/rbar;
    H  = M_BH*rbar/(rbar2+a2*SQR(k2));
    double C = 2.*H;
    double A = 1./(1+C*(SQR(k0)+SQR(k1)+SQR(k2)));
    
    KSalpha[ijk] = 1/sqrt(1+C);
    KSB_D0[ijk]  = C*k0;
    KSB_D1[ijk]  = C*k1;
    KSB_D2[ijk]  = C*k2;
    
    KSgamma_D0D0[ijk] = 1.+C*k0*k0;
    KSgamma_D0D1[ijk] = C*k0*k1;
    KSgamma_D0D2[ijk] = C*k0*k2;
    KSgamma_D1D1[ijk] = 1+C*k1*k1;
    KSgamma_D1D2[ijk] = C*k1*k2;
    KSgamma_D2D2[ijk] = 1+C*k2*k2;
    
    KSgammaI_U0U0[ijk] = A*(1+C*(SQR(k1)+SQR(k2)));
    KSgammaI_U0U1[ijk] = -A*C*k0*k1;
    KSgammaI_U0U2[ijk] = -A*C*k0*k2;
    KSgammaI_U1U1[ijk] = A*(1+C*(SQR(k0)+SQR(k2)));
    KSgammaI_U1U2[ijk] = -A*C*k1*k2;
    KSgammaI_U2U2[ijk] = A*(1+C*(SQR(k0)+SQR(k1)));
    
    /* quick test check KSgamma * KSgammaI = delta */
    if (0)
    {
        double delta_U0D0 = 
      KSgammaI_U0U0[ijk]*KSgamma_D0D0[ijk] + KSgammaI_U0U1[ijk]*
      KSgamma_D0D1[ijk] + KSgammaI_U0U2[ijk]*KSgamma_D0D2[ijk];

        double delta_U0D1 = 
      KSgammaI_U0U0[ijk]*KSgamma_D0D1[ijk] + KSgammaI_U0U1[ijk]*
      KSgamma_D1D1[ijk] + KSgammaI_U0U2[ijk]*KSgamma_D1D2[ijk];

        double delta_U0D2 = 
      KSgammaI_U0U0[ijk]*KSgamma_D0D2[ijk] + KSgammaI_U0U1[ijk]*
      KSgamma_D1D2[ijk] + KSgammaI_U0U2[ijk]*KSgamma_D2D2[ijk];

        double delta_U1D2 = 
      KSgammaI_U0U1[ijk]*KSgamma_D0D2[ijk] + KSgammaI_U1U1[ijk]*
      KSgamma_D1D2[ijk] + KSgammaI_U1U2[ijk]*KSgamma_D2D2[ijk];

        double delta_U1D0 = 
      KSgammaI_U0U1[ijk]*KSgamma_D0D0[ijk] + KSgammaI_U1U1[ijk]*
      KSgamma_D0D1[ijk] + KSgammaI_U1U2[ijk]*KSgamma_D0D2[ijk];

       double delta_U1D1 = 
      KSgammaI_U0U1[ijk]*KSgamma_D0D1[ijk] + KSgammaI_U1U1[ijk]*
      KSgamma_D1D1[ijk] + KSgammaI_U1U2[ijk]*KSgamma_D1D2[ijk];

        double delta_U2D2 = 
      KSgammaI_U0U2[ijk]*KSgamma_D0D2[ijk] + KSgammaI_U1U2[ijk]*
      KSgamma_D1D2[ijk] + KSgammaI_U2U2[ijk]*KSgamma_D2D2[ijk];

        double delta_U2D0 = 
      KSgammaI_U0U2[ijk]*KSgamma_D0D0[ijk] + KSgammaI_U1U2[ijk]*
      KSgamma_D0D1[ijk] + KSgammaI_U2U2[ijk]*KSgamma_D0D2[ijk];

        double delta_U2D1 = 
      KSgammaI_U0U2[ijk]*KSgamma_D0D1[ijk] + KSgammaI_U1U2[ijk]*
      KSgamma_D1D1[ijk] + KSgammaI_U2U2[ijk]*KSgamma_D1D2[ijk];

      if(!EQL(delta_U1D1,1))  abortEr("KSgammaI is not correct!\n");
      if(!EQL(delta_U0D1,0))  abortEr("KSgammaI is not correct!\n");
      if(!EQL(delta_U0D2,0))  abortEr("KSgammaI is not correct!\n");
      if(!EQL(delta_U1D2,0))  abortEr("KSgammaI is not correct!\n");
      if(!EQL(delta_U0D0,1))  abortEr("KSgammaI is not correct!\n");
      if(!EQL(delta_U2D1,0))  abortEr("KSgammaI is not correct!\n");
      if(!EQL(delta_U2D2,1))  abortEr("KSgammaI is not correct!\n");
      if(!EQL(delta_U2D0,0))  abortEr("KSgammaI is not correct!\n");
      if(!EQL(delta_U1D0,0))  abortEr("KSgammaI is not correct!\n");

    }/* end of if(0 or 1) */
    
  }/* end of for (ijk = 0; ijk < nn; ++ijk) */
  
}

/* populating Kerr Schild Chirstoffer symbols */
static void populating_KSGamma(Patch_T *const patch)
{
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  ADD_FIELD_NoMem(dKSgamma_D0D0D1)
  ADD_FIELD_NoMem(dKSgamma_D0D0D0)
  ADD_FIELD_NoMem(dKSgamma_D2D2D2)
  ADD_FIELD_NoMem(dKSgamma_D0D0D2)
  ADD_FIELD_NoMem(dKSgamma_D0D2D1)
  ADD_FIELD_NoMem(dKSgamma_D1D1D0)
  ADD_FIELD_NoMem(dKSgamma_D1D1D2)
  ADD_FIELD_NoMem(dKSgamma_D2D2D0)
  ADD_FIELD_NoMem(dKSgamma_D2D2D1)
  ADD_FIELD_NoMem(dKSgamma_D0D1D0)
  ADD_FIELD_NoMem(dKSgamma_D0D1D1)
  ADD_FIELD_NoMem(dKSgamma_D0D1D2)
  ADD_FIELD_NoMem(dKSgamma_D0D2D0)
  ADD_FIELD_NoMem(dKSgamma_D1D1D1)
  ADD_FIELD_NoMem(dKSgamma_D1D2D1)
  ADD_FIELD_NoMem(dKSgamma_D1D2D2)
  ADD_FIELD_NoMem(dKSgamma_D1D2D0)
  ADD_FIELD_NoMem(dKSgamma_D0D2D2)
  
  /* make the following block local scope */
  {
    DECLARE_FIELD(dKSgamma_D0D0D1)
    DECLARE_FIELD(dKSgamma_D0D0D0)
    DECLARE_FIELD(dKSgamma_D2D2D2)
    DECLARE_FIELD(dKSgamma_D0D0D2)
    DECLARE_FIELD(dKSgamma_D0D2D1)
    DECLARE_FIELD(dKSgamma_D1D1D0)
    DECLARE_FIELD(dKSgamma_D1D1D2)
    DECLARE_FIELD(dKSgamma_D2D2D0)
    DECLARE_FIELD(dKSgamma_D2D2D1)
    DECLARE_FIELD(dKSgamma_D0D1D0)
    DECLARE_FIELD(dKSgamma_D0D1D1)
    DECLARE_FIELD(dKSgamma_D0D1D2)
    DECLARE_FIELD(dKSgamma_D0D2D0)
    DECLARE_FIELD(dKSgamma_D1D1D1)
    DECLARE_FIELD(dKSgamma_D1D2D1)
    DECLARE_FIELD(dKSgamma_D1D2D2)
    DECLARE_FIELD(dKSgamma_D1D2D0)
    DECLARE_FIELD(dKSgamma_D0D2D2)
    
    DECLARE_FIELD(KSgamma_D2D2)
    DECLARE_FIELD(KSgamma_D0D2)
    DECLARE_FIELD(KSgamma_D0D0)
    DECLARE_FIELD(KSgamma_D0D1)
    DECLARE_FIELD(KSgamma_D1D2)
    DECLARE_FIELD(KSgamma_D1D1)

    
    /* filling the values */
    dKSgamma_D0D0D0->v = Partial_Derivative(KSgamma_D0D0,"x");
    dKSgamma_D0D0D1->v = Partial_Derivative(KSgamma_D0D0,"y");
    dKSgamma_D0D0D2->v = Partial_Derivative(KSgamma_D0D0,"z");
    dKSgamma_D1D2D2->v = Partial_Derivative(KSgamma_D1D2,"z");
    dKSgamma_D1D1D0->v = Partial_Derivative(KSgamma_D1D1,"x");
    dKSgamma_D0D1D1->v = Partial_Derivative(KSgamma_D0D1,"y");
    dKSgamma_D0D1D0->v = Partial_Derivative(KSgamma_D0D1,"x");
    dKSgamma_D1D1D2->v = Partial_Derivative(KSgamma_D1D1,"z");
    dKSgamma_D0D1D2->v = Partial_Derivative(KSgamma_D0D1,"z");
    dKSgamma_D0D2D2->v = Partial_Derivative(KSgamma_D0D2,"z");
    dKSgamma_D1D2D1->v = Partial_Derivative(KSgamma_D1D2,"y");
    dKSgamma_D0D2D0->v = Partial_Derivative(KSgamma_D0D2,"x");
    dKSgamma_D0D2D1->v = Partial_Derivative(KSgamma_D0D2,"y");
    dKSgamma_D2D2D2->v = Partial_Derivative(KSgamma_D2D2,"z");
    dKSgamma_D2D2D0->v = Partial_Derivative(KSgamma_D2D2,"x");
    dKSgamma_D1D2D0->v = Partial_Derivative(KSgamma_D1D2,"x");
    dKSgamma_D1D1D1->v = Partial_Derivative(KSgamma_D1D1,"y");
    dKSgamma_D2D2D1->v = Partial_Derivative(KSgamma_D2D2,"y");
  }
  
  ADD_FIELD(KSGamma_U2D1D1)
  ADD_FIELD(KSGamma_U2D1D2)
  ADD_FIELD(KSGamma_U0D1D1)
  ADD_FIELD(KSGamma_U2D0D2)
  ADD_FIELD(KSGamma_U2D2D2)
  ADD_FIELD(KSGamma_U0D1D2)
  ADD_FIELD(KSGamma_U0D0D2)
  ADD_FIELD(KSGamma_U0D0D1)
  ADD_FIELD(KSGamma_U0D0D0)
  ADD_FIELD(KSGamma_U1D2D2)
  ADD_FIELD(KSGamma_U2D0D1)
  ADD_FIELD(KSGamma_U0D2D2)
  ADD_FIELD(KSGamma_U2D0D0)
  ADD_FIELD(KSGamma_U1D0D2)
  ADD_FIELD(KSGamma_U1D1D2)
  ADD_FIELD(KSGamma_U1D0D0)
  ADD_FIELD(KSGamma_U1D0D1)
  ADD_FIELD(KSGamma_U1D1D1)

  GET_FIELD(KSGamma_U2D1D1)
  GET_FIELD(KSGamma_U2D1D2)
  GET_FIELD(KSGamma_U0D1D1)
  GET_FIELD(KSGamma_U2D0D2)
  GET_FIELD(KSGamma_U2D2D2)
  GET_FIELD(KSGamma_U0D1D2)
  GET_FIELD(KSGamma_U0D0D2)
  GET_FIELD(KSGamma_U0D0D1)
  GET_FIELD(KSGamma_U0D0D0)
  GET_FIELD(KSGamma_U1D2D2)
  GET_FIELD(KSGamma_U2D0D1)
  GET_FIELD(KSGamma_U0D2D2)
  GET_FIELD(KSGamma_U2D0D0)
  GET_FIELD(KSGamma_U1D0D2)
  GET_FIELD(KSGamma_U1D1D2)
  GET_FIELD(KSGamma_U1D0D0)
  GET_FIELD(KSGamma_U1D0D1)
  GET_FIELD(KSGamma_U1D1D1)

  GET_FIELD(dKSgamma_D0D0D1)
  GET_FIELD(dKSgamma_D0D0D0)
  GET_FIELD(dKSgamma_D2D2D2)
  GET_FIELD(dKSgamma_D0D0D2)
  GET_FIELD(dKSgamma_D0D2D1)
  GET_FIELD(dKSgamma_D1D1D0)
  GET_FIELD(dKSgamma_D1D1D2)
  GET_FIELD(dKSgamma_D2D2D0)
  GET_FIELD(dKSgamma_D2D2D1)
  GET_FIELD(dKSgamma_D0D1D0)
  GET_FIELD(dKSgamma_D0D1D1)
  GET_FIELD(dKSgamma_D0D1D2)
  GET_FIELD(dKSgamma_D0D2D0)
  GET_FIELD(dKSgamma_D1D1D1)
  GET_FIELD(dKSgamma_D1D2D1)
  GET_FIELD(dKSgamma_D1D2D2)
  GET_FIELD(dKSgamma_D1D2D0)
  GET_FIELD(dKSgamma_D0D2D2)
  
  GET_FIELD(KSgammaI_U0U2)
  GET_FIELD(KSgammaI_U0U0)
  GET_FIELD(KSgammaI_U0U1)
  GET_FIELD(KSgammaI_U1U2)
  GET_FIELD(KSgammaI_U1U1)
  GET_FIELD(KSgammaI_U2U2)
  
  for(ijk = 0; ijk < nn; ++ijk)
  {
    double GAMMA_U1D0D0 = 
0.5*dKSgamma_D0D0D0[ijk]*KSgammaI_U0U1[ijk] - 0.5*KSgammaI_U1U1[ijk]*
(dKSgamma_D0D0D1[ijk] - 2*dKSgamma_D0D1D0[ijk]) - 0.5*KSgammaI_U1U2[ijk]*
(dKSgamma_D0D0D2[ijk] - 2*dKSgamma_D0D2D0[ijk]);

    double GAMMA_U0D2D2 = 
0.5*dKSgamma_D2D2D2[ijk]*KSgammaI_U0U2[ijk] + 0.5*KSgammaI_U0U0[ijk]*(2*
dKSgamma_D0D2D2[ijk] - dKSgamma_D2D2D0[ijk]) + 0.5*KSgammaI_U0U1[ijk]*(2*
dKSgamma_D1D2D2[ijk] - dKSgamma_D2D2D1[ijk]);

    double GAMMA_U1D0D2 = 
0.5*dKSgamma_D0D0D2[ijk]*KSgammaI_U0U1[ijk] + 0.5*dKSgamma_D2D2D0[ijk]*
KSgammaI_U1U2[ijk] + 0.5*KSgammaI_U1U1[ijk]*(dKSgamma_D0D1D2[ijk] - 
dKSgamma_D0D2D1[ijk] + dKSgamma_D1D2D0[ijk]);

    double GAMMA_U2D2D2 = 
0.5*dKSgamma_D2D2D2[ijk]*KSgammaI_U2U2[ijk] + 0.5*KSgammaI_U0U2[ijk]*(2*
dKSgamma_D0D2D2[ijk] - dKSgamma_D2D2D0[ijk]) + 0.5*KSgammaI_U1U2[ijk]*(2*
dKSgamma_D1D2D2[ijk] - dKSgamma_D2D2D1[ijk]);

    double GAMMA_U2D1D1 = 
0.5*dKSgamma_D1D1D1[ijk]*KSgammaI_U1U2[ijk] + 0.5*KSgammaI_U0U2[ijk]*(2*
dKSgamma_D0D1D1[ijk] - dKSgamma_D1D1D0[ijk]) - 0.5*KSgammaI_U2U2[ijk]*
(dKSgamma_D1D1D2[ijk] - 2*dKSgamma_D1D2D1[ijk]);

    double GAMMA_U1D0D1 = 
0.5*dKSgamma_D0D0D1[ijk]*KSgammaI_U0U1[ijk] + 0.5*dKSgamma_D1D1D0[ijk]*
KSgammaI_U1U1[ijk] + 0.5*KSgammaI_U1U2[ijk]*(-dKSgamma_D0D1D2[ijk] + 
dKSgamma_D0D2D1[ijk] + dKSgamma_D1D2D0[ijk]);

    double GAMMA_U1D2D2 = 
0.5*dKSgamma_D2D2D2[ijk]*KSgammaI_U1U2[ijk] + 0.5*KSgammaI_U0U1[ijk]*(2*
dKSgamma_D0D2D2[ijk] - dKSgamma_D2D2D0[ijk]) + 0.5*KSgammaI_U1U1[ijk]*(2*
dKSgamma_D1D2D2[ijk] - dKSgamma_D2D2D1[ijk]);

    double GAMMA_U2D1D2 = 
0.5*dKSgamma_D1D1D2[ijk]*KSgammaI_U1U2[ijk] + 0.5*dKSgamma_D2D2D1[ijk]*
KSgammaI_U2U2[ijk] + 0.5*KSgammaI_U0U2[ijk]*(dKSgamma_D0D1D2[ijk] + 
dKSgamma_D0D2D1[ijk] - dKSgamma_D1D2D0[ijk]);

    double GAMMA_U2D0D2 = 
0.5*dKSgamma_D0D0D2[ijk]*KSgammaI_U0U2[ijk] + 0.5*dKSgamma_D2D2D0[ijk]*
KSgammaI_U2U2[ijk] + 0.5*KSgammaI_U1U2[ijk]*(dKSgamma_D0D1D2[ijk] - 
dKSgamma_D0D2D1[ijk] + dKSgamma_D1D2D0[ijk]);

    double GAMMA_U2D0D1 = 
0.5*dKSgamma_D0D0D1[ijk]*KSgammaI_U0U2[ijk] + 0.5*dKSgamma_D1D1D0[ijk]*
KSgammaI_U1U2[ijk] + 0.5*KSgammaI_U2U2[ijk]*(-dKSgamma_D0D1D2[ijk] + 
dKSgamma_D0D2D1[ijk] + dKSgamma_D1D2D0[ijk]);

    double GAMMA_U2D0D0 = 
0.5*dKSgamma_D0D0D0[ijk]*KSgammaI_U0U2[ijk] - 0.5*KSgammaI_U1U2[ijk]*
(dKSgamma_D0D0D1[ijk] - 2*dKSgamma_D0D1D0[ijk]) - 0.5*KSgammaI_U2U2[ijk]*
(dKSgamma_D0D0D2[ijk] - 2*dKSgamma_D0D2D0[ijk]);

    double GAMMA_U0D0D2 = 
0.5*dKSgamma_D0D0D2[ijk]*KSgammaI_U0U0[ijk] + 0.5*dKSgamma_D2D2D0[ijk]*
KSgammaI_U0U2[ijk] + 0.5*KSgammaI_U0U1[ijk]*(dKSgamma_D0D1D2[ijk] - 
dKSgamma_D0D2D1[ijk] + dKSgamma_D1D2D0[ijk]);

    double GAMMA_U1D1D1 = 
0.5*dKSgamma_D1D1D1[ijk]*KSgammaI_U1U1[ijk] + 0.5*KSgammaI_U0U1[ijk]*(2*
dKSgamma_D0D1D1[ijk] - dKSgamma_D1D1D0[ijk]) - 0.5*KSgammaI_U1U2[ijk]*
(dKSgamma_D1D1D2[ijk] - 2*dKSgamma_D1D2D1[ijk]);

    double GAMMA_U0D1D1 = 
0.5*dKSgamma_D1D1D1[ijk]*KSgammaI_U0U1[ijk] + 0.5*KSgammaI_U0U0[ijk]*(2*
dKSgamma_D0D1D1[ijk] - dKSgamma_D1D1D0[ijk]) - 0.5*KSgammaI_U0U2[ijk]*
(dKSgamma_D1D1D2[ijk] - 2*dKSgamma_D1D2D1[ijk]);

    double GAMMA_U0D0D0 = 
0.5*dKSgamma_D0D0D0[ijk]*KSgammaI_U0U0[ijk] - 0.5*KSgammaI_U0U1[ijk]*
(dKSgamma_D0D0D1[ijk] - 2*dKSgamma_D0D1D0[ijk]) - 0.5*KSgammaI_U0U2[ijk]*
(dKSgamma_D0D0D2[ijk] - 2*dKSgamma_D0D2D0[ijk]);

    double GAMMA_U1D1D2 = 
0.5*dKSgamma_D1D1D2[ijk]*KSgammaI_U1U1[ijk] + 0.5*dKSgamma_D2D2D1[ijk]*
KSgammaI_U1U2[ijk] + 0.5*KSgammaI_U0U1[ijk]*(dKSgamma_D0D1D2[ijk] + 
dKSgamma_D0D2D1[ijk] - dKSgamma_D1D2D0[ijk]);

    double GAMMA_U0D0D1 = 
0.5*dKSgamma_D0D0D1[ijk]*KSgammaI_U0U0[ijk] + 0.5*dKSgamma_D1D1D0[ijk]*
KSgammaI_U0U1[ijk] + 0.5*KSgammaI_U0U2[ijk]*(-dKSgamma_D0D1D2[ijk] + 
dKSgamma_D0D2D1[ijk] + dKSgamma_D1D2D0[ijk]);

    double GAMMA_U0D1D2 = 
0.5*dKSgamma_D1D1D2[ijk]*KSgammaI_U0U1[ijk] + 0.5*dKSgamma_D2D2D1[ijk]*
KSgammaI_U0U2[ijk] + 0.5*KSgammaI_U0U0[ijk]*(dKSgamma_D0D1D2[ijk] + 
dKSgamma_D0D2D1[ijk] - dKSgamma_D1D2D0[ijk]);


    /* populating: */
    KSGamma_U2D1D1[ijk] = GAMMA_U2D1D1;
    KSGamma_U2D1D2[ijk] = GAMMA_U2D1D2;
    KSGamma_U0D1D1[ijk] = GAMMA_U0D1D1;
    KSGamma_U2D0D2[ijk] = GAMMA_U2D0D2;
    KSGamma_U2D2D2[ijk] = GAMMA_U2D2D2;
    KSGamma_U0D1D2[ijk] = GAMMA_U0D1D2;
    KSGamma_U0D0D2[ijk] = GAMMA_U0D0D2;
    KSGamma_U0D0D1[ijk] = GAMMA_U0D0D1;
    KSGamma_U0D0D0[ijk] = GAMMA_U0D0D0;
    KSGamma_U1D2D2[ijk] = GAMMA_U1D2D2;
    KSGamma_U2D0D1[ijk] = GAMMA_U2D0D1;
    KSGamma_U0D2D2[ijk] = GAMMA_U0D2D2;
    KSGamma_U2D0D0[ijk] = GAMMA_U2D0D0;
    KSGamma_U1D0D2[ijk] = GAMMA_U1D0D2;
    KSGamma_U1D1D2[ijk] = GAMMA_U1D1D2;
    KSGamma_U1D0D0[ijk] = GAMMA_U1D0D0;
    KSGamma_U1D0D1[ijk] = GAMMA_U1D0D1;
    KSGamma_U1D1D1[ijk] = GAMMA_U1D1D1;
  }/*end of for(ijk = 0; ijk < nn; ++ijk)*/
}

/* populate KS trKij */    
static void populate_KS_trKij(Patch_T *const patch)
{
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  GET_FIELD(K)
  GET_FIELD(KSalpha)
  GET_FIELD(dKSB_D2D0)
  GET_FIELD(dKSB_D0D0)
  GET_FIELD(dKSB_D2D2)
  GET_FIELD(dKSB_D0D1)
  GET_FIELD(dKSB_D0D2)
  GET_FIELD(dKSB_D1D2)
  GET_FIELD(dKSB_D1D1)
  GET_FIELD(dKSB_D1D0)
  GET_FIELD(dKSB_D2D1)
  GET_FIELD(KSGamma_U2D1D1)
  GET_FIELD(KSGamma_U2D1D2)
  GET_FIELD(KSGamma_U0D1D1)
  GET_FIELD(KSGamma_U2D0D2)
  GET_FIELD(KSGamma_U2D2D2)
  GET_FIELD(KSGamma_U0D1D2)
  GET_FIELD(KSGamma_U0D0D2)
  GET_FIELD(KSGamma_U0D0D1)
  GET_FIELD(KSGamma_U0D0D0)
  GET_FIELD(KSGamma_U1D2D2)
  GET_FIELD(KSGamma_U2D0D1)
  GET_FIELD(KSGamma_U0D2D2)
  GET_FIELD(KSGamma_U2D0D0)
  GET_FIELD(KSGamma_U1D0D2)
  GET_FIELD(KSGamma_U1D1D2)
  GET_FIELD(KSGamma_U1D0D0)
  GET_FIELD(KSGamma_U1D0D1)
  GET_FIELD(KSGamma_U1D1D1)
  GET_FIELD(KSgammaI_U0U2)
  GET_FIELD(KSgammaI_U0U0)
  GET_FIELD(KSgammaI_U0U1)
  GET_FIELD(KSgammaI_U1U2)
  GET_FIELD(KSgammaI_U1U1)
  GET_FIELD(KSgammaI_U2U2)
  GET_FIELD(KSB_D0)
  GET_FIELD(KSB_D1)
  GET_FIELD(KSB_D2)
  
    for(ijk = 0; ijk < nn; ++ijk)
    {
    double DB_D0D2 = 
-KSB_D0[ijk]*KSGamma_U0D0D2[ijk] - KSB_D1[ijk]*KSGamma_U1D0D2[ijk] - 
KSB_D2[ijk]*KSGamma_U2D0D2[ijk] + dKSB_D0D2[ijk];

    double DB_D0D0 = 
-KSB_D0[ijk]*KSGamma_U0D0D0[ijk] - KSB_D1[ijk]*KSGamma_U1D0D0[ijk] - 
KSB_D2[ijk]*KSGamma_U2D0D0[ijk] + dKSB_D0D0[ijk];

    double DB_D0D1 = 
-KSB_D0[ijk]*KSGamma_U0D0D1[ijk] - KSB_D1[ijk]*KSGamma_U1D0D1[ijk] - 
KSB_D2[ijk]*KSGamma_U2D0D1[ijk] + dKSB_D0D1[ijk];

    double DB_D2D0 = 
-KSB_D0[ijk]*KSGamma_U0D0D2[ijk] - KSB_D1[ijk]*KSGamma_U1D0D2[ijk] - 
KSB_D2[ijk]*KSGamma_U2D0D2[ijk] + dKSB_D2D0[ijk];

    double DB_D2D1 = 
-KSB_D0[ijk]*KSGamma_U0D1D2[ijk] - KSB_D1[ijk]*KSGamma_U1D1D2[ijk] - 
KSB_D2[ijk]*KSGamma_U2D1D2[ijk] + dKSB_D2D1[ijk];

    double DB_D2D2 = 
-KSB_D0[ijk]*KSGamma_U0D2D2[ijk] - KSB_D1[ijk]*KSGamma_U1D2D2[ijk] - 
KSB_D2[ijk]*KSGamma_U2D2D2[ijk] + dKSB_D2D2[ijk];

    double DB_D1D2 = 
-KSB_D0[ijk]*KSGamma_U0D1D2[ijk] - KSB_D1[ijk]*KSGamma_U1D1D2[ijk] - 
KSB_D2[ijk]*KSGamma_U2D1D2[ijk] + dKSB_D1D2[ijk];

    double DB_D1D1 = 
-KSB_D0[ijk]*KSGamma_U0D1D1[ijk] - KSB_D1[ijk]*KSGamma_U1D1D1[ijk] - 
KSB_D2[ijk]*KSGamma_U2D1D1[ijk] + dKSB_D1D1[ijk];

    double DB_D1D0 = 
-KSB_D0[ijk]*KSGamma_U0D0D1[ijk] - KSB_D1[ijk]*KSGamma_U1D0D1[ijk] - 
KSB_D2[ijk]*KSGamma_U2D0D1[ijk] + dKSB_D1D0[ijk];

    double K_D2D2 = 
DB_D2D2/KSalpha[ijk];

    double K_D2D0 = 
0.5*(DB_D0D2 + DB_D2D0)/KSalpha[ijk];

    double K_D2D1 = 
0.5*(DB_D1D2 + DB_D2D1)/KSalpha[ijk];

    double K_D0D0 = 
DB_D0D0/KSalpha[ijk];

    double K_D0D1 = 
0.5*(DB_D0D1 + DB_D1D0)/KSalpha[ijk];

    double K_D1D1 = 
DB_D1D1/KSalpha[ijk];

    double K_D1D0 = 
0.5*(DB_D0D1 + DB_D1D0)/KSalpha[ijk];

    double K_D0D2 = 
0.5*(DB_D0D2 + DB_D2D0)/KSalpha[ijk];

    double K_D1D2 = 
0.5*(DB_D1D2 + DB_D2D1)/KSalpha[ijk];

    double trK = 
KSgammaI_U0U0[ijk]*K_D0D0 + KSgammaI_U0U1[ijk]*K_D0D1 +
KSgammaI_U0U1[ijk]*K_D1D0 + KSgammaI_U0U2[ijk]*K_D0D2 +
KSgammaI_U0U2[ijk]*K_D2D0 + KSgammaI_U1U1[ijk]*K_D1D1 +
KSgammaI_U1U2[ijk]*K_D1D2 + KSgammaI_U1U2[ijk]*K_D2D1 +
KSgammaI_U2U2[ijk]*K_D2D2;

      K[ijk] = trK;
    }/*end of for(ijk = 0; ijk < nn; ++ijk)*/
}
