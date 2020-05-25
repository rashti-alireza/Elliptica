/*
// Alireza Rashti
// July 2019
*/

#include "bbn_free_data.h"

/* populating the free data part of initial data that we chose ourself */
void bbn_populate_free_data(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Populating free data and related ...\n");
  
  /* if it is ready */
  if (0 && Pgeti("use_previous_data"))
  {
    printf("~> Using the free data of the previous grid.\n");
    printf("} Populating free data and related ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
    return;
  }
  
  /* populate conformal metric and its inverse */
  bbn_free_data_gammas(grid);
  printf("Conformal metric and its inverse ~> Done.\n");
  
  /* Christoffer symbols made up of conformal metric */
  bbn_free_data_Gamma(grid);
  printf("Christoffel symbols (_Gamma)     ~> Done.\n");
  
  /* partial derivtive of _Gamma, used in covariant derivative and _R */
  bbn_free_data_dGamma(grid);
  printf("Partial derivatives of _Gamma    ~> Done.\n");
  
  /* Ricci scalar made up of conformal metric _gamma */
  bbn_free_data_Ricci(grid);
  printf("Ricci scalar (_R)                ~> Done.\n");
  
  /* trace of Kerr Schild extrinsic curvature */
  bbn_free_data_tr_KSKij(grid);
  printf("Trace of Kerr Schild BH:tr(K_ij) ~> Done.\n");
  
  printf("} Populating free data and related ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* partial derivtive of _Gamma, used in covariant derivative and _R */
void bbn_free_data_dGamma(Grid_T *const grid)
{
  const unsigned np = grid->np;
  unsigned p;

  OpenMP_Patch_Pragma(omp parallel for)
  for(p = 0; p < np; ++p)
  {
    Patch_T *patch = grid->patch[p];
  
    /* partial derivative of _Gamma */
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D2D2D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D2D2D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D2D2D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D0D0D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D1D2D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D0D0D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D1D2D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D1D2D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D1D1D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D0D0D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D0D2D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D0D1D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D0D1D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D0D1D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D2D2D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D0D0D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D0D0D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D0D0D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D1D2D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D1D2D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D1D2D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D0D2D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D0D2D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D0D2D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D0D2D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D1D1D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D0D2D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D0D1D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D2D2D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D1D2D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D0D1D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D0D1D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D0D1D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D1D2D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D1D1D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D1D1D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D1D1D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D2D2D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D1D1D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D1D1D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U1D1D1D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D0D1D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D0D2D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D0D2D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D0D2D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D1D2D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D0D0D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D0D0D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D0D0D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U2D1D1D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D2D2D0)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D2D2D1)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D2D2D2)
    DECLARE_AND_EMPTY_FIELD(_dGamma_U0D0D1D2)
    
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
void bbn_preparing_conformal_metric_derivatives(Patch_T *const patch)
{
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

/* freeing conformal metric derivatives */
void bbn_free_conformal_metric_derivatives(Patch_T *const patch)
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
void bbn_free_data_gammas(Grid_T *const grid)
{
  /* roll off distance at exp(-(r/r0)^4)  */
  const double r0          = Pgetd("BH_KerrSchild_RollOff");
  const double BH_center_x = Pgetd("BH_center_x");
  const double BH_center_y = Pgetd("BH_center_y");
  const double BH_center_z = Pgetd("BH_center_z");
  const double M_BH        = Pgetd("BH_irreducible_mass");
  const double a_BH        = Pgetd("BH_net_spin");
  unsigned p;
  
  /* populate tB tR */
  Transformation_T *tB = initialize_transformation();
  Transformation_T *tR = initialize_transformation();
  bbn_transform_populate_boost_rotation(tB,tR);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn,ijk;
    nn = patch->nn;
    
    REALLOC_v_WRITE_v(_gamma_D2D2)
    REALLOC_v_WRITE_v(_gamma_D0D2)
    REALLOC_v_WRITE_v(_gamma_D0D0)
    REALLOC_v_WRITE_v(_gamma_D0D1)
    REALLOC_v_WRITE_v(_gamma_D1D2)
    REALLOC_v_WRITE_v(_gamma_D1D1)
    REALLOC_v_WRITE_v(_gammaI_U0U2)
    REALLOC_v_WRITE_v(_gammaI_U0U0)
    REALLOC_v_WRITE_v(_gammaI_U0U1)
    REALLOC_v_WRITE_v(_gammaI_U1U2)
    REALLOC_v_WRITE_v(_gammaI_U1U1)
    REALLOC_v_WRITE_v(_gammaI_U2U2)
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x,y,z,r,H,k0,k1,k2,kt;
  
      x = patch->node[ijk]->x[0]-BH_center_x;
      y = patch->node[ijk]->x[1]-BH_center_y;
      z = patch->node[ijk]->x[2]-BH_center_z;
      r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
      bbn_transform_get_k_and_H_KerrSchild
        (x,y,z,a_BH,M_BH,tB,tR,&kt,&k0,&k1,&k2,&H);
      
      double e = exp(-pow(r/r0,4));
      double C = 2.*H*e;
      double A = 1./(1+C*(Pow2(k0)+Pow2(k1)+Pow2(k2)));
      
      _gamma_D0D0[ijk] = 1.+C*k0*k0;
      _gamma_D0D1[ijk] = C*k0*k1;
      _gamma_D0D2[ijk] = C*k0*k2;
      _gamma_D1D1[ijk] = 1+C*k1*k1;
      _gamma_D1D2[ijk] = C*k1*k2;
      _gamma_D2D2[ijk] = 1+C*k2*k2;
      
      _gammaI_U0U0[ijk] = A*(1+C*(Pow2(k1)+Pow2(k2)));
      _gammaI_U0U1[ijk] = -A*C*k0*k1;
      _gammaI_U0U2[ijk] = -A*C*k0*k2;
      _gammaI_U1U1[ijk] = A*(1+C*(Pow2(k0)+Pow2(k2)));
      _gammaI_U1U2[ijk] = -A*C*k1*k2;
      _gammaI_U2U2[ijk] = A*(1+C*(Pow2(k0)+Pow2(k1)));
      
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

        if(!EQL(delta_U1D1,1)||!isfinite(delta_U1D1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D1,0)||!isfinite(delta_U0D1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D2,0)||!isfinite(delta_U0D2))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U1D2,0)||!isfinite(delta_U1D2))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D0,1)||!isfinite(delta_U0D0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D1,0)||!isfinite(delta_U2D1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D2,1)||!isfinite(delta_U2D2))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D0,0)||!isfinite(delta_U2D0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U1D0,0)||!isfinite(delta_U1D0))  Error0("_gammaI is not correct!\n");

      }
      
    }
  }
  free_transformation(tB);
  free_transformation(tR);
}

/* trace of Kerr Schild extrinsic curvature */
void bbn_free_data_tr_KSKij(Grid_T *const grid)
{
  const unsigned np = grid->np;
  unsigned p;

  OpenMP_Patch_Pragma(omp parallel for)
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
    bbn_free_data_KS_trKij(patch);
    
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
  ADD_FIELD(dKSB_D2D0)
  ADD_FIELD(dKSB_D0D0)
  ADD_FIELD(dKSB_D2D2)
  ADD_FIELD(dKSB_D0D1)
  ADD_FIELD(dKSB_D0D2)
  ADD_FIELD(dKSB_D1D2)
  ADD_FIELD(dKSB_D1D1)
  ADD_FIELD(dKSB_D1D0)
  ADD_FIELD(dKSB_D2D1)
  
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
  const double M_BH        = Pgetd("BH_irreducible_mass");
  const double a_BH        = Pgetd("BH_net_spin");
  const double BH_center_x = Pgetd("BH_center_x");
  const double BH_center_y = Pgetd("BH_center_y");
  const double BH_center_z = Pgetd("BH_center_z");
  const unsigned nn        = patch->nn;
  unsigned ijk;
  
  /* populate tB tR */
  Transformation_T *tB = initialize_transformation();
  Transformation_T *tR = initialize_transformation();
  bbn_transform_populate_boost_rotation(tB,tR);
  
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
  REALLOC_v_WRITE_v(KSgamma_D2D2)
  REALLOC_v_WRITE_v(KSgamma_D0D2)
  REALLOC_v_WRITE_v(KSgamma_D0D0)
  REALLOC_v_WRITE_v(KSgamma_D0D1)
  REALLOC_v_WRITE_v(KSgamma_D1D2)
  REALLOC_v_WRITE_v(KSgamma_D1D1)
  REALLOC_v_WRITE_v(KSgammaI_U0U2)
  REALLOC_v_WRITE_v(KSgammaI_U0U0)
  REALLOC_v_WRITE_v(KSgammaI_U0U1)
  REALLOC_v_WRITE_v(KSgammaI_U1U2)
  REALLOC_v_WRITE_v(KSgammaI_U1U1)
  REALLOC_v_WRITE_v(KSgammaI_U2U2)
  
  /* Kerr Schild lapse and shift*/
  ADD_FIELD(KSalpha)
  REALLOC_v_WRITE_v(KSalpha)
  
  ADD_FIELD(KSB_D0)
  ADD_FIELD(KSB_D1)
  ADD_FIELD(KSB_D2)
  REALLOC_v_WRITE_v(KSB_D0)
  REALLOC_v_WRITE_v(KSB_D1)
  REALLOC_v_WRITE_v(KSB_D2)
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double x,y,z,H,k0,k1,k2,kt;
 
    x = patch->node[ijk]->x[0]-BH_center_x;
    y = patch->node[ijk]->x[1]-BH_center_y;
    z = patch->node[ijk]->x[2]-BH_center_z;
    
    bbn_transform_get_k_and_H_KerrSchild
      (x,y,z,a_BH,M_BH,tB,tR,&kt,&k0,&k1,&k2,&H);
    
    double C = 2.*H;
    double A = 1./(1+C*(Pow2(k0)+Pow2(k1)+Pow2(k2)));
    
    KSalpha[ijk] = 1/sqrt(1+C*kt*kt);
    KSB_D0[ijk]  = C*k0*kt;
    KSB_D1[ijk]  = C*k1*kt;
    KSB_D2[ijk]  = C*k2*kt;
    
    KSgamma_D0D0[ijk] = 1.+C*k0*k0;
    KSgamma_D0D1[ijk] = C*k0*k1;
    KSgamma_D0D2[ijk] = C*k0*k2;
    KSgamma_D1D1[ijk] = 1+C*k1*k1;
    KSgamma_D1D2[ijk] = C*k1*k2;
    KSgamma_D2D2[ijk] = 1+C*k2*k2;
    
    KSgammaI_U0U0[ijk] = A*(1+C*(Pow2(k1)+Pow2(k2)));
    KSgammaI_U0U1[ijk] = -A*C*k0*k1;
    KSgammaI_U0U2[ijk] = -A*C*k0*k2;
    KSgammaI_U1U1[ijk] = A*(1+C*(Pow2(k0)+Pow2(k2)));
    KSgammaI_U1U2[ijk] = -A*C*k1*k2;
    KSgammaI_U2U2[ijk] = A*(1+C*(Pow2(k0)+Pow2(k1)));
    
    /* quick test check KSgamma * KSgammaI = delta */
    if (1)
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

      if(!EQL(delta_U1D1,1)||!isfinite(delta_U1D1))  Error0("_gammaI is not correct!\n");
      if(!EQL(delta_U0D1,0)||!isfinite(delta_U0D1))  Error0("_gammaI is not correct!\n");
      if(!EQL(delta_U0D2,0)||!isfinite(delta_U0D2))  Error0("_gammaI is not correct!\n");
      if(!EQL(delta_U1D2,0)||!isfinite(delta_U1D2))  Error0("_gammaI is not correct!\n");
      if(!EQL(delta_U0D0,1)||!isfinite(delta_U0D0))  Error0("_gammaI is not correct!\n");
      if(!EQL(delta_U2D1,0)||!isfinite(delta_U2D1))  Error0("_gammaI is not correct!\n");
      if(!EQL(delta_U2D2,1)||!isfinite(delta_U2D2))  Error0("_gammaI is not correct!\n");
      if(!EQL(delta_U2D0,0)||!isfinite(delta_U2D0))  Error0("_gammaI is not correct!\n");
      if(!EQL(delta_U1D0,0)||!isfinite(delta_U1D0))  Error0("_gammaI is not correct!\n");

    }/* end of if(0 or 1) */
    
  }/* end of for (ijk = 0; ijk < nn; ++ijk) */
  free_transformation(tB);
  free_transformation(tR);
}

/* populating Kerr Schild Chirstoffer symbols */
static void populating_KSGamma(Patch_T *const patch)
{
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  ADD_FIELD(dKSgamma_D0D0D1)
  ADD_FIELD(dKSgamma_D0D0D0)
  ADD_FIELD(dKSgamma_D2D2D2)
  ADD_FIELD(dKSgamma_D0D0D2)
  ADD_FIELD(dKSgamma_D0D2D1)
  ADD_FIELD(dKSgamma_D1D1D0)
  ADD_FIELD(dKSgamma_D1D1D2)
  ADD_FIELD(dKSgamma_D2D2D0)
  ADD_FIELD(dKSgamma_D2D2D1)
  ADD_FIELD(dKSgamma_D0D1D0)
  ADD_FIELD(dKSgamma_D0D1D1)
  ADD_FIELD(dKSgamma_D0D1D2)
  ADD_FIELD(dKSgamma_D0D2D0)
  ADD_FIELD(dKSgamma_D1D1D1)
  ADD_FIELD(dKSgamma_D1D2D1)
  ADD_FIELD(dKSgamma_D1D2D2)
  ADD_FIELD(dKSgamma_D1D2D0)
  ADD_FIELD(dKSgamma_D0D2D2)
  
  /* force the following block to be local scope */
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

  REALLOC_v_WRITE_v(KSGamma_U2D1D1)
  REALLOC_v_WRITE_v(KSGamma_U2D1D2)
  REALLOC_v_WRITE_v(KSGamma_U0D1D1)
  REALLOC_v_WRITE_v(KSGamma_U2D0D2)
  REALLOC_v_WRITE_v(KSGamma_U2D2D2)
  REALLOC_v_WRITE_v(KSGamma_U0D1D2)
  REALLOC_v_WRITE_v(KSGamma_U0D0D2)
  REALLOC_v_WRITE_v(KSGamma_U0D0D1)
  REALLOC_v_WRITE_v(KSGamma_U0D0D0)
  REALLOC_v_WRITE_v(KSGamma_U1D2D2)
  REALLOC_v_WRITE_v(KSGamma_U2D0D1)
  REALLOC_v_WRITE_v(KSGamma_U0D2D2)
  REALLOC_v_WRITE_v(KSGamma_U2D0D0)
  REALLOC_v_WRITE_v(KSGamma_U1D0D2)
  REALLOC_v_WRITE_v(KSGamma_U1D1D2)
  REALLOC_v_WRITE_v(KSGamma_U1D0D0)
  REALLOC_v_WRITE_v(KSGamma_U1D0D1)
  REALLOC_v_WRITE_v(KSGamma_U1D1D1)

  READ_v(dKSgamma_D0D0D1)
  READ_v(dKSgamma_D0D0D0)
  READ_v(dKSgamma_D2D2D2)
  READ_v(dKSgamma_D0D0D2)
  READ_v(dKSgamma_D0D2D1)
  READ_v(dKSgamma_D1D1D0)
  READ_v(dKSgamma_D1D1D2)
  READ_v(dKSgamma_D2D2D0)
  READ_v(dKSgamma_D2D2D1)
  READ_v(dKSgamma_D0D1D0)
  READ_v(dKSgamma_D0D1D1)
  READ_v(dKSgamma_D0D1D2)
  READ_v(dKSgamma_D0D2D0)
  READ_v(dKSgamma_D1D1D1)
  READ_v(dKSgamma_D1D2D1)
  READ_v(dKSgamma_D1D2D2)
  READ_v(dKSgamma_D1D2D0)
  READ_v(dKSgamma_D0D2D2)
  
  READ_v(KSgammaI_U0U2)
  READ_v(KSgammaI_U0U0)
  READ_v(KSgammaI_U0U1)
  READ_v(KSgammaI_U1U2)
  READ_v(KSgammaI_U1U1)
  READ_v(KSgammaI_U2U2)
  
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

/* ->return value: r funciont in Kerr-Schild coords */
double bbn_KerrSchild_r(const double x,const double y,const double z,const double a)
{
  const double r2 = Pow2(x)+Pow2(y)+Pow2(z);
  const double a2 = Pow2(a);
  
  return sqrt(0.5*(r2-a2+sqrt(Pow2(r2-a2)+4*a2*Pow2(z))));
}

/* ->return value: H function in Kerr-Schild coords */
double bbn_KerrSchild_H(const double M_BH,const double rbar,const double a,const double z)
{
  /* double lambda      = 0; */
  const double lambda= 1.;
  const double k2    = z/rbar;
  const double a2    = Pow2(a);
  const double rbar2 = Pow2(rbar);
  
  /* which metric specified
    if (Pcmps("BH_NS_free_data_metric","conformally_flat_metric"))
    {
      lambda = 0;
    }
    else if (Pcmps("BH_NS_free_data_metric","Boosted_KerrSchild_metric"))
    {
      lambda = 1;
    }
    else
      Error0(NO_OPTION);
  */
  return lambda*M_BH*rbar/(rbar2+a2*Pow2(k2));
}

/* transforming 4-vector 'in' to 'out' by the rotation tR 
// followed by the boost tB, i.e. out = (tB * tR) in.
// if IsInverse is 1 then it does inverse transformation,
// i.e. out = (tB * tR)^-1 in = tR^-1 * tB ^-1 out.
// Note: We need specific order for rotation of BH which is
// out = Ry*Rz in.*/
static void execute_boost_and_rotation(Transformation_T *const tB,
                                  Transformation_T *const tR,
                                  const int IsInverse,
                                  const double in[4],
                                  double out[4])
{
  
  if (IsInverse)
  {
    double outB[4] ={0};
    double outRy[4]={0};
    
    /* inverse */
    tB->boost->Bx *= -1.;
    tB->boost->By *= -1.;
    tB->boost->Bz *= -1.;
    /* first boost transform */
    boost_transformation(tB,in,outB);
    /* bring back */
    tB->boost->Bx *= -1.;
    tB->boost->By *= -1.;
    tB->boost->Bz *= -1.;
    
    /* second rotation transfrom, the oreder used was RyRz,
    // thus, the inverse is Rz^-1 Ry^-1 */
    Transformation_T *tIR = initialize_transformation();
    tIR->rotation->active = tR->rotation->active;
    assert(EQL(tR->rotation->Rx,0));
    
    tIR->rotation->Rx = 0;
    tIR->rotation->Ry = -tR->rotation->Ry;
    tIR->rotation->Rz = 0;
    rotation_transformation(tIR,outB,outRy);
    
    tIR->rotation->Rx = 0;
    tIR->rotation->Ry = 0;
    tIR->rotation->Rz = -tR->rotation->Rz;
    rotation_transformation(tIR,outRy,out);
    
    free_transformation(tIR);
  }
  else
  {
    double outRz[4]={0};
    double outRy[4]={0};
    
    /* rotation: out = Ry Rz in */
    Transformation_T *tIR = initialize_transformation();
    tIR->rotation->active = tR->rotation->active;
    assert(EQL(tR->rotation->Rx,0));
    /* Rz */
    tIR->rotation->Rx = 0;
    tIR->rotation->Ry = 0;
    tIR->rotation->Rz = tR->rotation->Rz;
    rotation_transformation(tIR,in,outRz);
    
    /* Ry */
    tIR->rotation->Rx = 0;
    tIR->rotation->Ry = tR->rotation->Ry;
    tIR->rotation->Rz = 0;
    rotation_transformation(tIR,outRz,outRy);
    
    free_transformation(tIR);
    
    /* boost */
    boost_transformation(tB,outRy,out);
  }
}

/* populating tB(boost) and tR(rotation) for boost and rotation of BH. */
void bbn_transform_populate_boost_rotation(Transformation_T *const tB,
                        Transformation_T *const tR)
{
  const double BH_center_x = Pgetd("BH_center_x");
  const double BH_center_y = Pgetd("BH_center_y");
  //const double BH_center_z = Pgetd("BH_center_z");
  const double chi_U0   = Pgetd("BH_chi_U0");
  const double chi_U1   = Pgetd("BH_chi_U1");
  const double chi_U2   = Pgetd("BH_chi_U2");
  const double y_CM = Pgetd("y_CM");
  const double x_CM = Pgetd("x_CM");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double chi = sqrt(Pow2(chi_U0)+Pow2(chi_U1)+Pow2(chi_U2));
  double phiy,phiz;
  double Bx,By,Bz;/* B = v/c */
  
  /* boost */
  Bx = -Omega_BHNS*(BH_center_y-y_CM);
  By =  Omega_BHNS*(BH_center_x-x_CM);
  Bz = Pgetd("BH_Vz");
  tB->boost->Bx = Bx;
  tB->boost->By = By;
  tB->boost->Bz = Bz;
  tB->boost->B2 = Pow2(Bx)+Pow2(By)+Pow2(Bz);
  
  /* rotation */
  if (!EQL(chi,0))
  {
    phiz = arctan(chi_U1,chi_U0);
    phiy = acos(chi_U2/chi);
    assert(isfinite(phiy));
    tR->rotation->Rx = 0;
    tR->rotation->Ry = phiy;
    tR->rotation->Rz = phiz;
  }
} 

/* given x,y,z and transformation info of non inertial frame, 
// it gets the correspoinding k0,k1,k2 and H in Kerr-Schild context */
void bbn_transform_get_k_and_H_KerrSchild(const double x,const double y,const double z,
                                          const double a,const double M_BH,
                                          Transformation_T *const tB,
                                          Transformation_T *const tR,
                                          double *const kt,double *const k0,
                                          double *const k1,double *const k2,
                                          double *const H)

{
  const double a2 = Pow2(a);
  /*notation:
  // all of the variables with '_' prefix are in inetrial frame and
  // all of the variables without '_' prefix, are in transformed frame. */
  double _k0,_k1,_k2,_kt;/* in ds^2 = (delta_ij+2*H*ki*kj)dx^i*dx^j */
  double _r,_r2,_x,_y,_z;
  double x_mu[4] = {0/* time component */,x,y,z};
  double _x_mu[4];
      
  /* _x_mu = T^-1 x_mu */
  execute_boost_and_rotation(tB,tR,1,x_mu,_x_mu);
  
  _x  = _x_mu[1];
  _y  = _x_mu[2];
  _z  = _x_mu[3];
  _r  = bbn_KerrSchild_r(_x,_y,_z,a);
  _r2 = Pow2(_r);
  _k0 = (_r*_x+a*_y)/(_r2+a2);
  _k1 = (_r*_y-a*_x)/(_r2+a2);
  _k2 = _z/_r;
  _kt = 1;
  
  double _k_mu[4] = {_kt,_k0,_k1,_k2};
  double k_mu[4];/* Lorentz *k^mu */
  
  /* k_mu = T _k_mu */
  execute_boost_and_rotation(tB,tR,0,_k_mu,k_mu);
  
  *kt = k_mu[0];
  *k0 = k_mu[1];
  *k1 = k_mu[2];
  *k2 = k_mu[3];
  *H  = bbn_KerrSchild_H(M_BH,_r,a,_z);
}






