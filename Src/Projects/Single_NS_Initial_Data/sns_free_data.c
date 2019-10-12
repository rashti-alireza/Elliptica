/*
// Alireza Rashti
// October 2019
*/

#include "sns_free_data.h"

/* populating the free data part of initial data that we chose ourself */
void sns_populate_free_data(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Populating free data and related ...\n");
  
  /* populate conformal metric and its inverse */
  sns_free_data_gammas(grid);
  printf("Conformal metric and its inverse ~> Done.\n");
  
  /* Christoffer symbols made up of conformal metric */
  sns_free_data_Gamma(grid);
  printf("Christoffer symbols (_Gamma)     ~> Done.\n");
  
  /* partial derivtive of _Gamma, used in covariant derivative and _R */
  sns_free_data_dGamma(grid);
  printf("Partial derivatives of _Gamma    ~> Done.\n");
  
  /* Ricci scalar made up of conformal metric _gamma */
  sns_free_data_Ricci(grid);
  printf("Ricci scalar (_R)                ~> Done.\n");
  
  printf("Populating free data and related ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* partial derivtive of _Gamma, used in covariant derivative and _R */
static void sns_free_data_dGamma(Grid_T *const grid)
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
void sns_preparing_conformal_metric_derivatives(Patch_T *const patch)
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
void sns_free_conformal_metric_derivatives(Patch_T *const patch)
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
void sns_free_data_gammas(Grid_T *const grid)
{
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    nn = patch->nn;
    
    PREP_FIELD(_gamma_D2D2)
    PREP_FIELD(_gamma_D0D2)
    PREP_FIELD(_gamma_D0D0)
    PREP_FIELD(_gamma_D0D1)
    PREP_FIELD(_gamma_D1D2)
    PREP_FIELD(_gamma_D1D1)
    PREP_FIELD(_gammaI_U0U2)
    PREP_FIELD(_gammaI_U0U0)
    PREP_FIELD(_gammaI_U0U1)
    PREP_FIELD(_gammaI_U1U2)
    PREP_FIELD(_gammaI_U1U1)
    PREP_FIELD(_gammaI_U2U2)
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      
      _gamma_D0D0[ijk] = 1;
      _gamma_D0D1[ijk] = 0;
      _gamma_D0D2[ijk] = 0;
      _gamma_D1D1[ijk] = 1;
      _gamma_D1D2[ijk] = 0;
      _gamma_D2D2[ijk] = 1;
      
      _gammaI_U0U0[ijk] = 1;
      _gammaI_U0U1[ijk] = 0;
      _gammaI_U0U2[ijk] = 0;
      _gammaI_U1U1[ijk] = 1;
      _gammaI_U1U2[ijk] = 0;
      _gammaI_U2U2[ijk] = 1;
      
    }
  }
}
