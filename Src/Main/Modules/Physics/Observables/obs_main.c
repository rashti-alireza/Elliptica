/*
// Alireza Rashti
// September 2019
*/

/* synopsis:
// =========
//
// * initialize observable *
// Observable_T *obs = init_observable(object,string_quantity);
// to see the list of string_quantity see "obs_integrals.c"
//
// * after initialization calculate the observable. example:*
// double Px_ADM = obs->Px(obs);# x component
// double Py_ADM = obs->Py(obs);# y component
// double Pz_ADM = obs->Pz(obs);# z component
// double Jx_ADM = obs->Jx(obs);# x component of angular momentum
// double Jy_ADM = obs->Jy(obs);# y component of angular momentum
// double Jz_ADM = obs->Jz(obs);# z component of angular momentum
// double M_ADM  = obs->M(obs) ;# a specifed mass for example ADM mass
//
// * free *
// free_observable(obs);
*/

#include "obs_main.h"

/* initialzing stuct Observable_T for sq look explanation on top. */
Observable_T *init_observable(Physics_T *const obj,const char *const sq)
{
  Observable_T *const obs = calloc(1,sizeof(*obs));
  IsNull(obs);

  obs->obj      = obj;
  obs->grid     = obj->grid;
  sprintf(obs->quantity,"%s|%s",sq,obj->stype);
  
  obs_plan(obs);
  
  return obs;
}


/* free stuct Observable_T and items */
void free_observable(Observable_T *obs)
{
  if (!obs)
    return;
    
  if (obs->grid->kind == Grid_SplitCubedSpherical_BHNS)
  {  
  
  struct items_S **adm = obs->items;
  unsigned i;
  
  for (i = 0; i < obs->Nitems; ++i)
  {
    Patch_T *patch = adm[i]->patch;
    if (patch)
    {
      if (_Ind("ADM_integrand_P_U0") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_P_U0);
        REMOVE_FIELD(ADM_integrand_P_U0);
      }
      if (_Ind("ADM_integrand_P_U1") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_P_U1);
        REMOVE_FIELD(ADM_integrand_P_U1);
      }
      if (_Ind("ADM_integrand_P_U2") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_P_U2);
        REMOVE_FIELD(ADM_integrand_P_U2);
      }
      if (_Ind("ADM_integrand_G_U0") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_G_U0);
        REMOVE_FIELD(ADM_integrand_G_U0);
      }
      if (_Ind("ADM_integrand_G_U1") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_G_U1);
        REMOVE_FIELD(ADM_integrand_G_U1);
      }
      if (_Ind("ADM_integrand_G_U2") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_G_U2);
        REMOVE_FIELD(ADM_integrand_G_U2);
      }
      
      if (_Ind("ADM_integrand_xiP_U0") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiP_U0);
        REMOVE_FIELD(ADM_integrand_xiP_U0);
      }
      if (_Ind("ADM_integrand_xiP_U2") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiP_U2);
        REMOVE_FIELD(ADM_integrand_xiP_U2);
      }
      if (_Ind("ADM_integrand_xiP_U1") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiP_U1);
        REMOVE_FIELD(ADM_integrand_xiP_U1);
      }
      
      if (_Ind("ADM_integrand_xiG_U0") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiG_U0);
        REMOVE_FIELD(ADM_integrand_xiG_U0);
      }
      if (_Ind("ADM_integrand_xiG_U2") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiG_U2);
        REMOVE_FIELD(ADM_integrand_xiG_U2);
      }
      if (_Ind("ADM_integrand_xiG_U1") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiG_U1);
        REMOVE_FIELD(ADM_integrand_xiG_U1);
      }
    }
    _free(adm[i]->g00);
    _free(adm[i]->g01);
    _free(adm[i]->g02);
    _free(adm[i]->g11);
    _free(adm[i]->g12);
    _free(adm[i]->g22);
    _free(adm[i]->n_U0);
    _free(adm[i]->n_U1);
    _free(adm[i]->n_U2);
    
    free(adm[i]);
  }
  _free(adm);
  
  }
  else
    Error0(NO_OPTION);
}

