/*
// Alireza Rashti
// September 2019
*/

/* synopsis:
// =========
//
// * initialize observable *
// observe(phys,sq,save); # the observed value saves in save 
//                          # and sq is one of the followings:
//
// * list of quantities *
// "ADM(P,J)|BHNS" #=> compute P and J ADM for the system 
// "ADM(P,J)|NS"   #=> compute P and J ADM for single NS 
// "ADM(P,J)|BH"   #=> compute P and J ADM for single BH
// "Kommar(M)|BHNS" #=> compute Kommar mass for the system 
// "Kommar(M)|NS"  #=> compute kommar mass for NS 
// "Kommar(M)|BH"  #=> compute Kommar mass for BH
// "ADM(M)|BHNS"   #=> compute ADM mass for the system 
// "ADM(M)|NS"     #=> compute ADM mass for NS 
// "ADM(M)|BH"     #=> compute ADM mass for BH
// "CM"          #=> compute the center of mass of phys->type
// "Spin|method" #=> compute spin of object phys->type
//                       with the specified method below:
//
// spin calculation methods:
// Campanelli: gr-qc/0612076v4
// JRB:        Phys. Rev. D 100, 124046
// AKV:        Phys.Rev.D78:084017,2008
//
//
*/


#include "obs_main.h"

/* calculate the quantity of interest sq based on given physics 
// and then return in ret variable. 
// NOTE: one must provide enough memory for return value in ret; 
// for instace for "ADM(P,J)" one need 6 double type memory.
// NOTE: the order of population for "ADM(P,J)" is:
// Px = ret[0], py = ret[1], pz = ret[2]
// Jx = ret[3], Jy = ret[4], Jz = ret[5]
// and generally for index quantities (x,y,z) fills with order (0,1,2). */
int observe(Physics_T *const phys,const char *const sq,double *const ret)
{
  Observe_T *const obs = calloc(1,sizeof(*obs));
  IsNull(obs);

  obs->phys  = phys;
  obs->grid  = phys->grid;
  obs->ret   = ret;
  
  sprintf(obs->quantity,"%s|%s",sq,phys->stype);
  
  obs_calculate(obs);
  free_obs(obs);
  
  return EXIT_SUCCESS;
}


/* free stuct Observe_T and items */
static void free_obs(Observe_T *obs)
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

