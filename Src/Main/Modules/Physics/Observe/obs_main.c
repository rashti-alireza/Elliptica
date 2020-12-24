/*
// Alireza Rashti
// September 2019
*/

/* synopsis:
// =========
//
// * initialize observable *
// observe(phys,sq,ret); # the observed value saves in ret 
//                       # and sq is one of the followings:
//
// * list of quantities *
// ======================
//
// "ADM(P,J)" #=> compute P and J ADM for the given physics
//        Note: definition of ADM quantities are only held when
//        asymptotic flat gauge conditions are satisfied.
//        (arXiv:gr-qc/0703035v1) a counter example is Painleve-Gullstrand.
//        the used method is very sensitive to fluctuation of the 
//        volume integrand. But, if the metric is conformally flat
//        and the slice is maximal there would be no error.
//        Note: ADM(P,J)|BH is a very crude approximation.
//        Note: for the best result multiple split at the outermost 
//        patches recommended.
// ----------------------------------------------------------------------
//
// "Komar(M)" #=> compute Komar mass for the given physics
// ----------------------------------------------------------------------
// "ADM(M)|method"   #=> compute Komar ADM for the given physics
//        Note: definition of ADM quantities are only held when
//        asymptotic flat gauge conditions are satisfied.
//        (arXiv:gr-qc/0703035v1) a counter example is Painleve-Gullstrand.
//
//        methods: [S_inf,S+V]
//        S_inf: carry out integral at a S2 surface at inf not very accurate.
//               Notes:
//               this method utilizes conformal factor and conformal
//               metric, thus if conformal factor is const and metric 
//               is flat this gives zero.
//
//        S+V  : use Gauss lema to carry out the integral in a 
//               volume and a closer surface than inf. 
//               Notes:
//               this is a gauge dependent and to be valid g_{ij} 
//               must be decreasing as O(r^-1) and Kij as O(r^-2). 
//               examples that not working with this method:
//               Painleve-Gullstrand BH which fails g_{ij} condition, 
//               furthermore, this method depends on conformal 
//               factor psi thus if psi = const => 0.
// ----------------------------------------------------------------------
//
// "Irreducible(M)" #=> irreducible mass for the givne physics (BH)
// ----------------------------------------------------------------------
//
// "CM"          #=> compute the center of mass for the given physics
// ----------------------------------------------------------------------
//
// "Spin|method" #=> compute spin for the given physics
//       methods:
//       Campanelli: gr-qc/0612076v4
//       JRB:        Phys. Rev. D 100, 124046
//       AKV:        Phys.Rev.D78:084017,2008
// ----------------------------------------------------------------------
//
*/


#include "obs_main.h"

/* calculate the quantity of interest sq based on given physics 
// and then return in ret variable. 
// NOTE: one must provide enough memory for return value in ret; 
// for instace for "ADM(P,J)" one needs 6 double type memory.
// NOTE: the order of population for "ADM(P,J)" is:
// Px = ret[0], py = ret[1], pz = ret[2]
// Jx = ret[3], Jy = ret[4], Jz = ret[5]
// and generally for index quantities (x,y,z) fills with order (0,1,2). 
// NOTE: for observe(bh,"Irreducible(M)",ret) we have:
// ret[0] = bh irreducible mass and ret[1] = AH physical(proper) area. */
int observe(Physics_T *const phys,const char *const sq,double *const ret)
{
  /* set default parameters: */
  /* the proportion of (patch->n[2]-1) for surface and 
  // volume integral used in ADM(P,J) at outermost patches. */
  Pset_default(P_"ADM_PJ_outermost_factor","0.5");
  
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
    
  struct items_S **adm = obs->items;
  Uint i;
  
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
    Free(adm[i]->g00);
    Free(adm[i]->g01);
    Free(adm[i]->g02);
    Free(adm[i]->g11);
    Free(adm[i]->g12);
    Free(adm[i]->g22);
    Free(adm[i]->n_U0);
    Free(adm[i]->n_U1);
    Free(adm[i]->n_U2);
    
    Free(adm[i]);
  }
  Free(adm);
  free(obs); 
}

