/*
// Alireza Rashti
// September 2019
*/

/* synopsis:
// =========
//
// * initialize observable *
// observe(phys,quantity,method,ret); # the observed value saves in ret 
//
// * for list of methods, see function set_observe_params *
//
// * list of quantities *
// ======================
//
// "ADM(P)" and "ADM(J)" #=> compute P and J ADM for the given physics
//        NOTE: the definition of ADM momentum is only held when
//        asymptotic flat gauge conditions are satisfied.
//        (arXiv:gr-qc/0703035v1) a counter example is Painleve-Gullstrand.
//        NOTE: the definition of ADM angular momentum is only held when
//        quasi-isotropic and asymptotically maximal gauge are satisfied.
//        asymptotic flat gauge conditions are satisfied.
//        (arXiv:gr-qc/0703035v1) a counter example is Painleve-Gullstrand.
//        the used method is very sensitive to fluctuation of the 
//        volume integrand. But, if the metric is conformally flat
//        and the slice is maximal there would be no error.
//        NOTE: ADM(P,J)|BH is a very crude approximation.
//        NOTE: for the best result multiple split at the outermost 
//        patches recommended.
// ----------------------------------------------------------------------
//
// "Komar(M)" #=> compute Komar mass for the given physics
//            NOTE: the definition of Komar mass is based upon a time 
//            like Killing vector, thus, if we don't this symmetry
//            we cannot use Komar mass. Also, in the integral surface
//            formula, it is assumed Killing^mu = alpha n^mu + beta^mu
//            which might not be true close to compact objects, however,
//            this assumption mostly correct when the surface integral
//            is quite far away from the compact objects, since the space-time
//            becomes Minkowski and this holds.
//            NOTE: the accuracy decreases if the integrating surface
//            is not a diffeomorphism (C^inf continuous) to S2.
//            NOTE: ADM mass is generally not equal to Komar mass unless
//            lapse -> 1 and shift -> 0 at spatial infinity in the foliation.
//
// ----------------------------------------------------------------------
// "ADM(M)"   #=> compute ADM mass for the given physics
//        NOTE: definition of ADM quantities are only held when
//        asymptotic flat gauge conditions are satisfied.
//        (arXiv:gr-qc/0703035v1) a counter example is Painleve-Gullstrand.
//
//        some explanation about parameters:
//        S_inf: carry out integral at a S2 surface at inf.
//
//        Error analysis: ~ M/R where R is outer boundary radius and M is total mass
//        this was found by testing the analytic solution of the 
//        Christodoulou mass (~4.2) of a single Kerr-Schild BH with arbitary spin 
//        direction where the outer boundary was at 1*10^5 solar_mass and we got
//        at res. 20x20x20, E_rel ~ 4*10^-5.
//        when the outer radius increases to 10^6 the error becomes E_rel ~ 10^-6.
//        Moreover, this error depends on the total mass of the system, since by increasing 
//        the BH mass (with arbitrary spin) to ~ 10 the E_rel ~ 10^-4 at outer bound R = 10^5
//        and for outer bound R = 10^6, E_rel ~ 10^-5.
//               NOTE:
//               this method utilizes conformal factor and conformal
//               metric, thus if conformal factor is const and metric 
//               is flat this gives zero(like Painleve-Gullstrand BH).
//
//        S+V  : use Gauss lema to carry out the integral in a 
//               volume and a closer surface than inf. 
//               NOTE:(ref. arXiv:gr-qc/0703035v1)
//               this is a gauge dependent and to be valid gConf_{ij} 
//               must be decreasing as O(r^-2) and K as O(r^-3). 
//               this gauges are called quasi-isotropic and 
//               asymptotically maximal gauge, respectively. as one notices
//               these are stronger than asymptotically flat gauges.
//               examples that not working with this method:
//               Painleve-Gullstrand BH and Kerr-Schild which fail
//               g_{ij} condition and the latter fails K too.
//               furthermore, this method depends on conformal 
//               factor psi thus if psi = const => 0.
//        a rough approximation of the absolute error for this method
//        (S+V,default) at res. 20x20x20 is:
//        E_abs = 6*10^-4 * total_mass.
//        this was found by testing the Christodoulou mass of a single 
//        Schwarzschild (in isotropic coords), and the source of the error, 
//        is due to resolution, since at res 40x40x40 the absolut error becomes:
//        2*10^-7*total_mass.
// ----------------------------------------------------------------------
//
// "Irreducible(M)" #=> irreducible mass for the givne physics (BH)
// ----------------------------------------------------------------------
//
// "Baryonic(M)" #=> compute baryonic mass of NS
// ----------------------------------------------------------------------
//
// "CM" #=> compute the center of mass for the given physics,
//          NOTE:CM of an object measured with respect to the system CM.
// ----------------------------------------------------------------------
//
// "Spin" #=> compute spin for the given physics
// ----------------------------------------------------------------------
//
*/


#include "obs_main.h"

/* main function for mainly initialization */
int observe_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case OBSERVE_SET_PARAMS:
      ret = set_observe_params(phys);
    break;
    
    case OBSERVE_ADD_FIELDS:
      ret = add_observe_fields(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* calculate the quantity of interest sq based on given physics 
// and then return in ret variable. 
// NOTE: one must provide enough memory for return value in ret; 
// for instace for "ADM(P)" or "ADM(J)" one needs 3 double type memory.
// NOTE: the order of population for "ADM(P)" and "ADM(J)" are:
// Px = ret[0], py = ret[1], pz = ret[2]
// Jx = ret[0], Jy = ret[1], Jz = ret[2]
// and generally for index quantities (x,y,z) fills with order (0,1,2). 
// NOTE: for observe(bh,"Irreducible(M)",ret) we have:
// ret[0] = bh irreducible mass and ret[1] = AH physical(proper) area. */
int observe(Physics_T *const phys,const char *const sq,
            const char *const method,double *const ret)
{
  if (!method)
    Error0("No method was given!\n");
    
  Observe_T *const obs = calloc(1,sizeof(*obs));
  IsNull(obs);

  obs->phys  = phys;
  obs->grid  = phys->grid;
  obs->ret   = ret;
  obs->method= method;
  
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

/* set default paramters
// NOTE: since these parameters are soft, they are prefixed with stype.
// and they are supposed to set in paramter file. 
// these soft parameters are only used by this project so, they are 
// prefixed with P_. */
static int set_observe_params(Physics_T *const phys)
{
  FUNC_TIC
  
  /* followings are soft parameters to be set in parameter file.
  // they are PREFIXED with phys->stype to instruct how calculations
  // carried out, ex: BH1_observe_Komar_M.  */
  // the method is suffixed and separated with a comma,
  // ex: method = "S_inf,default" 
  // thus for instance param "BH1_observe_Komar_M" = "S_inf,default". */
  
  /* how to compute Komar mass:
  // param:
  // ======
  // "observe_Komar_M"
  //
  // methods:
  // ========
  // default  : used arXiv:gr-qc/0703035v1
  // conformal: used Eq.7.69 at arXiv:gr-qc/0703035v1.
  //
  // options:
  // ========
  // S_inf: on a surface at infinity. 
  // S_obj: over the surface of compact object (for single physics)
  // V_obj: over the volume of compact object (for single physics)
  // S+V  : over all space and on BH surface if any. */
  
  /* how to compute ADM mass: 
  // param:
  // ======
  // "observe_ADM_M"
  //
  // methods:
  // ========
  // default: used arXiv:gr-qc/0703035v1
  //          note: for S+V this method assumes isotropic gauge.
  // 
  // options:
  // ========
  // S_inf: on a surface at infinity. 
  // S_obj: over the surface of compact object (for single physics)
  // V_obj: over the volume of compact object (for single physics) (arXiv:1910.09690v1)
  // S+V  : over all space and on BH surface if any. */
  
  /* how to compute irreducible mass:
  // param:
  // ======
  // "observe_irreducible_M"
  //
  // methods:
  // ========
  // default: used standard definition
  //
  // options:
  // ========
  // S_obj: over the surface of compact object (for single physics) */
  
  /* how to compute baryonic mass:
  // param:
  // ======
  // "observe_baryonic_M"
  //
  // methods:
  // ========
  // default: used standard definition
  //
  // options:
  // ========
  // V_obj: over the volume of compact object (for NS physics) */
  
  /* how to compute ADM angular momentum: 
  // param:
  // ======
  // "observe_ADM_J"
  //
  // methods:
  // ========
  // default    : used arXiv:gr-qc/0703035v1
  // Ossokine   : used arXiv:1506.01689
  // constraint : used Stokes and assuming momentum constraints are met.
  //
  // options:
  // ========
  // S_inf: on a surface at infinity. 
  // S_obj: over the surface of compact object (for single physics)
  // S_obj1+S_obj2: over surfaces of compact object 1 and 2 */
  
  /* how to compute ADM momentum:
  // param:
  // ======
  // "observe_ADM_P"
  //
  // methods:
  // ========
  // default    : used arXiv:gr-qc/0703035v1
  // Ossokine   : used arXiv:1506.01689
  // Rashti     : used full Stokes theorem
  // constraint : used Stokes and assuming momentum constraints are met.
  //
  // options:
  // ========
  // S_inf: on a surface at infinity. 
  // S_obj: over the surface of compact object (for single physics)
  // S+V  : over surfaces of compact object 1 and 2 and the rest of space.
  // S_obj1+S_obj2: over surfaces of compact object 1 and 2 */
  
  /* how to compute Spin:
  // param:
  // ======
  // "Observe_Spin"
  //
  // methods:
  // ========
  // Campanelli: gr-qc/0612076v4
  // JRB:        Phys. Rev. D 100, 124046
  // AKV:        Phys.Rev.D78:084017,2008
  //
  // options:
  // ========
  // S_obj: over the surface of compact object (for single physics) */
  
  /* how to compute CM of an object
  // param:
  // ======
  // "observe_CM"
  // 
  // methods:
  // default: using classical mechanics definition
  // 
  // options:
  // ========
  // S_obj: over the surface of the object (mainly for BH)
  // V_obj: over the volume of the object (mainly for NS)
  */

  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;  
}

/* add fields */
static int add_observe_fields(Physics_T *const phys)
{
  FUNC_TIC
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;  
}

