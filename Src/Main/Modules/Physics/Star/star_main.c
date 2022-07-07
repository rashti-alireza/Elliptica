/*
// Alireza Rashti
// November 2020
*/

/* star general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. */

#include "star_main.h"

/* main function to issue commands */
int star_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case STAR_TUNE_EULER_CONST:
      AssureType(phys->ctype == NS);
      ret = tune_star_Euler_constant(phys);
    break;
    
    case STAR_EXTRAPOLATE_MATTERS:
      AssureType(phys->ctype == NS);
      ret = extrapolate_star_matter(phys);
    break;
    
    case STAR_FIND_SURFACE:
      AssureType(phys->ctype == NS);
      ret = find_star_surface(phys);
    break;
    
    case STAR_TUNE_FORCE_BALANCE:
      AssureType(phys->ctype == NS);
      ret = tune_star_force_balance_equation(phys);
    break;
    
    case STAR_TUNE_CENTER:
      AssureType(phys->ctype == NS);
      ret = tune_star_center(phys);
    break;
    
    case STAR_SET_PARAMS:
      ret = set_star_params(phys);
    break;
    
    case STAR_ADD_FIELDS:
      ret = add_star_fields(phys);
    break;
    
    case STAR_START:
      AssureType(phys->ctype == NS);
      ret = start_off_star(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* tune Euler constant in fluid eqs. */
static int tune_star_Euler_constant(Physics_T *const phys)
{
  FUNC_TIC
  
  int ret = EXIT_SUCCESS;
  
  if (Pcmps(P_"type","NS"))
  {
    if (Pcmps(P_"NS_fluid","ideal_fluid") && 
        Pcmps(P_"NS_gConf","general"))
    {
      IF_sval("Euler_const_criteria","fix_baryonic_mass")
        ret = star_NS_ifluid_gConf_find_EulerC_fix_baryon_mass(phys);
      else
        Error0(NO_OPTION);
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
 
  FUNC_TOC
  return ret;
}

/* extrapolate matter where there is no matter! */
static int extrapolate_star_matter(Physics_T *const phys)
{
  FUNC_TIC
  
  int ret = EXIT_SUCCESS;
  
  if (Pcmps(P_"type","NS"))
  {
    if (Pcmps(P_"NS_fluid","ideal_fluid") && 
        Pcmps(P_"NS_gConf","general"))
    {
      ret = star_NS_idealfluid_extrapolate_matter_fields(phys);
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
 
  FUNC_TOC
  return ret;
}

/*  find star surface */
static int find_star_surface(Physics_T *const phys)
{
  FUNC_TIC
  
  int ret = EXIT_SUCCESS;
  
  if (Pcmps(P_"type","NS"))
  {
    if (Pcmps(P_"NS_fluid","ideal_fluid") && 
        Pcmps(P_"NS_gConf","general"))
    {
      ret = star_NS_find_star_surface(phys);
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
 
  FUNC_TOC
  return ret;
}

/* tune force balance eq. */
static int tune_star_force_balance_equation(Physics_T *const phys)
{
  FUNC_TIC
  
  int ret = EXIT_SUCCESS;
  
  if (Pcmps(P_"type","NS"))
  {
    if (Pcmps(P_"NS_fluid","ideal_fluid") && 
        Pcmps(P_"NS_gConf","general"))
    {
      ret = star_NS_idealfluid_gConf_force_balance(phys);
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
 
  FUNC_TOC
  return ret;
}

/* keep NS center fixed */
static int tune_star_center(Physics_T *const phys)
{
  FUNC_TIC
  
  int ret = EXIT_SUCCESS;
  
  if (Pcmps(P_"type","NS"))
  {
    if (Pcmps(P_"NS_fluid","ideal_fluid") && 
        Pcmps(P_"NS_gConf","general"))
    {
      ret = star_NS_keep_center_fixed(phys);
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
 
  FUNC_TOC
  return ret;
}

/* star add fields */
static int add_star_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  int ret = EXIT_SUCCESS;
  
  assert(phys->grid);
  
  if (Pcmps(P_"type","NS"))
  {
    if (Pcmps(P_"NS_fluid","ideal_fluid") && 
        Pcmps(P_"NS_gConf","general"))
    {
      star_NS_idealfluid_gConf_add_fields(phys->grid);
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
 
  FUNC_TOC
  return ret;
}

/* set default paramters */
static int set_star_params(Physics_T *const phys)
{
  FUNC_TIC
  
  /* star type:
  // options:
  // NS: neutron star  */
  Pset_default(P_"type","NS");

  
  if (Pcmps(P_"type","NS"))
  {
    /* fluid type:
    // options:
    // ideal_fluid: like: Phys. Rev. D 100, 124046  */
    Pset_default(P_"NS_fluid","ideal_fluid");
    
    /* conformal metric type: 
    // options:
    // flat   :  => gConf  = delta_{ij},
    // general:  => general gConf. */
    Pset_default(P_"NS_gConf","general");
    
    /* soft parameters: */
    
    /* how to start off at the very beginning:
    //
    // param:
    // ======
    // start_off
    //
    // options:
    // ========
    // TOV : using a TOV star. */
    
    /* how NS surface looks like:
    //
    // param:
    // ======
    // surface_type
    //
    // options:
    // ========
    // perfect_s2  : pefect sphere 
    // topology_s2 : a general surface with s2 topology. */
    
    /* l max for Ylm expansion
    //
    // param:
    // ======
    // surface_Ylm_max_l .*/
    
    /* tracking the number of surface finder failings. 
    // this stops the code if something is wrong!
    //
    // param:
    // ======
    // surface_max_fail: max allowed surface fails.
    // surface_num_fail: number of surface fails.
    */
    
    /* how to extrapolate matter fields outside the NS :
    //
    // param:
    // ======
    // extrapolate_matter_fields
    //
    // options:
    // ========
    // poly2: C^2 continuity across the boundary using a 2nd-order polynomial.
    // exp2 : C^2 continuity across the boundary using a exponential.
    // inverse_r2 : C^2 continuity across the boundary using a+b/r+c/r^2.
    // inverse_r2_expmr: C^2 continuity across the boundary using (a+b/r+c/r^2)*exp(-r/r0).
    // inverse_r2_expmAr: C^2 continuity across the boundary using (a+b/r+c/r^2)*exp(-Att*(r-r0)).
    // inverse_r_expmAr: C^1 continuity across the boundary using (a+b/r)*exp(-Att*(r-r0)).
    // inverse_r_expmr: C^1 continuity across the boundary using (a+b/r)*exp(-Att*(r/r0)).
    // expmr: C^0 continuity across the boundary using exp(-att*(r-r0)).(not implemented separately)
    // enthalpy_expmr_phi_inverse_r2: use expmr for enthalpy and inverse_r2 for phi. */
    
    /* which root finder to be used to find NS surface:
    //
    // param:
    // ======
    // surface_finder:
    //
    // options:
    // ========
    // bisection: using bisection root finder to find enthalpy = 1. 
    // note: if the surface_type is perfect_s2 the surface is set 
    // by its parameters and not a root finder. */
    
    /* force balance:
    // 
    // param:
    // ======
    // force_balance_equation
    // 
    // options:
    // ========
    // none: do nothing
    // d/dx:x_CM:   adjust x_CM such that dh/dx|NS_center = 0.
    // d/dy:x_CM:   adjust x_CM such that dh/dy|NS_center = 0.
    // d/dz:x_CM:   adjust x_CM such that dh/dz|NS_center = 0.
    // d/dx:y_CM:   adjust y_CM such that dh/dx|NS_center = 0.
    // d/dy:y_CM:   adjust y_CM such that dh/dy|NS_center = 0.
    // d/dz:y_CM:   adjust y_CM such that dh/dz|NS_center = 0.
    // d/dCM:Omega: adjust Omega such that dh/dCM|NS_center = 0.
    // d/dx:Omega:  adjust Omega such that dh/dx|NS_center = 0.
    // d/dy:Omega:  adjust Omega such that dh/dy|NS_center = 0.
    // d/dz:Omega:  adjust Omega such that dh/dz|NS_center = 0.
    // exmaple: adjust(d/dy:Omega).
    */
  }
  else
    Error0(NO_OPTION);
 
  
  UNUSED(phys);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* start off a star according to the parameter */
static int start_off_star(Physics_T *const phys)
{
  FUNC_TIC
  
  IF_sval("start_off","TOV")
    star_start_off_TOV(phys);
  else
    Error0(NO_OPTION);

  FUNC_TOC
  return EXIT_SUCCESS;
}


