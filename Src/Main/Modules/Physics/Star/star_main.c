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
      ret = extrapolate_matter(phys);
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
      ret = star_add_fields(phys);
    break;
    
    //case STAR_START:
      //ret = star_NS_start_off(phys);
    //break;
    
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
      ret = star_NS_idealfluid_gConf_find_Euler_const(phys);
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
static int extrapolate_matter(Physics_T *const phys)
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
static int star_add_fields(Physics_T *const phys)
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

/* set defualt paramters */
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
    
    /* how to extrapolate matter fields outside the NS :
    // options:
    // slop_method: required to have C^2 field across the boundary. */
    Pset_default(P_"NS_extrapolate_matter_fields","poly2");
    
    /* which root finder to be used to find NS surface :
    // options:
    // slop_method: required to have C^2 field across the boundary. */
    Pset_default(P_"NS_surface_finder","bisection");
  }
  else
    Error0(NO_OPTION);
 
  
  UNUSED(phys);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

