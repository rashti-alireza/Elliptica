/*
// Alireza Rashti
// November 2020
*/

/* star general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. 
// also in a project first mount must be called and then update. */

#include "star_main.h"

/* main function to issue commands */
int star_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  assert(phys->ctype == NS);
  
  if (Pcmps("star_type","NS"))
  {
    if (Pcmps("star_NS_fluid","ideal_fluid") && 
        Pcmps("star_NS_gConf","non_flat"))
    {
      switch (phys->cmd)
      {
        case TUNE_EULER_CONST:
          ret = star_NS_idealfluid_gConf_find_Euler_const(phys);
        break;
        case EXTRAPOLATE_MATTERS:
          ret = star_NS_idealfluid_extrapolate_matter_fields(phys);
        break;
        case FIND_STAR_SURFACE:
          ret = star_NS_find_star_surface(phys);
        break;
        case TUNE_FORCE_BALANCE:
          ret = star_NS_idealfluid_gConf_force_balance(phys);
        break;
        case TUNE_NS_CENTER:
          ret = star_NS_keep_center_fixed(phys);
        break;
        
        default:
          Error0(NO_OPTION);
      }
    }
    else
    {
      Error0(NO_OPTION);
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
    
  return ret;
}

/* adding default parameters and fields. */
int star_mount(Grid_T *const grid)
{
  /* star type:
  // NS: only neutron star,
  // WD: only white dwarf,
  // BS: only bosonic star,
  // mix(NS+BS): mix of NS and BS.
  // one can generalize this and add these parameters here 
  // and in the file for other systems. */
  Pset_default("star_type","NS");
  
  if (Pcmps("star_type","NS"))
  {
    /* fluid type:
    // ideal_fluid: like: Phys. Rev. D 100, 124046  */
    Pset_default("star_NS_fluid","ideal_fluid");
    
    /* conformal metric type: 
    // flat     => gConf  = delta_{ij},
    // non_flat => gConf != delta_{ij}. */
    Pset_default("star_NS_gConf","non_flat");
    
    /* how to extrapolate matter fields outside the NS :
    // slop_method: required to have C^2 field across the boundary. */
    Pset_default("star_NS_extrapolate_matter_fields","poly2");
    
    /* which root finder to be used to find NS surface :
    // slop_method: required to have C^2 field across the boundary. */
    Pset_default("star_NS_surface_finder","bisection");
    
    if(Pcmps("star_NS_fluid","ideal_fluid") &&
       Pcmps("star_NS_gConf","non_flat"))
    {
      star_NS_idealfluid_gConf_add_fields(grid);
    }
    else
    {
      Error0(NO_OPTION);
    }
  }
  /* one can for instance add various params here for example:
  // else if (Pcmps("star_type","BS+WD")) */
  else
    Error0(NO_OPTION);
  
  return EXIT_SUCCESS;
}

