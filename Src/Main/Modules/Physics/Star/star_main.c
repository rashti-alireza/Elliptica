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
        
        case NS_ADD_PARAMS:
          ret = star_NS_add_params(phys);
        break;
        
        case NS_ADD_FIELDS:
          ret = star_NS_add_fields(phys);
        break;
        
        //case NS_START:
          //ret = star_NS_start_off(phys);
        //break;
        
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

