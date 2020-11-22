/*
// Alireza Rashti
// November 2020
*/

/* star general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. 
// also in a project first mount must be called and then update. */

#include "star_main.h"

/* update stress energy tensor */
int star_tune(Physics_T *const phys)
{

  if (Pcmps("star_type","NS")           &&
      Pcmps("star_fluid","ideal_fluid") && 
      Pcmps("star_gConf","non_flat"))
  {
    switch (phys->cmd)
    {
      case EULER_CONST:
        star_idealfluid_NS_gConf_find_Euler_const(phys);
      break;
      default:
        Error0(NO_OPTION);
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
    
  return EXIT_SUCCESS;
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
  
  /* fluid type:
  // ideal_fluid: like: Phys. Rev. D 100, 124046  */
  Pset_default("star_fluid","ideal_fluid");
  
  /* conformal metric type: 
  // flat     => gConf  = delta_{ij},
  // non_flat => gConf != delta_{ij}. */
  Pset_default("star_gConf","non_flat");
  
  if (Pcmps("star_type","NS")           &&
      Pcmps("star_fluid","ideal_fluid") && 
      Pcmps("star_gConf","non_flat"))
  {
   star_NS_idealfluid_gConf_add_fields(grid);
  }
  /* one can for instance add various params here for example:
  // else if (Pcmps("star_type","BS+WD")) */
  else
    Error0(NO_OPTION);
  
  return EXIT_SUCCESS;
}

