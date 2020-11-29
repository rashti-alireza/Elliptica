/*
// Alireza Rashti
// November 2020
*/

/* Tij general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. 
// also in a project first mount must be called and then update. */

#include "Tij_main.h"

/* main function to issue command */
int Tij_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case UPDATE_STRESS_ENERGY:
      ret = update_stress_energy_tensor(phys);
    break;
    
    case STRESS_ENERGY_ADD_PARAMS:
      ret = add_stress_energy_parameters(phys);
    break;
    
    case STRESS_ENERGY_ADD_FIELDS:
      ret = add_stress_energy_tensor_fields(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}


/* add default parameters */
static int add_stress_energy_parameters(Physics_T *const phys)
{
  FUNC_TIC
  
  /* fluid type:
  // options:
  // NS_ideal_fluid: like: Phys. Rev. D 100, 124046  */
  Pset_default("Tij_fluid","NS_ideal_fluid");
  
  if (Pcmps("Tij_fluid","NS_ideal_fluid"))
  {
    /* decomposition type:
    // options:
    // CTS (conformal thin sandwich): like: Phys. Rev. D 100, 124046  */
    Pset_default("Tij_NS_decomposition","CTS");
    
    /* conformal metric type: 
    // options:
    // flat    => gConf  = delta_{ij},
    // general => a general gConf. */
    Pset_default("Tij_NS_gConf","general");
  }
  else
    Error0(NO_OPTION);
    
  UNUSED(phys);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* add fields parameters */
static int add_stress_energy_tensor_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  if (Pcmps("Tij_NS_decomposition","CTS") &&
      Pcmps("Tij_NS_gConf","general"))
  {
    Tij_NS_idealfluid_CTS_gConf_add_fields(phys->grid);
  }
  else
    Error0(NO_OPTION);
    
  FUNC_TOC
  return EXIT_SUCCESS; 
}

/* update stress energy tensor */
static int update_stress_energy_tensor(Physics_T *const phys)
{
  FUNC_TIC
  if (Pcmps("Tij_fluid","NS_ideal_fluid"))
  {
    if(Pcmps("Tij_NS_decomposition","CTS") &&
       Pcmps("Tij_NS_gConf","general"))
    {
      Tij_NS_idealfluid_CTS_gConf_update(phys);
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

