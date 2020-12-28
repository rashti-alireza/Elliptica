/*
// Alireza Rashti
// November 2020
*/

/* Tij general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. */

#include "Tij_main.h"

/* main function to issue command */
int Tij_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case STRESS_ENERGY_UPDATE:
      AssureType(phys->ctype == NS);
      ret = update_stress_energy_tensor(phys);
    break;
    
    case STRESS_ENERGY_SET_PARAMS:
      ret = set_stress_energy_params(phys);
    break;
    
    case STRESS_ENERGY_ADD_FIELDS:
      ret = add_stress_energy_tensor_fields(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}


/* set default paramters */
static int set_stress_energy_params(Physics_T *const phys)
{
  FUNC_TIC
  
  /* fluid type:
  // options:
  // NS_ideal_fluid: like: Phys. Rev. D 100, 124046  */
  Pset_default(P_"fluid","NS_ideal_fluid");
  
  if (Pcmps(P_"fluid","NS_ideal_fluid"))
  {
    /* decomposition type:
    // options:
    // XCTS: extended conformal thin sandwich(Phys. Rev. D 100, 124046)
    //       mostly because we need lapse in enthalpy calculations.  */
    Pset_default(P_"NS_decomposition","XCTS");
    
    /* conformal metric type: 
    // options:
    // flat    => gConf  = delta_{ij},
    // general => a general gConf. */
    Pset_default(P_"NS_gConf","general");
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
  
  assert(phys->grid);
  
  if (Pcmps(P_"NS_decomposition","XCTS") &&
      Pcmps(P_"NS_gConf","general"))
  {
    Tij_NS_idealfluid_XCTS_gConf_add_fields(phys->grid);
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
  if (Pcmps(P_"fluid","NS_ideal_fluid"))
  {
    if(Pcmps(P_"NS_decomposition","XCTS") &&
       Pcmps(P_"NS_gConf","general"))
    {
      Tij_NS_idealfluid_XCTS_gConf_update(phys);
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

