/*
// Alireza Rashti
// November 2020
*/

/* ADM general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. */

#include "adm_main.h"

/* main function to issue command */
int adm_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case ADM_ADD_PARAMS:
      ret = add_adm_params(phys);
    break;
    
    case ADM_ADD_FIELDS:
      ret = add_adm_fields(phys);
    break;
    
    case ADM_COMPUTE_CONSTRAINTS:
      ret = compute_ham_and_mom_constraints(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}


/* add default parameters */
static int add_adm_params(Physics_T *const phys)
{
  FUNC_TIC
  
  /* how to compute constraints:
  // options:
  // from_scratch: using beta,psi,eta and trK to make Kij, 
  //               then using psi, gConf to make g
  //               and then finally make R and using sources 
  //               to calucalte constraints.
  // from_identities: using AConfIJ and various other identities
  //                  to calculate constraints. */
  Pset_default(P_"constraints_method","from_scratch");
  
  
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* add fields parameters */
static int add_adm_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  assert(phys->grid);
  
  adm_add_3plus1_fields(phys->grid);
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}

/* compute constraints */
static int compute_ham_and_mom_constraints(Physics_T *const phys)
{
  FUNC_TIC
  
  if(strstr_i(Pgets(P_"constraints_method"),"from_scratch"))
    adm_compute_constraints(phys,".*","from_scratch"  ,"ham1","mom1");
    
  if(strstr_i(Pgets(P_"constraints_method"),"from_identities"))
    adm_compute_constraints(phys,".*","from_identities","ham2","mom2");
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}


