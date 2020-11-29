/*
// Alireza Rashti
// November 2020
*/

/* physical system general affairs. one can add new different 
// function readily by adding new parameter and the name of 
// the function as shown. */

#include "sys_main.h"

/* update stress energy tensor */
int sys_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case TUNE_SYS_P_ADM:
      ret = tune_system_ADM_momenta(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* tune adm momenta */
static int tune_system_ADM_momenta(Physics_T *const phys)
{
  FUNC_TIC
  
  int ret = EXIT_SUCCESS;
  
  ret = sys_tune_ADM_momenta(phys);
  
  FUNC_TOC
  return ret;
}

