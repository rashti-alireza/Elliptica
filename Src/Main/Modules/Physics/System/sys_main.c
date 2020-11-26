/*
// Alireza Rashti
// November 2020
*/

/* physical system general affairs. one can add new different 
// function readily by adding new parameter and the name of 
// the function as shown.
// note: in a project first mount must be called and then update. */

#include "sys_main.h"

/* update stress energy tensor */
int sys_main(Physics_T *const phys)
{
  int ret = -1;
  
  switch (phys->cmd)
  {
    case TUNE_P_ADM:
      ret = sys_tune_ADM_momenta(phys);
    break;
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* adding default parameters and fields. */
int sys_mount(Grid_T *const grid)
{
  UNUSED(grid);
  return EXIT_SUCCESS;
}

