/*
// Alireza Rashti
// December 2020
*/

/* equations manager, to set, to solve etc. */

#include "eq_main.h"

/* main function to issue commands */
int eq_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case EQ_SET_PARAMS:
      ret = set_equation_params(phys);
    break;
    
    case EQ_ADD_FIELDS:
      ret = add_equation_fields(phys);
    break;
    
    case EQ_SOLVE:
      ret = solve_equation(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* set default parameters. */
static int set_equation_params(Physics_T *const phys)
{
  FUNC_TIC
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* adding fields. */
static int add_equation_fields(Physics_T *const phys)
{
  FUNC_TIC
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* solve fields. */
static int solve_equation(Physics_T *const phys)
{
  FUNC_TIC
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}


