/*
// Alireza Rashti
// November 2020
*/

/* frd general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. */

#include "frd_main.h"

/* main function to issue command */
int frd_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case FREE_DATA_ADD_PARAMS:
      ret = add_free_data_params(phys);
    break;
    
    case FREE_DATA_ADD_FIELDS:
      ret = add_free_data_fields(phys);
    break;
    
    //case FREE_DATA_POPULATE:
      //ret = populate_free_data(phys);
    //break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}


/* add default parameters */
static int add_free_data_params(Physics_T *const phys)
{
  FUNC_TIC
  
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* add fields parameters */
static int add_free_data_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS; 
}

