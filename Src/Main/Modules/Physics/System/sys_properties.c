/*
// Alireza Rashti
// December 2020
*/


/* properties of black hole, spin, mass, etc */

#include "sys_properties.h"


void 
sys_print_properties
  (Physics_T *const phys,
  FILE *const file,
  const int pr_screen)
{
  if (!phys || !file)
    return;
    
  fprintf(file,"\n");
  
  PR_PROPERTY_IN_FILE_d("x_CM", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("y_CM", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("z_CM", file, pr_screen)
}

