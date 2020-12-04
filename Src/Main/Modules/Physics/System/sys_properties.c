/*
// Alireza Rashti
// December 2020
*/


/* properties of black hole, spin, mass, etc */

#include "sys_properties.h"

/* print system properties.
// Note: the physics ctype must be of system type
// arguments:
// ==========
// file: pointer to where writing properties
// pr_screen: if 1, it ALSO prints in standard output, 
//            otherwise only in file. */
void 
sys_print_properties
  (Physics_T *const phys,
  FILE *const file,
  const int pr_screen)
{
  if (!phys || !file)
    return;
  
  AssureType (phys->ctype == BHNS ||
              phys->ctype == BHBH ||
              phys->ctype == SBH  ||
              phys->ctype == NSNS ||
              phys->ctype == SNS);
  
  
  PR_PROPERTY_IN_FILE_d("x_CM", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("y_CM", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("z_CM", file, pr_screen)
  
  fprintf(file,"\n");
}

