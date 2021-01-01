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
// params: print these parameters (comma separated)
// file: pointer to where writing properties
// pr_screen: if 1, it ALSO prints in standard output, 
//            otherwise only in file. */
void 
sys_print_properties
  (Physics_T *const phys,
  const char *const params,
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
  
  
  char **param = read_separated_items_in_string(params,',');
  Uint p;
  
  p = 0;
  while(param[p])
  {
    PR_PROPERTY_IN_FILE_d(param[p], file, pr_screen);
    p++;
  }
  
  fprintf(file,"\n");
  
  free_2d(param);
}

