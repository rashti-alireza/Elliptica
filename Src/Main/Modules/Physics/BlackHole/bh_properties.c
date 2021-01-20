/*
// Alireza Rashti
// December 2020
*/


/* properties of black hole, spin, mass, etc */

#include "bh_properties.h"


/* print black-hole properties.
// Note: the physics ctype must be of BH type
// arguments:
// ==========
// params: print these parameters (comma separated)
// file: pointer to where writing properties
// pr_screen: if 1, it ALSO prints in standard output, 
//            otherwise only in file. */
void 
bh_print_properties
  (Physics_T *const phys,
  const char *const params,
  FILE *const file,
  const int pr_screen)
{
  if (!phys || !file)
    return;
   
  AssureType (phys->ctype == BH);
  
  char **param = read_separated_items_in_string(params,',');
  Uint p;
  
  p = 0;
  while (param[p])
  {
    PR_PROPERTY_IN_FILE_d(param[p], file, pr_screen);
    p++;
  }
  
  free_2d(param);
  fprintf(file,"\n");
}

