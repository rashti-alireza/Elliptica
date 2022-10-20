/*
// Alireza Rashti
// January 2021
*/


/* properties of star, spin, mass, etc */

#include "star_properties.h"


/* print star properties.
// arguments:
// ==========
// params: print these parameters (comma separated)
// file: pointer to where writing properties
// pr_screen: if 1, it ALSO prints in standard output, 
//            otherwise only in file. */
void 
star_print_properties
  (Physics_T *const phys,
  const char *const params,
  FILE *const file,
  const int pr_screen)
{
  if (!phys || !file)
    return;
   
  AssureType (phys->ctype == NS);
  
  char **param = read_separated_items_in_string(params,',');
  Uint p;
  
  p = 0;
  while (param[p])
  {
    PR_PROPERTY_IN_FILE_d(param[p], file, pr_screen);
    p++;
  }
  free_2d(param);
  
  PR_PROPERTY_IN_FILE_s("EoS_description",file, pr_screen)
  PR_PROPERTY_IN_FILE_s("EoS_type",file, pr_screen)
  PR_PROPERTY_IN_FILE_s("EoS_unit",file, pr_screen)
  PR_PROPERTY_IN_FILE_s("EoS_K0",file, pr_screen)
  PR_PROPERTY_IN_FILE_s("EoS_Gamma",file, pr_screen)
  PR_PROPERTY_IN_FILE_s("EoS_rho0_th",file, pr_screen)
  
  fprintf(file,"\n");
}

