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
  
  /*PR_PROPERTY_IN_FILE_d("irreducible_mass_current", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("ADM_mass", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Komar_mass", file, pr_screen)
  
  PR_PROPERTY_IN_FILE_d("AH_area", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Sx_Campanelli", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Sy_Campanelli", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Sz_Campanelli", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Sx_JRP", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Sy_JRP", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Sz_JRP", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("chi_x", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("chi_y", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("chi_z", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Omega_U0", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Omega_U1", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Omega_U2", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Rcenter_x", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Rcenter_y", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Rcenter_z", file, pr_screen)*/
 
}

