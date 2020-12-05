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
// file: pointer to where writing properties
// pr_screen: if 1, it ALSO prints in standard output, 
//            otherwise only in file. */
void 
bh_print_properties
  (Physics_T *const phys,
  FILE *const file,
  const int pr_screen)
{
  if (!phys || !file)
    return;
   
  AssureType (phys->ctype == BH);
  
  PR_PROPERTY_IN_FILE_d("center_x", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("center_y", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("center_z", file, pr_screen)
  
  PR_PROPERTY_IN_FILE_d("max_radius", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("min_radius", file, pr_screen)
 
  PR_PROPERTY_IN_FILE_d("irreducible_mass", file, pr_screen)
  /*PR_PROPERTY_IN_FILE_d("irreducible_mass_current", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("ADM_mass", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Kommar_mass", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Christodoulou_mass", file, pr_screen)
  
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
  PR_PROPERTY_IN_FILE_d("Rcenter_z", file, pr_screen)

  PR_PROPERTY_IN_FILE_d("Px_ADM", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Py_ADM", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Pz_ADM", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Jx_ADM", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Jy_ADM", file, pr_screen)
  PR_PROPERTY_IN_FILE_d("Jz_ADM", file, pr_screen)*/
 
  fprintf(file,"\n");
}

