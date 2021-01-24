#ifndef bh_LIB_H
#define bh_LIB_H

#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_EoS_lib.h"
#include "physics_observe_lib.h"
#include "physics_lib.h"
#include "fields_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_spectral_methods_lib.h"

/* prefix internal parameters of this project, 
// PLEASE keep it capitalized.
// NOTE: it must be different with "BH" to distinguish 
// between SBH or BHNS systems. */
#define P_ "BlackHole_"

int bh_fill_inside_black_hole(Physics_T *const phys);
double bh_bhf_poly_smoother
(const double r,const double rmax,const double rmin);

void bh_bhf_ChebTn_extrapolate
(double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill,const Uint N);


void bh_tune_BH_radius_irreducible_mass_perfect_s2(Physics_T *const phys);
void bh_find_bh_surface_perfect_s2(Physics_T *const phys);
void bh_start_off_KerrSchild_perfect_s2(Physics_T *const phys);
int bh_add_patch_inside_black_hole(Physics_T *const phys,const char *const region);

void 
bh_interpolating_fields_on_a_line
  (
  Physics_T *const phys/* physics of interest */,
  const char *const sfields_name/* comma separated fields */,
  const char *const dir/* output directory */,
  const char *const stem_g/* if stem of a metric given => test det(g) > 0 */
  );

void 
bh_print_properties
  (Physics_T *const phys,
  const char *const params,
  FILE *const file,
  const int pr_screen);  
  
void bh_start_off_IsoSchild_perfect_s2(Physics_T *const phys);
void bh_start_off_PGSchild_perfect_s2(Physics_T *const phys);
void bh_start_off_KerrSchild_general_s2(Physics_T *const phys);
void bh_find_bh_surface_KerrSchild_s2(Physics_T *const phys);
void bh_add_fields(Grid_T *const grid);
void bh_update_sConf_dsConf(Physics_T *const phys);
double bh_calculate_expansion_on_AH(Physics_T *const phys);
void bh_update_inner_BC(Physics_T *const phys);
void bh_start_off_CloseKerrSchild_perfect_s2(Physics_T *const phys);
void bh_tune_BH_chi_simple(Physics_T *const phys);

#endif

