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


int bh_fill_inside_black_hole(Physics_T *const phys);
double bh_bhf_poly_smoother
(const double r,const double rmax,const double rmin);

void bh_bhf_ChebTn_extrapolate
(double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill,const unsigned N);


void bh_tune_BH_radius_irreducible_mass_perfect_s2(Physics_T *const phys);
void bh_find_bh_surface_perfect_s2(Physics_T *const phys);
void bh_start_off_KerrSchild_perfect_s2(Physics_T *const phys);



#endif

