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


int bh_tune_black_hole_radius(Physics_T *const phys);
int bh_find_black_hole_surface(Physics_T *const phys);
int bh_fill_inside_black_hole(Physics_T *const phys);

#endif

