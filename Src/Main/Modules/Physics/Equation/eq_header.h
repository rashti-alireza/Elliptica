#ifndef eq_LIB_H
#define eq_LIB_H

#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_EoS_lib.h"
#include "physics_lib.h"
#include "fields_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_spectral_methods_lib.h"

/* prefix of internal parameters of this project */
#define P_  "Eq_"

/* backup field name prefix */
#define P_Backup_  P_"backup_"


void eq_solve_elliptic_equation(Physics_T *const phys);


#endif

