#ifndef sys_LIB_H
#define sys_LIB_H

#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_EoS_lib.h"
#include "physics_observe_lib.h"
#include "physics_lib.h"
#include "physics_freedata_lib.h"
#include "fields_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_spectral_methods_lib.h"

/* parameter prefix */
#define P_ "sys_"

int sys_tune_ADM_momenta(Physics_T *const phys);

#endif

