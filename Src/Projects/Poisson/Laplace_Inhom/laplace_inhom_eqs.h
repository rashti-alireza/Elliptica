#include "core_lib.h"
#include "maths_calculus_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "coordinates_lib.h"
#include "maths_solvings_lib.h"
#include "maths_general_lib.h"
#include "memory_managing_lib.h"

void Laplace_Inhom_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_eq);
static void *eq_alpha(void *vp1,void *vp2);
static void *bc_alpha(void *vp1,void *vp2);
static void *jacobian_alpha(void *vp1,void *vp2);
