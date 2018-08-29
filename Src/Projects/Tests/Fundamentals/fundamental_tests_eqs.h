#include "core_lib.h"
#include "maths_solvings_lib.h"

void fundamental_tests_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_eq);
static void *eq_alpha(void *vp1,void *vp2);
static void *bc_alpha(void *vp1,void *vp2);
static void *jacobian_alpha(void *vp1,void *vp2);