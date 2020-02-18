#include "core_lib.h"
#include "maths_calculus_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "manifold_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_general_lib.h"
#include "fields_lib.h"

void mt_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq);
int mt_initial_data_alpha(Grid_T *const grid);
static void *eq_alpha(void *vp1,void *vp2);
static void *bc_alpha(void *vp1,void *vp2);
static void *jacobian_eq_alpha(void *vp1,void *vp2);
static void *jacobian_bc_alpha(void *vp1,void *vp2);
