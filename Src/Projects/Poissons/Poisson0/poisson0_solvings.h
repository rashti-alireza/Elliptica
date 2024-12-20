#include "core_lib.h"
#include "maths_equation_solvings_lib.h"
#include "poisson0_equations_lib.h"
#include "fields_lib.h"

int poisson0_solve_eq(Grid_T *const grid);
int poisson0_initial_data_alpha(Grid_T *const grid);
void poisson0_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq);
