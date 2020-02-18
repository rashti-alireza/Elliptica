#include "core_lib.h"
#include "maths_equation_solvings_lib.h"
#include "memory_managing_lib.h"
#include "laplace_inhom_equations_lib.h"
#include "fields_lib.h"

int Laplace_Inhom_solve_eq(Grid_T *const grid);
int Laplace_Inhom_initial_data_alpha(Grid_T *const grid);
void Laplace_Inhom_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq);
