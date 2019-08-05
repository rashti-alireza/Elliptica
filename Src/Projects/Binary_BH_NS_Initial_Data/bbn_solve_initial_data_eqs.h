#include "core_lib.h"
#include "maths_calculus_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "coordinates_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_general_lib.h"
#include "memory_managing_lib.h"
#include "bbn_XCTS_equations_lib.h"

void bbn_solve_initial_data_eqs(Grid_T *const grid);
static void bbn_XCTS_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq, sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq);

                          
