#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "coordinates_lib.h"
#include "memory_managing_lib.h"
#include "maths_calculus_lib.h"
#include "maths_approximation_lib.h"
#include "maths_matrix_solvers_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_general_lib.h"
#include "maths_linear_algebra_lib.h"
#include "physics_EoS_lib.h"

#define FT_OpenMP(x) _Pragma(#x)

int Fundamental_Tests(void);
Grid_T *fundamental_tests_make_grid(void);
void fundamental_tests_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq);
int fundamental_test_initial_data_alpha(Grid_T *const grid);
