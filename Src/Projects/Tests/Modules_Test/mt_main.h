#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "manifold_lib.h"
#include "maths_calculus_lib.h"
#include "maths_approximation_lib.h"
#include "maths_matrix_solvers_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_general_lib.h"
#include "maths_linear_algebra_lib.h"
#include "physics_EoS_lib.h"

#define FT_OpenMP(x) _Pragma(#x)
#define STR_LEN_MAX 400


int Modules_Test(void);
Grid_T *mt_make_grid(void);
void mt_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq);
int mt_initial_data_alpha(Grid_T *const grid);





