#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "maths_matrix_solvers_lib.h"
#include "maths_equation_solvings_lib.h"

#define FORMAT_ER_PAR "Each field(s) for Solving_Order parameter "\
                          "must be written in a curly brackets {}\n"
void *init_eq(void);
void add_eq(sEquation_T ***const data_base, fEquation_T *const eq,const char *const name);
void initialize_solving_man(Grid_T *const grid,sEquation_T **const field_eq,sEquation_T **const bc_eq,sEquation_T **const jacobian_field_eq, sEquation_T **const jacobian_bc_eq);
void enable_fields(Grid_T *const grid);
fEquation_T *get_field_eq(const char *const name, sEquation_T **const db);
fEquation_Solver_T *get_solver_method(const char *const solver);
void free_db_eqs(sEquation_T **db);

