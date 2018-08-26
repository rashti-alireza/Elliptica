#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "memory_managing_lib.h"
#include "maths_solvers_lib.h"

#define FORMAT_ER_PAR "Each field(s) for Solving_Order parameter "\
                          "must be written in a curly brackets {}\n"
void *init_eq(void);
void add_eq(sEquation_T ***const data_base, fEquation_T *const eq,const char *const name);
void populate_solution_man(Grid_T *const grid,sEquation_T **const field_eq,sEquation_T **const bc_eq,sEquation_T **const jacobian_eq);
void enable_fields(Grid_T *const grid);
fEquation_T *get_field_eq(const char *const name, sEquation_T **const db);
fEquation_Solver_T *get_solver_method(const char *const solver);
static void fill_solve(Grid_T *const grid,char **const group,const unsigned ng,sEquation_T **const field_eq,sEquation_T **const bc_eq,sEquation_T **const jacobian_eq);
static Field_T *prepare_field(const char *const name,const char *const attr,Patch_T *const patch);

