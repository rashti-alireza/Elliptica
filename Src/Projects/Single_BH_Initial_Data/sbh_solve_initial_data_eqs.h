#include "sbh_headers.h"
#include "maths_equation_solvings_lib.h"
#include "sbh_XCTS_equations_lib.h"

void sbh_solve_initial_data_eqs(Grid_T *const grid);
static void sbh_XCTS_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq, sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq);
void sbh_SolveEqs_FieldUpdate(Patch_T *const patch,const char *const name);
int sbh_stop_criteria(Grid_T *const grid,const char *const name);
static void sbh_backtrack(Grid_T *const grid,const char *const name);


                          
