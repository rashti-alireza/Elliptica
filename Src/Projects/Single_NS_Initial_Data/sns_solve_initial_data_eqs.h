#include "sns_headers.h"
#include "maths_equation_solvings_lib.h"
#include "sns_XCTS_equations_lib.h"

void sns_solve_initial_data_eqs(Grid_T *const grid);
static void sns_XCTS_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq, sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq);
static Grid_T *sns_phi_grid(Grid_T *const grid);
static void sns_phi_grid_CS(Grid_T *const phi_grid,Grid_T *const grid);
void sns_SolveEqs_FieldUpdate(Patch_T *const patch,const char *const name);
static void sns_free_phi_grid(Grid_T *grid);
void sns_SolveEqs_SourceUpdate(Grid_T *const grid,const char *const name);
int sns_stop_criteria(Grid_T *const grid,const char *const name);
static void sns_backtrack(Grid_T *const grid,const char *const name);


                          
