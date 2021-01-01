#include "eq_header.h"
#define STR_LEN (999)

/* define some external functions */
extern fFunc_stop_criteria_T(*eq_stop_criteria);/* when to stop solve */
extern fFunc_source_update_T(*eq_source_update);/* how to update source after solve */
extern fFunc_field_update_T(*eq_field_update);/* how to update field after solve */
extern fFunc_analyze_solution_T(*eq_analyze_solution);/* analyzing solution after solve */


/* global variables */
/* equation data base */
extern sEquation_T **eq_global_field_eq;/* field equation */
extern sEquation_T **eq_global_bc_eq;/* B.C. for the field */
extern sEquation_T **eq_global_jacobian_field_eq;/* jacobian for eq_global_field_eq */
extern sEquation_T **eq_global_jacobian_bc_eq;/* jacobian for eq_global_bc_eq */


void eq_solve_elliptic_equation(Physics_T *const phys);
static Grid_T **set_equation_grid(Physics_T *const phys,
                                   Solve_Equations_T *const SolveEqs);
static void free_equation_grid(Grid_T **lgrid);
static void backup_fields(Physics_T *const phys);
static void update_fields_relaxed_scheme(Physics_T *const phys);


