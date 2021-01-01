#include "eq_header.h"
#include "eq_db_header.h"

#define STR_LEN (999)


/* define some external functions */
fFunc_stop_criteria_T(*eq_stop_criteria);/* when to stop solve */
fFunc_source_update_T(*eq_source_update);/* how to update source after solve */
fFunc_field_update_T(*eq_field_update);/* how to update field after solve */
fFunc_analyze_solution_T(*eq_analyze_solution);/* analyzing solution after solve */


/* global variables */
/* equation data base */
sEquation_T **eq_global_field_eq = 0;/* field equation */
sEquation_T **eq_global_bc_eq = 0;/* B.C. for the field */
sEquation_T **eq_global_jacobian_field_eq = 0;/* jacobian for eq_global_field_eq */
sEquation_T **eq_global_jacobian_bc_eq = 0;/* jacobian for eq_global_bc_eq */


int eq_main(Physics_T *const phys);
static int set_equation_params(Physics_T *const phys);
static int add_equation_fields(Physics_T *const phys);
static int solve_equation(Physics_T *const phys);


