#include "eq_header.h"
#include "eq_db_header.h"

/* equation data base */
static sEquation_T **field_eq = 0;/* field equation */
static sEquation_T **bc_eq = 0;/* B.C. for the field */
static sEquation_T **jacobian_field_eq = 0;/* jacobian for field_eq */
static sEquation_T **jacobian_bc_eq = 0;/* jacobian for bc_eq */


int eq_main(Physics_T *const phys);
static int set_equation_params(Physics_T *const phys);
static int add_equation_fields(Physics_T *const phys);
static int solve_equation(Physics_T *const phys);



