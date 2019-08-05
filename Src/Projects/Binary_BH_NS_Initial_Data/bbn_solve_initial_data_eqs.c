/*
// Alireza Rashti
// August 2019
*/

#include "bbn_solve_initial_data_eqs.h"

/* solving initial data equations for the given grid */
void bbn_solve_initial_data_eqs(Grid_T *const grid)
{
  pr_line_custom('='); 
  printf("Solving initial data equations for Binary BH and NS ...\n");
  
  sEquation_T **field_eq/* field equation */,
              **bc_eq/* B.C. for the field */,
              **jacobian_field_eq/* jacobian for field equation */,
              **jacobian_bc_eq/* jacobian for B.C. */;

  /* filling db equations of XCTS */
  bbn_XCTS_fill_db_eqs(&field_eq,&bc_eq,&jacobian_field_eq,&jacobian_bc_eq);

  printf("Solving initial data equations for Binary BH and NS ==> Done.\n");
  pr_clock();
  pr_line_custom('='); 
}

static void bbn_XCTS_fill_db_eqs(sEquation_T ***const field_eq, 
                          sEquation_T ***const bc_eq, 
                          sEquation_T ***const jacobian_field_eq,
                          sEquation_T ***const jacobian_bc_eq)
{
  /* adding field and boundary condition equations of XCTS formalism */
  *field_eq      = init_eq();
  *bc_eq         = init_eq();
  *jacobian_field_eq   = init_eq();
  *jacobian_bc_eq      = init_eq();

  /* psi equations */
  add_eq(field_eq,eq_psi,"eq_psi");
  add_eq(bc_eq   ,bc_psi,"bc_psi");
  add_eq(jacobian_field_eq,jacobian_eq_psi,"jacobian_eq_psi");
  add_eq(jacobian_bc_eq   ,jacobian_bc_psi,"jacobian_bc_psi");
}
