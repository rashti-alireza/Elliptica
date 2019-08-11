/*
// Alireza Rashti
// August 2019
*/

#include "bbn_solve_initial_data_eqs.h"

/* solving initial data equations for the given grid */
void bbn_solve_initial_data_eqs(Grid_T *const grid)
{
  UNUSED(grid);
  
  pr_line_custom('='); 
  printf("Solving initial data equations for Binary BH and NS ...\n");
  
  sEquation_T **field_eq/* field equation */,
              **bc_eq/* B.C. for the field */,
              **jacobian_field_eq/* jacobian for field equation */,
              **jacobian_bc_eq/* jacobian for B.C. */;

  /* filling db equations of XCTS */
  bbn_XCTS_fill_db_eqs(&field_eq,&bc_eq,&jacobian_field_eq,&jacobian_bc_eq);
  
  /* populating solution managing */
  initialize_solving_man(grid,field_eq,bc_eq,jacobian_field_eq,jacobian_bc_eq);
  
  /* solving equation(s) */
  solve_eqs(grid);

  /* free data base of equations */
  free_db_eqs(field_eq);
  free_db_eqs(bc_eq);
  free_db_eqs(jacobian_field_eq);
  free_db_eqs(jacobian_bc_eq);
  
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
  
  /* eta equations */
  add_eq(field_eq,eq_eta,"eq_eta");
  add_eq(bc_eq   ,bc_eta,"bc_eta");
  add_eq(jacobian_field_eq,jacobian_eq_eta,"jacobian_eq_eta");
  add_eq(jacobian_bc_eq   ,jacobian_bc_eta,"jacobian_bc_eta");
  
  /* Beta_U0 equations */
  add_eq(field_eq,eq_Beta_U0,"eq_Beta_U0");
  add_eq(bc_eq   ,bc_Beta_U0,"bc_Beta_U0");
  add_eq(jacobian_field_eq,jacobian_eq_Beta_U0,"jacobian_eq_Beta_U0");
  add_eq(jacobian_bc_eq   ,jacobian_bc_Beta_U0,"jacobian_bc_Beta_U0");
  
  /* Beta_U1 equations */
  add_eq(field_eq,eq_Beta_U1,"eq_Beta_U1");
  add_eq(bc_eq   ,bc_Beta_U1,"bc_Beta_U1");
  add_eq(jacobian_field_eq,jacobian_eq_Beta_U1,"jacobian_eq_Beta_U1");
  add_eq(jacobian_bc_eq   ,jacobian_bc_Beta_U1,"jacobian_bc_Beta_U1");
  
  /* Beta_U2 equations */
  add_eq(field_eq,eq_Beta_U2,"eq_Beta_U2");
  add_eq(bc_eq   ,bc_Beta_U2,"bc_Beta_U2");
  add_eq(jacobian_field_eq,jacobian_eq_Beta_U2,"jacobian_eq_Beta_U2");
  add_eq(jacobian_bc_eq   ,jacobian_bc_Beta_U2,"jacobian_bc_Beta_U2");
  
}
