/*
// Alireza Rashti
// August 2018
*/

#include "poisson0_solvings.h"

/* solving inhomogeneous Laplace equation
// ->return value: EXIT_SUCCEESS if succeeds
*/
int poisson0_solve_eq(Grid_T *const grid)
{
  Solve_Equations_T *SolveEqs = init_solve_equations(grid);/* initialization */
  sEquation_T **field_eq/* field equation */,
              **bc_eq/* B.C. for the field */,
              **jacobian_field_eq/* jacobian for field equation */,
              **jacobian_bc_eq/* jacobian for B.C. */;
  
  /* fill data base of equations */
  poisson0_fill_db_eqs(&field_eq,&bc_eq,&jacobian_field_eq,&jacobian_bc_eq);
  
  /* initializing and solving */
  initialize_solving_man(grid,field_eq,bc_eq,jacobian_field_eq,jacobian_bc_eq);/* populating solution managing */
  enable_fields(grid);/* allocating required fields in patch->pool */
  poisson0_initial_data_alpha(grid);/* initial data for field alpha */
  
  SolveEqs->solving_order = Pgets("Solving_Order");
  solve_eqs(SolveEqs);/* solving equation(s) */
  
  /* freeing */
  free_solve_equations(SolveEqs);
  free_db_eqs(field_eq);
  free_db_eqs(bc_eq);
  free_db_eqs(jacobian_field_eq);
  free_db_eqs(jacobian_bc_eq);
  
  return EXIT_SUCCESS;
}

/* adding Laplace Inhomogeneous equations to data base */
void poisson0_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_field_eq,sEquation_T ***const jacobian_bc_eq)
{
  /* adding field and boundary condition equations */
  *field_eq      = init_eq();
  *bc_eq         = init_eq();
  *jacobian_field_eq   = init_eq();
  *jacobian_bc_eq      = init_eq();

  add_eq(field_eq,eq_alpha,"eq_alpha");
  add_eq(bc_eq   ,bc_alpha,"bc_alpha");
  add_eq(jacobian_field_eq,jacobian_eq_alpha,"jacobian_eq_alpha");
  add_eq(jacobian_bc_eq   ,jacobian_bc_alpha,"jacobian_bc_alpha");
}

