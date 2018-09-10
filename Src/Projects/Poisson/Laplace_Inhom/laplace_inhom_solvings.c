/*
// Alireza Rashti
// August 2018
*/

#include "laplace_inhom_solvings.h"

/* solving inhomogeneous Laplace equation
// ->return value: EXIT_SUCCEESS if succeeds
*/
int Laplace_Inhom_solve_eq(Grid_T *const grid)
{
  sEquation_T **field_eq,**bc_eq,**jacobian_eq;/* data base of equations */
  
  /* fill data base of equations */
  Laplace_Inhom_fill_db_eqs(&field_eq,&bc_eq,&jacobian_eq);
  
  /* initializing and solving */
  populate_solution_man(grid,field_eq,bc_eq,jacobian_eq);/* populating solution managing */
  enable_fields(grid);/* allocating required fields in patch->pool */
  Laplace_Inhom_initial_data_alpha(grid);/* initial data for field alpha */
  solve_eqs(grid);/* solving equation(s) */
  
  free_db_eqs(field_eq);
  free_db_eqs(bc_eq);
  free_db_eqs(jacobian_eq);
  
  return EXIT_SUCCESS;
}
