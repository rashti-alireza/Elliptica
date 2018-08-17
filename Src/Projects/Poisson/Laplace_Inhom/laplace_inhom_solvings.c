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
  //test jacobian
  const char *const types[] = {"J_x","J_xx",0};

  test_make_jacobian_eq(grid,types);
  //end test jacobian
  return EXIT_SUCCESS;
}
