/*
// Alireza Rashti
// June 2018
*/

#include "laplace_inhom.h"

int Laplace_Inhom(void)
{
  Grid_T *grid;
  char *test;
  /* print clock */
  pr_clock();
  
  grid = Laplace_Inhom_make_grid();// making grid
  //Laplace_Inhom_solve_eq(grid);// solving laplace eq
  //Laplace_Inhom_pr_answer(grid);// printing found answer
  //Laplace_Inhom_clean_up(grid);// cleaning up
  
  test = get_parameter_value_S("test_derivative",0);
  if (strcmp_i(test,"yes"))
    DerivativeTest(grid);
    
  pr_clock();
  
  return EXIT_SUCCESS;
}
