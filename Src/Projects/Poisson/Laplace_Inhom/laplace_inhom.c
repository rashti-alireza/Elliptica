/*
// Alireza Rashti
// June 2018
*/

#include "laplace_inhom.h"

int Laplace_Inhom(void)
{
  Grid_T *grid;
  
  /* print clock */
  pr_clock();
  
  grid = Laplace_Inhom_make_grid();// making grid
  //Laplace_Inhom_solve_eq(grid);// solving laplace eq
  //Laplace_Inhom_pr_answer(grid);// printing found answer
  //Laplace_Inhom_clean_up(grid);// cleaning up
  
  pr_clock();
  
  return EXIT_SUCCESS;
}
