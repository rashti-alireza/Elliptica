/*
// Alireza Rashti
// June 2018
*/

#include "lalpace_inhom.h"

int Laplace_Inhom(void)
{
  Grid_T *grid;
  
  /* print clock */
  pr_clock();
  
  grid = laplace_inhom_make_grid();// making grid
  //laplace_inhom_solve_eq(grid);// solving laplace eq
  //laplace_inhom_pr_answer(grid);// printing found answer
  //laplace_inhom_clean_up(grid);// cleaning up
  
  pr_clock();
  
  return EXIT_SUCCESS;
}
