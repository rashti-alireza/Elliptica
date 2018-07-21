/*
// Alireza Rashti
// June 2018
*/

#include "laplace_inhom.h"

  double *sinxyz_f(Patch_T *const patch);
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
  if (strstr_i(test,"yes"))
    DerivativeTest(grid);
  
  //test
  //Patch_T *patch = grid->patch[0];
  //Field_T *df_num = add_field("Numerica_derivative","(3dim)",patch,NO);
  //df_num->v = sinxyz_f(patch);
  //make_coeffs_1d(df_num,_N0_);
  //make_coeffs_3d(df_num);
  //end
    
  pr_clock();
  
  return EXIT_SUCCESS;
}
