/*
// Alireza Rashti
// August 2018
*/

#include "fundamental_tests.h"

/* make sure different routines and algorithms are properly works.
// tests will be done according to input file.
// ->return value: EXIT_SUCCESS
*/
int Fundamental_Tests(void)
{
  Grid_T *grid = fundamental_tests_make_grid();/* making grid */
  
  if (strstr_i(GetParameterS("test_derivative"),"yes"))
    DerivativeTest(grid);
    
  if (strcmp_i(GetParameterS("Test_Jacobian_Elements_Js"),"yes"))
  {
    const char *const types[] = {"J_x","J_xx","J_y","J_yy","J_z","J_zz",0};
    const double start = get_time_sec();
    test_make_Js_jacobian_eq(grid,types);
    pr_spent_time(start,"Making Jacobian");
  }
  
  return EXIT_SUCCESS;
}
