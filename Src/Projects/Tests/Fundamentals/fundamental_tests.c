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
    
  if (strcmp_i(GetParameterS("Test_Jacobian_Elements_Js_Values"),"yes"))
  {
    const char *const types[] = {"J_x","J_xx","J_y","J_yy","J_z","J_zz",0};
    const double start = get_time_sec();
    test_make_Js_jacobian_eq(grid,types);
    pr_spent_time(start,"Making Jacobian");
  }
  
  if (strcmp_i(GetParameterS("Test_Jacobian_Elements_Js_Consistency"),"yes"))
  {
    sEquation_T **field_eq,**bc_eq,**jacobian_eq;/* data base of equations */
    const char *types[] = {"j_x","j_y",0};
    unsigned p;
    
    /* fill data base of equations */
    fundamental_tests_fill_db_eqs(&field_eq,&bc_eq,&jacobian_eq);

    /* initializing and solving */
    populate_solution_man(grid,field_eq,bc_eq,jacobian_eq);/* populating solution_man */
    
    double start = get_time_sec();
    FT_OpenMP(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *patch = grid->patch[p];
      prepare_Js_jacobian_eq(patch,types);
    }
    pr_spent_time(start,"Making Jacobian");
    
    const char *types2[] = {"j_z","j_y",0};
    start = get_time_sec();
    FT_OpenMP(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *patch = grid->patch[p];
      prepare_Js_jacobian_eq(patch,types2);
    }
    pr_spent_time(start,"Making Jacobian");
    
    free_db_eqs(field_eq);
    free_db_eqs(bc_eq);
    free_db_eqs(jacobian_eq);
  }
  
  if (strcmp_i(GetParameterS("Test_Matrix_Consistency"),"yes"))
  {
    matrix_tests();
  }
  
  return EXIT_SUCCESS;
}
