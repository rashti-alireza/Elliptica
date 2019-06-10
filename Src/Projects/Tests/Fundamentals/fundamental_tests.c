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
  sEquation_T **field_eq,**bc_eq,
              **jacobian_field_eq,**jacobian_bc_eq;/* data base of equations */
  
  if (strstr_i(GetParameterS("Test_EoS"),"yes"))
  {
    test_EoS(grid);
  }
  
  
  if (strstr_i(GetParameterS("Test_Schur_Complement"),"yes"))
  {
     /* fill data base of equations */
    fundamental_tests_fill_db_eqs(&field_eq,&bc_eq,&jacobian_field_eq,&jacobian_bc_eq);
    /* initializing and solving */
    initialize_solving_man(grid,field_eq,bc_eq,jacobian_field_eq,jacobian_bc_eq);/* populating solution_man */
    /* allocating alpha field */
    enable_fields(grid);
    /* initial data for alpha field */
    fundamental_test_initial_data_alpha(grid);
    /* testing schur complement method */
    test_solve_ddm_schur_complement(grid);
    
    free_db_eqs(field_eq);
    free_db_eqs(bc_eq);
    free_db_eqs(jacobian_field_eq);
    free_db_eqs(jacobian_bc_eq);
  }
    
  if (strcmp_i(GetParameterS("Test_Jacobian_Elements_Js_Values"),"yes"))
  {
    test_dfs_df_values(grid);
  }
  
  if (strcmp_i(GetParameterS("Test_d(interp_a)/df"),"yes"))
  {
    test_dInterp_a_df(grid);
  }
  
  if (strstr_i(GetParameterS("Test_Derivative"),"yes"))
    derivative_tests(grid);
    
  if (strstr_i(GetParameterS("Test_Integral"),"yes"))
    integration_tests(grid);
  
  if (strstr_i(GetParameterS("Test_dNi/dxj"),"yes"))
    test_dNi_dxj(grid);
  
  if (strcmp_i(GetParameterS("Test_Jacobian_Elements_Js_Consistency"),"yes"))
  {
    const char *types[] = {"dfx_df","dfy_df",0};
    unsigned p;
    
    /* fill data base of equations */
    fundamental_tests_fill_db_eqs(&field_eq,&bc_eq,&jacobian_field_eq,&jacobian_bc_eq);

    /* initializing and solving */
    initialize_solving_man(grid,field_eq,bc_eq,jacobian_field_eq,jacobian_bc_eq);/* populating solution_man */
    
    double start = get_time_sec();
    FT_OpenMP(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *patch = grid->patch[p];
      prepare_Js_jacobian_eq(patch,types);
    }
    pr_spent_time(start,"Making Jacobian");
    
    const char *types2[] = {"dfz_df","dfy_df",0};
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
    free_db_eqs(jacobian_field_eq);
    free_db_eqs(jacobian_bc_eq);
  }
  
  if (strcmp_i(GetParameterS("Test_Matrix_Consistency"),"yes"))
  {
    matrix_tests();
  }
  
  if (strcmp_i(GetParameterS("Test_Solver_Consistency"),"yes"))
  {
    solver_tests();
  }
  
  if (strstr_i(GetParameterS("Test_Interpolation"),"yes"))
  {
    interpolation_tests(grid);
  }
  
  if (strstr_i(GetParameterS("Test_Math_General"),"yes"))
  {
    summation_tests();
  }
  
  
  return EXIT_SUCCESS;
}
