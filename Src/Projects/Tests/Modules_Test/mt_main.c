/*
// Alireza Rashti
// August 2018
*/

#include "mt_main.h"

/* make sure different routines and algorithms are properly works.
// tests will be done according to input file.
// ->return value: EXIT_SUCCESS
*/
int Modules_Test(void *vp)
{
  /* making output directory for this project */
  char folder[STR_LEN_MAX] = {'\0'};
  char *outdir = 0;
  sprintf(folder,"%s",Pgets("parameter_file_name_stem"));
  outdir = make_directory(Pgets("relative_root_path"),folder);
  add_parameter("output_directory_path",outdir);
  free(outdir);
  
  Grid_T *grid = mt_make_grid();/* making grid */
  sEquation_T **field_eq,**bc_eq,
              **jacobian_field_eq,**jacobian_bc_eq;/* data base of equations */
  
  if (strstr_i(PgetsEZ("Test_EoS"),"yes"))
  {
    test_EoS(grid);
  }
  
  if (strcmp_i(PgetsEZ("Test_Jacobian_Elements_Js_Values"),"yes"))
  {
    test_dfs_df_values(grid);
  }
  
  if (strcmp_i(PgetsEZ("Test_d(interp_a)/df"),"yes"))
  {
    test_dInterp_a_df(grid);
  }
  
  if (strstr_i(PgetsEZ("Test_RootFinders"),"yes"))
    test_root_finders(grid);
  
  if (strstr_i(PgetsEZ("Test_Derivative"),"yes"))
    derivative_tests(grid);
  
  if (strstr_i(PgetsEZ("Test_FourierTransformation"),"yes"))
    fourier_transformation_tests(grid);
    
  if (strstr_i(PgetsEZ("Test_Ylm_Transformation"),"yes"))
     Ylm_transformation_tests(grid);
    
  if (strstr_i(PgetsEZ("Test_Integration"),"yes"))
    integration_tests(grid);
  
  if (strstr_i(PgetsEZ("Test_CubedSpherical_Coordinates"),"yes"))
    test_CubedSpherical_Coordinates(grid);
  
  if (strcmp_i(PgetsEZ("Test_Jacobian_Elements_Js_Consistency"),"yes"))
  {
    const char *types[] = {"dfx_df","dfy_df",0};
    unsigned p;
    
    /* fill data base of equations */
    mt_fill_db_eqs(&field_eq,&bc_eq,&jacobian_field_eq,&jacobian_bc_eq);

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
  
  if (strcmp_i(PgetsEZ("Test_Matrix_Consistency"),"yes"))
  {
    matrix_tests();
  }
  
  if (strcmp_i(PgetsEZ("Test_Solver_Consistency"),"yes"))
  {
    solver_tests();
  }
  
  if (strstr_i(PgetsEZ("Test_Interpolation"),"yes"))
  {
    interpolation_tests(grid);
  }
  
  if (strstr_i(PgetsEZ("Test_Math_General"),"yes"))
  {
    summation_tests();
  }
  
  
  UNUSED(vp);
  return EXIT_SUCCESS;
}
