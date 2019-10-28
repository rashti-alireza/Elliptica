/*
// Alireza Rashti
// October 2019
*/

#include "sbh_study_sbh.h"

/* study and analyze initial data */
void sbh_study_initial_data(Grid_T *const grid)
{
  pr_clock();
  pr_line_custom('=');
  printf("{ Studying Initial Data for Binary BH and NS ...\n");

  /* print fields */
  const char *path_par = GetParameterS_E("iteration_output");
  char *folder         = make_directory(path_par,"output_4d");
  sbh_print_fields(grid,(unsigned)GetParameterI_E("iteration_number"),folder);
  free(folder);
  
  /* convergence test */
  sbh_convergence_test(grid);
  
  printf("} Studying Initial Data for Binary BH and NS ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
}

/* printing fields determined in parameter file */
void sbh_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder)
{
  pr_line_custom('=');
  printf("{ Printing Specified Fields for Binary BH and NS ...\n");

  Pr_Field_T *pr  = init_PrField(grid);
  pr->folder = folder;
  pr->par    = "print_fields_4d";
  pr->cycle  = (int)iteration;
  pr_fields(pr);
  free_PrField(pr);
  
  printf("} Printing Specified Fields for Binary BH and NS ==> Done.\n");
  pr_line_custom('=');
}

/* some convergence tests */
void sbh_convergence_test(const Grid_T *const grid)
{
  const unsigned Nf = 3;
  const char *const Beta_nu[3]  = {"Beta_U0","Beta_U1","Beta_U2"};/* numeric */
  const char *const Beta_an[3]  = {"_Beta_U0","_Beta_U1","_Beta_U2"};/* analytic */
  unsigned f;
  
  for (f = 0; f < Nf; f++)
  {
    pr_field_difference(grid,Beta_an[f],Beta_nu[f]);
    analytic_numeric_convergence_test(grid,Beta_an[f],Beta_nu[f]);
  }
}
