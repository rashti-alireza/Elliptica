/*
// Alireza Rashti
// October 2019
*/

#include "sns_study_sns.h"

/* study and analyze initial data */
void sns_study_initial_data(Grid_T *const grid)
{
  pr_clock();
  pr_line_custom('=');
  printf("{ Studying Initial Data for Single NS ...\n");

  /* print fields */
  const char *path_par = GetParameterS_E("iteration_output");
  char *folder         = make_directory(path_par,"output_4d");
  sns_print_fields(grid,(unsigned)GetParameterI_E("iteration_number"),folder);
  free(folder);
  
  printf("} Studying Initial Data for Single NS ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
}

/* printing fields determined in parameter file */
void sns_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder)
{
  pr_line_custom('=');
  printf("{ Printing Specified Fields for Single NS ...\n");

  Pr_Field_T *pr  = init_PrField(grid);
  pr->folder = folder;
  pr->par    = "print_fields_4d";
  pr->cycle  = (int)iteration;
  pr_fields(pr);
  free_PrField(pr);
  
  printf("} Printing Specified Fields for Single NS ==> Done.\n");
  pr_line_custom('=');
}
