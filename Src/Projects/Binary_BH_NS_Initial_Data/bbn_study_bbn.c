/*
// Alireza Rashti
// August 2019
*/

#include "bbn_study_bbn.h"

/* study and analyze initial data */
void bbn_study_initial_data(Grid_T *const grid)
{
  pr_clock();
  pr_line_custom('=');
  printf("{ Studing Initial Data for Binary BH and NS ...\n");

  /* print fields */
  print_fields(grid);
  
  printf("} Studing Initial Data for Binary BH and NS ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
}

/* printing fields determined in parameter file */
static void print_fields(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Printing Specified Fields for Binary BH and NS ...\n");

  Pr_Field_T *pr  = init_PrField(grid);
  const char *path_par = GetParameterS("iteration_output");
  char *folder = make_directory(path_par,"output_3D");
  
  pr->folder = folder;
  pr->par    = "print_fields_4d";
  pr_fields(pr);
  free_PrField(pr);
  free(folder);
  
  printf("} Printing Specified Fields for Binary BH and NS ==> Done.\n");
  pr_line_custom('=');
}
