/*
// Alireza Rashti
// August 2019
*/

#include "bbn_study_bbn.h"

/* study and analyze initial data. whenever it is called, 
// it increments the solving_iteration_number by 1. */
void bbn_study_initial_data(Grid_T *const grid)
{
  pr_clock();
  pr_line_custom('=');
  printf("{ Studying Initial Data for Binary BH and NS ...\n");
  
  const char *const folder = GetParameterS_E("Diagnostics");
  int solving_iter         = GetParameterI("solving_iteration_number");
  
  /* calculating the constraints */
  bbn_calculate_constraints(grid);
  bbn_print_fields(grid,(unsigned)solving_iter,folder);
  
  update_parameter_integer("solving_iteration_number",solving_iter++);

  printf("} Studying Initial Data for Binary BH and NS ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* printing fields determined in parameter file */
void bbn_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder)
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
