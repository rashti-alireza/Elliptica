/*
// Alireza Rashti
// June 2018
*/

#include "laplace_inhom.h"

//double *sinxyz_f(Patch_T *const patch);
int Laplace_Inhom(void)
{
  Grid_T *grid;
  
  /* print clock */
  pr_clock();
  
  grid = Laplace_Inhom_make_grid();/* making grid */
  Laplace_Inhom_solve_eq(grid);/* solving laplace eq */
  Laplace_Inhom_analyze_answer(grid);/* analyze the found answer */
  //Laplace_Inhom_pr_answer(grid);/* printing found answer */
  //Laplace_Inhom_clean_up(grid);/* cleaning up */
  
  /* print */
  Pr_Field_T *pr  = init_PrField(grid);
  const char *path_par = Pgets_E("output_directory_path");
  char *folder = make_directory(path_par,"output_4d");
  pr->folder = folder;
  pr->par    = "print_fields_4d";
  pr_fields(pr);
  free_PrField(pr);
  free(folder);
  
  pr_clock();
  
  return EXIT_SUCCESS;
}
