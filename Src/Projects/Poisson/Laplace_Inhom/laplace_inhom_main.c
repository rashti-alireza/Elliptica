/*
// Alireza Rashti
// June 2018
*/

#include "laplace_inhom_main.h"

//double *sinxyz_f(Patch_T *const patch);
int Laplace_Inhom(void)
{
  Grid_T *grid;
  
  /* print clock */
  pr_clock();
  
  /* making output directory for this project */
  char folder[STR_LEN_MAX] = {'\0'};
  char *outdir = 0;
  sprintf(folder,"%s",Pgets("parameter_file_name_stem"));
  outdir = make_directory(Pgets("relative_root_path"),folder);
  add_parameter("output_directory_path",outdir);
  free(outdir);

  
  grid = Laplace_Inhom_make_grid();/* making grid */
  Laplace_Inhom_solve_eq(grid);/* solving laplace eq */
  Laplace_Inhom_analyze_answer(grid);/* analyze the found answer */
  //Laplace_Inhom_pr_answer(grid);/* printing found answer */
  //Laplace_Inhom_clean_up(grid);/* cleaning up */
  
  /* print */
  Pr_Field_T *pr  = init_PrField(grid);
  outdir     = make_directory(Pgets("output_directory_path"),"output_3d");
  pr->folder = outdir;
  pr_fields(pr);
  free_PrField(pr);
  free(outdir);
  
  pr_clock();
  
  return EXIT_SUCCESS;
}
