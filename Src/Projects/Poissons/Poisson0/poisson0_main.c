/*
// Alireza Rashti
// June 2018
*/

#include "poisson0_main.h"

/* testing Poisson equations on cubed spherical and box kind grid */
int Poisson0(void *vp)
{
  Grid_T *grid;
  
  /* print clock */
  pr_clock();
  
  /* making output directory for this project */
  char folder[STR_LEN_MAX] = {'\0'};
  char *outdir = 0;
  sprintf(folder,"%s",Pgets("parameter_file_name_stem"));
  outdir = make_directory(Pgets("relative_root_path"),folder);
  add_parameter("top_directory",outdir);
  free(outdir);

  
  grid = poisson0_make_grid();/* making grid */
  poisson0_solve_eq(grid);/* solving laplace eq */
  poisson0_analyze_answer(grid);/* analyze the found answer */
  //poisson0_pr_answer(grid);/* printing found answer */
  //poisson0_clean_up(grid);/* cleaning up */
  
  /* print */
  Pr_Field_T *pr  = init_PrField(grid);
  outdir     = make_directory(Pgets("top_directory"),"output_3d");
  pr->folder = outdir;
  pr_fields(pr);
  free_PrField(pr);
  free(outdir);
  
  pr_clock();
  
  UNUSED(vp);
  return EXIT_SUCCESS;
}
