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
  
  /*
  //test print
  unsigned p;
  Field_T *x_ff,*y_ff,*z_ff;
  Field_T *sinx_ff;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    x_ff = add_field("x_ff","(3dim)",patch,NO);
    y_ff = add_field("y_ff","(3dim)",patch,NO);
    z_ff = add_field("z_ff","(3dim)",patch,NO);
    sinx_ff = add_field("sinx_ff","(3dim)",patch,NO);
    
    x_ff->v = x_f(patch);
    y_ff->v = y_f(patch);
    z_ff->v = z_f(patch);
    sinx_ff->v = sinx_f(patch);
  }
  
  const char *path_par = GetParameterS("output_directory_path");
  char *folder = make_directory(path_par,"Pr_Fields_4D");
  Pr_Field_T *pr = init_PrField(grid);
  pr->folder = folder;
  pr->par = "print_fields_4d";
  pr_fields(pr);
  free(folder);
  free_PrField(pr);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    sinx_ff = patch->pool[Ind("sinx_ff")];
    x_ff    = patch->pool[Ind("x_ff")];
    y_ff    = patch->pool[Ind("y_ff")];
    z_ff    = patch->pool[Ind("z_ff")];
    remove_field(sinx_ff);
    remove_field(x_ff);
    remove_field(y_ff);
    remove_field(z_ff);
  }
  //end test print*/
    
  pr_clock();
  
  return EXIT_SUCCESS;
}
