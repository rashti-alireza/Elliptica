/*
// Alireza Rashti
// October 2019
*/

#include "sbh_main.h"

/* constructing initial data for system of single black hole */
int Single_BH_Initial_Data(void)
{
  /* print some description */
  pr_clock();
  pr_line_custom('=');
  printf("Constructing Initial Data for Single BH ...\n");
  
  /***********************************************************/
  /* adding some parameters that are used in different parts: */
  
  /* -> BH_Omega, the angular frequency of the horizon,
  // is a free vector that determines the spin of BH
  // and it is related to the dimensionless spin by:
  // BH_X = 4*BH_mass*BH_Omega .
  // we only use U2 component, since we assume BH only has spin 
  // in +/- of z direction (PRD 86 084033) */
  const double BH_X_U2 = GetParameterD_E("BH_X_U2");
  const double BH_mass = GetParameterD_E("BH_mass");
  AddParameterDoubleF("BH_Omega_U2",BH_X_U2/(4*BH_mass));
  
  /***********************************************************/
  /* the outer most main algorithm: */
  Grid_T *grid_prev = 0, 
         *grid_next = 0, 
         *grid = 0;
  const unsigned N_iter     = total_iterations_ip();
  const unsigned N_iter_par = total_iterative_parameters_ip();
  const char *path_par = GetParameterS("output_directory_path");
  char folder_name[100] = {'\0'};
  char *folder_path;
  unsigned iter;
  
  /* iterate over all parameters specified in parameter file */
  for (iter = 0; iter < N_iter; ++iter)
  {
    unsigned i;
    
    /* updating some parameters for the new round of iteration */
    update_parameter_integer("iteration_number",(int)iter);
    
    /* making a directory for this iteration and save the path */
    sprintf(folder_name,"SBH_Iteration_%02d",iter);
    folder_path = make_directory(path_par,folder_name);
    update_parameter_string("iteration_output",folder_path);
    free(folder_path);
    
    /* update the parameter accoding to the iteration number */
    update_iterative_parameter_ip(iter);
    printf("{ Iteration %u for the parameter(s) below:\n",iter);
    pr_parameters();/* printing in the folder */
    for (i = 0; i < N_iter_par; ++i)
    {
      printf("%-10s = %-10s\n",par_name_ip(i),par_value_str_ip(i));
    }
    
    /* preparing fields and grid according to the given previous grid */
    grid_next = sbh_initialize_next_grid(grid_prev);
    
    /* free previous grid completely */
    free_grid(grid_prev);
    
    /* constructing ID for the given grid */
    sbh_solve_initial_data_eqs(grid_next);
    
    /* study and analyse the new grid */
    sbh_study_initial_data(grid_next);
    
    grid_prev = grid_next;
    
    printf("} Iteration %u for the parameter(s) below is done.\n",iter);
    for (i = 0; i < N_iter_par; ++i)
    {
      printf("%-10s = %-10s\n",par_name_ip(i),par_value_str_ip(i));
    }
  }
  grid = grid_next;/* final grid */
  
  /* free grid */
  free_grid(grid);
    
  /* print some description */
  printf("Constructing Initial Data for Single BH ==> Done. :)\n");
  pr_clock();
  pr_line_custom('=');

  return EXIT_SUCCESS;
}