/*
// Alireza Rashti
// June 2019
*/

#include "bbn_main.h"

/* constructing initial data for system of binary black hole neutron star */
int Binary_BH_NS_Initial_Data(void)
{
  /* print some description */
  pr_clock();
  pr_line_custom('=');
  printf("Constructing Initial Data for Binary BH and NS ...\n");
  
  /* the outer most iteration algorithm: */
  Grid_T *grid_prev = 0, 
         *grid_next = 0, 
         *grid = 0;
  const unsigned N_iter     = total_iterations_ip();
  const unsigned N_iter_par = total_iterative_parameters_ip();
  unsigned n[3];/* number of points */
  const char *path_par = GetParameterS("output_directory_path");
  char folder_name_next[1000] = {'\0'},
       folder_name_prev[1000] = {'\0'};
  char *folder_path;
  unsigned iter;
  
  /* iterate over all parameters specified in parameter file */
  for (iter = 0; iter < N_iter; ++iter)
  {
    unsigned i;
    
    /* updating some parameters for the new round of iteration */
    update_parameter_integer("iteration_number",(int)iter);
    
    /* update the parameter accoding to the iteration number */
    update_iterative_parameter_ip(iter);
    
    /* making a directory for this iteration and save the path */
    n[0] = (unsigned)GetParameterI("n_a");
    n[1] = (unsigned)GetParameterI("n_b");
    n[2] = (unsigned)GetParameterI("n_c");
    
    sprintf(folder_name_next,"BBN_%ux%ux%u",n[0],n[1],n[2]);
    if (strcmp(folder_name_next,folder_name_prev))/* if n is updated */
    {
      /* iteration number used in solving, reset this for each resolution */
      update_parameter_integer("solving_iteration_number",0);
      sprintf(folder_name_next,"BBN_%ux%ux%u",n[0],n[1],n[2]);
      sprintf(folder_name_prev,"BBN_%ux%ux%u",n[0],n[1],n[2]);
      folder_path = make_directory(path_par,folder_name_next);
      update_parameter_string("iteration_output",folder_path);
      free(folder_path);
    }
    
    printf("{ Iteration %u for the parameter(s) below:\n",iter);
    pr_parameters();/* printing in the folder */
    for (i = 0; i < N_iter_par; ++i)
    {
      printf("%-10s = %-10s\n",par_name_ip(i),par_value_str_ip(i));
    }
    
    /* preparing fields and grid according to the given previous grid */
    grid_next = bbn_initialize_next_grid(grid_prev);
    
    /* free previous grid completely */
    free_grid(grid_prev);
    
    /* constructing ID for the given grid */
    bbn_solve_initial_data_eqs(grid_next);
    
    /* calculating the constraints */
    bbn_calculate_constraints(grid_next);
    
    /* study and analyse the new grid */
    bbn_study_initial_data(grid_next);
    
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
  printf("Constructing Initial Data for Binary BH and NS ==> Done. :)\n");
  pr_clock();
  pr_line_custom('=');

  return EXIT_SUCCESS;
}
