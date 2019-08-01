/*
// Alireza Rashti
// June 2019
*/

#include "bbn_main.h"

/* constructing initial data for system of binary black hole neutron star */
int Binary_BH_NS_Initial_Data(void)
{
  Grid_T *grid_prev = 0, 
         *grid_next = 0, 
         *grid = 0;
  const unsigned N_iter     = total_iterations_ip();
  const unsigned N_iter_par = total_iterative_parameters_ip();
  unsigned iter;
  
  /* print some description */
  pr_clock();
  pr_line_custom('=');
  printf("Constructing Initial Data for Binary BH and NS ...\n");
  
  /* iterate over all parameters specified in parameter file */
  for (iter = 0; iter < N_iter; ++iter)
  {
    unsigned i;
    /* update the parameter accoding to the iteration number */
    update_iterative_parameter_ip(iter);
    
    printf("Iteration %u for the parameter(s) below:\n",iter);
    for (i = 0; i < N_iter_par; ++i)
    {
      printf("%-10s = %-10s\n",par_name_ip(i),par_value_str_ip(i));
    }
    /* preparing fields and grid according to the given previous grid */
    grid_next = bbn_initialize_next_grid(grid_prev);
    
    /* free grid prev */
    //if (grid_next != grid_prev)/* in some case grid_new is just a pointer copy of grid_prev */
      //free_grid(grid_prev);
    
    /* constructing ID for the given grid */
    //bbn_construct_initial_data(grid_next);
    
    /* study and analyse the new grid */
    //bbn_study_initial_data(grid_next);
    
    grid_prev = grid_next;
  }
  grid = grid_next;/* final grid */
  UNUSED(grid);
  
  /* free grid */
  /* free_grid(grid); */
    
  /* print some description */
  printf("Constructing Initial Data for Binary BH and NS ==> Done. :)\n");
  pr_clock();
  pr_line_custom('=');

  return EXIT_SUCCESS;
}
