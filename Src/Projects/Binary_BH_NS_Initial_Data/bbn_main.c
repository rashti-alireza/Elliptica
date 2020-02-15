/*
// Alireza Rashti
// June 2019
*/

#include "bbn_main.h"

/* constructing initial data for system of binary black hole neutron star */
int Binary_BH_NS_Initial_Data(void)
{
  if (strcmp_i(PgetsEZ("Elliptic_Convergence_Test"),"yes"))
  {
    Elliptic_Eqs_Convergence_Test_BBN();
    return EXIT_SUCCESS;
  }
  
  /* print some description */
  pr_clock();
  pr_line_custom('=');
  printf("{ Constructing Initial Data for Binary BH and NS ...\n\n");
  
  /* setting the default parameters */
  bbn_set_default_parameters();
  
  /* the outer most iteration algorithm: */
  Grid_T *grid_prev = 0, *grid_next = 0, *grid = 0;
  unsigned iter = 0;
    
  /* main iteration loop */
  while(!Pgeti("STOP"))
  {
    printf("{ Outermost iteration %u ...\n",iter);
    
    /* update parameters and directories */
    update_parameters_and_directories(iter);
  
    /* writing checkpoints */
    bbn_write_checkpoint(grid_prev);
    
    /* preparing fields and grid according to the given previous grid */
    grid_next = bbn_initialize_next_grid(grid_prev);
    
    /* free previous grid and related parameters */
    bbn_free_grid_and_its_parameters(grid_prev);
    
    /* solve the elliptic equations for the given grid */
    bbn_solve_elliptic_eqs(grid_next);
    
    /* study and analyse the new grid */
    bbn_study_initial_data(grid_next);
    
    /* extrapolate metric fields inside the BH */
    bbn_extrapolate_metric_fields_insideBH(grid_next);
    
    grid_prev = grid_next;
    
    iter++;
    
    printf("} Outermost iteration %u ==> Done.\n",iter);
  }
  grid = grid_next;/* final grid */
    
  /* free grid */
  free_grid(grid);
  
  /* print some description */
  printf("} Constructing Initial Data for Binary BH and NS ==> Done. :)\n");
  pr_clock();
  pr_line_custom('=');

  return EXIT_SUCCESS;
}

/* convergence test for binary black hole neutron star elliptic equations.
// for each resolution it starts with the same initial guess, 
// then solve the elliptic equations. finally, one can compare the residuals
// plot or the other plots to test the convergence. */
static void Elliptic_Eqs_Convergence_Test_BBN(void)
{
  /* print some description */
  pr_clock();
  pr_line_custom('=');
  printf("Convergence Test of Elliptic Equations for Binary BH-NS ...\n\n");
  
  /* setting the default parameters */
  bbn_set_default_parameters();
  
  /* the outer most iteration algorithm: */
  const unsigned N_iter = total_iterations_ip();
  Grid_T *grid_prev = 0, *grid_next = 0, *grid = 0;
  unsigned iter;
  
  /* iterate over all parameters specified in parameter file */
  for (iter = 0; iter < N_iter; ++iter)
  {
    printf("{ Outermost Iteration %u ...\n\n",iter);
    
    /* update iterative parameters and directories */
    update_parameters_and_directories(iter);
    
    /* preparing fields and grid according to the given previous grid */
    grid_next = bbn_initialize_next_grid(grid_prev);
    
    /* free previous grid completely */
    bbn_free_grid_and_its_parameters(grid_prev);
    
    /* solve the elliptic equations for the given grid */
    bbn_solve_elliptic_eqs(grid_next);
        
    /* study and analyse the new grid */
    bbn_study_initial_data(grid_next);
    
    grid_prev = 0;
  
    printf("} Outermost Iteration %u ==> Done.\n",iter);  
  }
  grid = grid_next;/* final grid */
  
  /* free grid */
  free_grid(grid);
  
  /* print some description */
  printf("\nConvergence Test of Elliptic Equations for Binary BH-NS ==> Done. :)\n");
  pr_clock();
  pr_line_custom('=');
}

/* updating iterative parameters, STOP parameter and output directories.
// new output directory is made based on changing of resolution. */
static void update_parameters_and_directories(const unsigned iter)
{
  const unsigned N_iter_main_loop = total_iterations_ip();
  const unsigned N_iter_par       = total_iterative_parameters_ip();
  unsigned n[3];/* number of points */
  const char *path_par = Pgets("output_directory_path");
  char folder_name_next[1000] = {'\0'},
       folder_name_prev[1000] = {'\0'};
  char *folder_path,*folder_path2;
  unsigned i;
  
  /* if exceeds total iteration => stop */
  if (iter >= N_iter_main_loop || Pgeti("STOP"))
  {
    Pseti("STOP",1);
    return;
  } 
  
  /* find the previous folder name */
  n[0] = (unsigned)PgetiEZ("n_a");
  n[1] = (unsigned)PgetiEZ("n_b");
  n[2] = (unsigned)PgetiEZ("n_c");
  sprintf(folder_name_prev,"BBN_%ux%ux%u",n[0],n[1],n[2]);  
  
  /* updating some parameters for the new round of iteration */
  Pseti("iteration_number",(int)iter);
  
  /* update the parameter accoding to the iteration number */
  update_iterative_parameter_ip(iter);
  
  /* print the iterative parameters */
  if (N_iter_par)
  {
    printf("Iterative parameters:\n");
    for (i = 0; i < N_iter_par; ++i)
      printf("|--> %-30s = %-15s\n",par_name_ip(i),par_value_str_ip(i));
    printf("\n");
  }
  
  /* find the name of next folder */
  n[0] = (unsigned)PgetiEZ("n_a");
  n[1] = (unsigned)PgetiEZ("n_b");
  n[2] = (unsigned)PgetiEZ("n_c");
  
  /* this parameter helps to use some of the previous grid data */
  Pseti("did_resolution_change?",0);
  
  sprintf(folder_name_next,"BBN_%ux%ux%u",n[0],n[1],n[2]);
  /* if the resolution isn't the same or it is the first iteration */
  if (strcmp(folder_name_next,folder_name_prev) || iter == 0)/* if n is updated */
  {
    /* iteration number used in solving, reset this for each resolution */
    Pseti("solving_iteration_number",0);
    sprintf(folder_name_next,"BBN_%ux%ux%u",n[0],n[1],n[2]);
    folder_path = make_directory(path_par,folder_name_next);
    Psets("iteration_output",folder_path);
    folder_path2 = make_directory(folder_path,"Diagnostics");
    Psets("Diagnostics",folder_path2);
    free(folder_path);
    free(folder_path2);
    
    /* => resolution changed */
    Pseti("did_resolution_change?",1);
  }
}
