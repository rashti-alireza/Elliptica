/*
// Alireza Rashti
// June 2019
*/

#include "bbn_initial_data_main.h"

/* constructing id for bbn */
void bbn_construct_id(void)
{
  /* print some description */
  pr_clock();
  pr_line_custom('=');
  printf("{ Constructing Initial Data for Binary BH and NS ...\n");
  
  /* making output directory for this project if needed */
  if (!Pcmps("BH_NS_initialization","checkpoint_file"))
  {
    char folder[STR_LEN_MAX] = {'\0'};
    char *outdir = 0;
    sprintf(folder,"%s",Pgets("parameter_file_name_stem"));
    outdir = make_directory(Pgets("relative_root_path"),folder);
    add_parameter("output_directory_path",outdir);
    free(outdir);
  }
  else
  {
    add_parameter("output_directory_path","NOT_SPECIFIED_YET");
  }

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
    
    /* preparing fields and grid according to the given previous grid */
    grid_next = bbn_initialize_next_grid(grid_prev);
    
    /* write checkpoint before updating the params for the next grid */
    bbn_write_checkpoint(grid_next);
    
    /* free previous grid and related parameters */
    bbn_free_grid_and_its_parameters(grid_prev);
    
    /* solve the elliptic equations for the given grid */
    bbn_solve_elliptic_eqs(grid_next);
    
    /* study and analyse the new grid */
    bbn_study_initial_data(grid_next);
    
    /* extrapolate metric fields inside the BH */
    bbn_extrapolate_metric_fields_insideBH(grid_next);
    
    grid_prev = grid_next;
    
    printf("} Outermost iteration %u ==> Done.\n",iter);
    
    iter++;
  }
  grid = grid_next;/* final grid */
    
  /* free grid */
  free_grid(grid);
  
  /* print some description */
  printf("} Constructing Initial Data for Binary BH and NS ==> Done. :)\n");
  pr_clock();
  pr_line_custom('=');
}

/* convergence test for binary black hole neutron star elliptic equations.
// for each resolution it starts with the same initial guess, 
// then solve the elliptic equations. finally, one can compare the residuals
// plot or the other plots to test the convergence. */
void bbn_elliptic_eqs_convergence_test(void)
{
  /* print some description */
  pr_clock();
  pr_line_custom('=');
  printf("Convergence Test of Elliptic Equations for Binary BH-NS ...\n");
  
  /* making output directory for this project */
  char folder[STR_LEN_MAX] = {'\0'};
  char *outdir = 0;
  sprintf(folder,"%s",Pgets("parameter_file_name_stem"));
  outdir = make_directory(Pgets("relative_root_path"),folder);
  add_parameter("output_directory_path",outdir);
  free(outdir);

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
    
    /* preparing fields and grid according to the given previous grid */
    grid_next = bbn_initialize_next_grid(grid_prev);
    
    /* write checkpoint before updating the params for the next grid */
    bbn_write_checkpoint(grid_next);
    
    /* free previous grid completely */
    bbn_free_grid_and_its_parameters(grid_prev);
    
    /* solve the elliptic equations for the given grid */
    bbn_solve_elliptic_eqs(grid_next);
    
    /* study and analyse the new grid */
    bbn_study_initial_data(grid_next);
    
    grid_prev = 0;
    
    iter++;
    
    printf("} Outermost iteration %u ==> Done.\n",iter);
  }
  grid = grid_next;/* final grid */
    
  /* free grid */
  free_grid(grid);
  
  /* print some description */
  printf("Convergence Test of Elliptic Equations for Binary BH-NS ==> Done. :)\n");
  pr_clock();
  pr_line_custom('=');
}

/* updating iterative parameters, STOP parameter and output directories.
// new output directory is made based on changing of resolution. */
static void update_parameters_and_directories(const unsigned main_loop_iter)
{
  const char *const FOLDER_NANE_FORMAT = "BBN_%s_%ux%ux%u";
  const unsigned total_iters = total_iterations_ip();
  const unsigned total_ipars = total_iterative_parameters_ip();
  const unsigned iter_n = (unsigned)Pgeti("iteration_number");
  unsigned iter;/* number of iterations have been performed for the simulation */
  unsigned n[3];/* number of points */
  const char *path_par = Pgets("output_directory_path");
  char folder_name_next[1000] = {'\0'},
       folder_name_prev[1000] = {'\0'};
  char *folder_path,*folder_path2;
  const char *const parfile_name = Pgets("parameter_file_name");
  const char *const parfile_stem = Pgets("parameter_file_name_stem");
  char cp_cmd[STR_LEN_MAX];
  char *str;
  unsigned i;
  
  /* if the initialization is from checkpoint_file do nothing */
  if (Pcmps("BH_NS_initialization","checkpoint_file"))
    return;
  
  /* when starting from checkpoint, iter_n > main_loop_iter 
  // so to avoid redo the simulations we set iter to the largest. */
  if (iter_n > main_loop_iter)
  {
    iter = iter_n;
    /* updating iterative parameters for the new round of iteration */
    Pseti("iteration_number",(int)iter+1);/* +1 is crucial */
  }
  else
  {
    iter = main_loop_iter;
    /* updating iterative parameters for the new round of iteration */
    Pseti("iteration_number",(int)iter);
  }
  
  /* if exceeds total iteration => stop */
  if (iter >= total_iters)
  {
    Pseti("STOP",1);
    Pseti("iteration_number",(int)iter);
    return;
  } 
  
  /* find the previous folder name */
  n[0] = (unsigned)PgetiEZ("n_a");
  n[1] = (unsigned)PgetiEZ("n_b");
  n[2] = (unsigned)PgetiEZ("n_c");
  sprintf(folder_name_prev,FOLDER_NANE_FORMAT,parfile_stem,n[0],n[1],n[2]);  
  
  /* update the parameter accoding to the iteration number */
  update_iterative_parameter_ip(iter);
  
  /* print the iterative parameters */
  if (total_ipars && main_loop_iter)
  {
    printf("Iterative parameters:\n");
    for (i = 0; i < total_ipars; ++i)
      printf("|--> %-30s = %-15s\n",par_name_ip(i),par_value_str_ip(i));
  }
  
  /* find the name of next folder */
  n[0] = (unsigned)PgetiEZ("n_a");
  n[1] = (unsigned)PgetiEZ("n_b");
  n[2] = (unsigned)PgetiEZ("n_c");
  
  /* this parameter helps to use some of the previous grid data */
  Pseti("did_resolution_change?",0);
  
  sprintf(folder_name_next,FOLDER_NANE_FORMAT,parfile_stem,n[0],n[1],n[2]);
  /* if the resolution isn't the same or it is the first iteration */
  if (strcmp(folder_name_next,folder_name_prev) || iter == 0)/* if n is updated */
  {
    /* iteration number used in solving, reset this for each resolution */
    Pseti("solving_iteration_number",0);
    sprintf(folder_name_next,FOLDER_NANE_FORMAT,parfile_stem,n[0],n[1],n[2]);
    folder_path = make_directory(path_par,folder_name_next);
    Psets("iteration_output",folder_path);
    folder_path2 = make_directory(folder_path,"Diagnostics");
    Psets("Diagnostics",folder_path2);
    
    /* copy the parameter file into the new directory */
    str = strrchr(folder_path,'/');
    assert(str);
    str++;
    sprintf(cp_cmd,"cp %s/%s %s/%s.par",
    Pgets("relative_root_path"),
    parfile_name,folder_path,str);
    shell_command(cp_cmd);
    
    free(folder_path);
    free(folder_path2);
    
    /* => resolution changed */
    Pseti("did_resolution_change?",1);
  }
}

