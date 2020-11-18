/*
// Alireza Rashti
// November 2020
*/

/* general function used in various Projects */

#include "project_functions.h"


/* updating iterative parameters, STOP parameter and output directories.
// new output directory is made based on changing of resolution. */
void 
update_parameters_and_directories
  (
   const unsigned main_loop_iter,
   const char *const dir_name_format/* eg: "BBN_%s_%ux%ux%u" */
  )
{
  assert(dir_name_format);
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
  if (Pcmps("project_initialization","checkpoint_file"))
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
  sprintf(folder_name_prev,dir_name_format,parfile_stem,n[0],n[1],n[2]);  
  
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
  
  sprintf(folder_name_next,dir_name_format,parfile_stem,n[0],n[1],n[2]);
  /* if the resolution isn't the same or it is the first iteration */
  if (strcmp(folder_name_next,folder_name_prev) || iter == 0)/* if n is updated */
  {
    /* iteration number used in solving, reset this for each resolution */
    Pseti("solving_iteration_number",0);
    sprintf(folder_name_next,dir_name_format,parfile_stem,n[0],n[1],n[2]);
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
