/*
// Alireza Rashti
// May 2018
*/

#include "main.h"

int main(int argn, char **argv, char **argv2)
{
  char *proj_name;
  Project_T *proj;
  
  global_variables_init(argv[argn-1]);// initiating global variables
  make_parameters(argv[argn-1]);// reading input file 
                                // and making parameters
  
  projects_data_base();// add all of the desired projects 
                       // to the project data base
  
  proj_name = get_parameter_value_S("Initial_Data",0);
  proj = get_project(proj_name);
  
  /* check if the parameter file is correct and this project exists */
  if (!proj)
    abortEr_s("There is no such %s project!\n",proj_name);
  
  project_execute(proj);
  
  //project_cleanup(proj1);
  //data_base_cleanup();
  
  return EXIT_SUCCESS;
}
