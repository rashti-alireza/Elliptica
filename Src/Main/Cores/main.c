/*
// Alireza Rashti
// May 2018
*/

#include "main.h"

int main(int argn, char **argv, char **argv2)
{
  char *proj_name;
  Project_T proj;
  
  global_variables_init(argv[argn-1]);// initiating global variables
  read_input_file(argv[argn-1]);// reading and populating parameters
  
  projects_data_base();// add all of the desired projects 
                       // to the project data base
  
  proj_name = get_parameter_value("Initial_Data",LITERAL,0);
  proj.func = get_project(proj_name);
  
  //project_execute(proj1);
  
  //project_cleanup(proj1);
  //data_base_cleanup();
  
  return EXIT_SUCCESS;
}
