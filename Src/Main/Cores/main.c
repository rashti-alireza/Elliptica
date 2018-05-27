/*
// Alireza Rashti
// May 2018
*/

#include "main.h"

int main(int argn, char **argv, char **argv2)
{
  read_input_file(argv[argn-1]);
  
  add_projets();//add all of the desired projects to the project data base
  
  proj1 = get_project("BNS_Initial_Data");
  project_execute(proj1);
  
  project_cleanup(proj1);
  data_base_cleanup();
  
  return EXIT_SUCCESS;
}
