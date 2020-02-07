/*
// Alireza Rashti
// June 2018
*/

#include "projects_data_base.h"

/* add all of the projects you demand to run */
int create_db_projects(void)
{
  /* add all project you want with the syntax of example below:
  // example:
  // in projects_data_base.c for adding:
  // add_project(TOV_star,"TOV_Star","guessing initial data for a NS");
  // in projects_data_base.h for deceleration:
  // int TOV_star(void);
  // note: that the deceleration of project functions are all the same. 
  */
  
  add_project(Laplace_Inhom,"Laplace_Inhom",
    "Solving a Laplace eq with a source for testing purposes");
    
  add_project(Modules_Test,"Modules_Test",
    "Testing the modules of Elliptica");
    
  add_project(Binary_BH_NS_Initial_Data,"Binary_BH_NS_Initial_Data",
  "Constructing Initial Data for Binary Black Hole Neutron Star");
  
  add_project(TOV_star,"TOV_Star","guessing initial data for a NS");
  
  //add_project(Single_NS_Initial_Data,"Single_NS_Initial_Data",
  //"Constructing Initial Data for Single Neutron Star");
  
  //add_project(Single_BH_Initial_Data,"Single_BH_Initial_Data",
  //"Constructing Initial Data for Single Black-Hole");
  
  return EXIT_SUCCESS;
}
