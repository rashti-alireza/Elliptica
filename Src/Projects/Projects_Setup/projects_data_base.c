/*
// Alireza Rashti
// June 2018
*/

#include "projects_data_base.h"

/* add all of the projects you demand to run */
int projects_data_base(void)
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
    
  add_project(Fundamental_Tests,"Fundamental_Tests",
    "Testing fundamentals routines and codes of abc");
    
  add_project(Binary_BH_NS_Initial_Data,"Binary_BH_NS_Initial_Data",
  "Constructing Initial Data for Binary Black Hole Neutron Star");
  
  add_project(TOV_star,"TOV_Star","guessing initial data for a NS");
  
  add_project(Single_NS_Initial_Data,"Single_NS_Initial_Data",
  "Constructing Initial Data for Single Neutron Star");
  
  return EXIT_SUCCESS;
}
