/*
// Alireza Rashti
// June 2019
*/

#include "bbn_main.h"

/* constructing initial data for system of binary black hole neutron star */
int Binary_BH_NS_Initial_Data(void)
{
  Grid_T *grid = 0;
  
  /* print clock */
  pr_clock();
  
  /* approximate this system with TOV star and Kerr-Shild blak hole */
  grid = bbn_initialize_fields_and_grid();
  
  /* constructing ID starting with the above approximation */
  UNUSED(grid);
  //grid = bbn_construct_initial_data(grid);
  
  /* print clock */
  pr_clock();
  
  return EXIT_SUCCESS;
}
