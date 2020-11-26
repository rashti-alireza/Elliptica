/*
// Alireza Rashti
// November 2020
*/

/* black hole general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. 
// also in a project first mount must be called and then update. */

#include "bh_main.h"

/* main function to issue commands */
int bh_main(Physics_T *const phys)
{
  int ret = -1;
  
  switch (phys->cmd)
  {
    case TUNE_AH_RADIUS:
      ret = bh_tune_apparent_horizon_radius(phys);
    break;
    //case BH_FILLER:
      //ret = bh_fill_black_hole_inside(phys);
    //break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* adding default parameters and fields. */
int bh_mount(Grid_T *const grid)
{
  UNUSED(grid);
  return EXIT_SUCCESS;
}

