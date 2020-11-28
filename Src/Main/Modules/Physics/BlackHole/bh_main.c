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
  int ret = EXIT_SUCCESS;
  
  assert(phys->ctype == BH);
  
  switch (phys->cmd)
  {
    case TUNE_BH_RADIUS:
      ret = bh_tune_black_hole_radius(phys);
    break;
    case FIND_BH_SURFACE:
      ret = bh_find_black_hole_surface(phys);
    break;
    case FILL_BH:
      ret = bh_fill_inside_black_hole(phys);
    break;
    case BH_GRID_INITIAL_PARAMS:
      ret = bh_set_initial_grid_parameters(phys);
    break;
    
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

