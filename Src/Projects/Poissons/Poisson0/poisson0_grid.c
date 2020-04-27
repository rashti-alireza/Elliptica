/*
// Alireza Rashti
// June 2018
*/

#include "poisson0_grid.h"

/* making grid for poisson0 project.
// ->return value: made grid.
*/
Grid_T *poisson0_make_grid(void)
{
  Grid_T *grid = alloc_grid();/* adding a new grid */
  grid_characteristics_example(grid);
  make_patches(grid);/* making patch(es) to cover the grid */
  realize_geometry(grid);/* realizing the geometry of whole grid
                     // including the way patches have been sewed,
                     // normal to the boundary, 
                     // outer-boundary, inner boundary and etc. */
  return grid;
}

