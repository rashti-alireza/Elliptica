/*
// Alireza Rashti
// June 2018
*/

#include "laplace_inhom_grid.h"

Grid_T *laplace_inhom_make_grid(void)
{
  Grid_T *grid;
  
  grid = add_grid();// adding a new grid
  make_patches(grid);// making patch(es) to cover the grid
  make_coords(grid);// making coord sys
  realize_geometry(grid);// realizing the geometry of whole grid 
                     // including the way patches have been sewed
  
  return grid;
}

/* adding a new grid to grid struct*/
static void *add_grid(void)
{
  
}
