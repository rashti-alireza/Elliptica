/*
// Alireza Rashti
// June 2018
*/

#include "laplace_inhom_grid.h"

Grid_T *Laplace_Inhom_make_grid(void)
{
  Grid_T *grid = alloc_grid();// adding a new grid
  make_patches(grid);// making patch(es) to cover the grid
  realize_geometry(grid);// realizing the geometry of whole grid 
                     // including the way patches have been sewed,
                     // normal to the boundary, 
                     // outer-boundary, inner boundary and etc.
  
  return grid;
}

