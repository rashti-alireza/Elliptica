/*
// Alireza Rashti
// September 2019
*/

#include "bbn_adjustments.h"

/* adjusting different quantities and then make new grid */
Grid_T *bbn_adjust_and_make_new_grid(Grid_T *const grid)
{
  Grid_T *new_grid = 0;
  /* adjust Euler constant to fix NS baryonic mass */
  /* adjust BH center to make P_ADM zero */
  /* adjust the BH radius to acquire the desired BH mass */
  /* adjust the Omega_BH to acquire the desired BH spin */
  /* adjust y_CM using force balance equation */
  /* adjust NS surface */
  /* make new grid with new parameters */
  /* use previous grid to interpolate values of fields for new grid */
  abortEr(NO_JOB);
  UNUSED(grid);
  return new_grid;
}
    