/*
// Alireza Rashti
// August 2018
*/

#include "poisson0_initial_data.h"

/* initial data for field alpha
// ->return value: EXIT_SUCCESS
*/
int poisson0_initial_data_alpha(Grid_T *const grid)
{
  Uint p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    double *alpha = patch->fields[Ind("alpha")]->v;
    Uint n;
    
    FOR_ALL_POINTS(n,patch)
      alpha[n] = Pow2(x_(n))+Pow2(y_(n))+Pow2(z_(n))+0.3;
  }
  
  return EXIT_SUCCESS;
}
