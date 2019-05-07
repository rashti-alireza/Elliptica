/*
// Alireza Rashti
// August 2018
*/

#include "laplace_inhom_initial_data.h"

/* initial data for field alpha
// ->return value: EXIT_SUCCESS
*/
int Laplace_Inhom_initial_data_alpha(Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    double *alpha = patch->pool[Ind("alpha")]->v;
    unsigned n;
    
    FOR_ALL_POINTS(n,patch)
      alpha[n] = SQR(x_(n))+SQR(y_(n))+SQR(z_(n))+0.3;
  }
  
  return EXIT_SUCCESS;
}
