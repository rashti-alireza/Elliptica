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
    {
      double x = patch->node[n]->x[0];
      double y = patch->node[n]->x[1];
      double z = patch->node[n]->x[2];
      double r2 = Pow2(x)+Pow2(y)+Pow2(z);
      /* true solution + a noise */
      alpha[n] = 1./sqrt(r2+1) + 0.1*exp(-sqrt(r2));
    }
  }
  
  return EXIT_SUCCESS;
}
