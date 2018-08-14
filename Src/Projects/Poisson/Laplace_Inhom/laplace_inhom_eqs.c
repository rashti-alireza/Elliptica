/*
// Alireza Rashti
// August 2018
*/

#include "laplace_inhom_eqs.h"

/* laplace equation for alpha
// \nabla^2(alpha) = sin(x y z)
*/
void laplace_inhom_eq_alpha(Field_T *const alpha,double *const F)
{
  double *alpha_xx = Partial_Derivative(alpha,"xx");
  double *alpha_yy = Partial_Derivative(alpha,"yy");
  double *alpha_zz = Partial_Derivative(alpha,"zz");
  Patch_T *const patch = alpha->patch;
  const unsigned nn = total_nodes_patch(patch);
  unsigned i;
  
  for(i = 0; i < nn; ++i)
    F[i] = alpha_xx[i]+alpha_yy[i]+alpha_zz[i] - sin(x_(i)*y_(i)*z_(i));
}
