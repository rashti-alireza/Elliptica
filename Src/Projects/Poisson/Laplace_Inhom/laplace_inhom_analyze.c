/*
// Alireza Rashti
// Feb 2019
*/

#include "laplace_inhom_analyze.h"

/* analyze the found answer, it calculates the difference between 
// the calculated field and analytic value. 
// ->return value: EXIT_SUCCESS. */
int Laplace_Inhom_analyze_answer(const Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Field_T *alpha_real;
    unsigned n;
    
    alpha_real = add_field("alpha_real",0,patch,YES);
    
    FOR_ALL_POINTS(n,patch)
      alpha_real->v[n] = SQR(x_(n))+SQR(y_(n))+SQR(z_(n));
  }
  
  pr_field_difference(grid,"alpha","alpha_real");
  
  analytic_numeric_convergence_test(grid,"alpha_real","alpha");
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Field_T *f = patch->pool[Ind("alpha_real")];
    remove_field(f);
  }
  
  return EXIT_SUCCESS;
}
