/*
// Alireza Rashti
// Feb 2019
*/

#include "poisson0_analyze.h"

/* analyze the found answer, it calculates the difference between 
// the calculated field and analytic value. 
// ->return value: EXIT_SUCCESS. */
int poisson0_analyze_answer(const Grid_T *const grid)
{
  Uint p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Field_T *alpha_real;
    Uint n;
    
    alpha_real = add_field("alpha_real",0,patch,YES);
    
    FOR_ALL_POINTS(n,patch)
    {
      double x = patch->node[n]->x[0];
      double y = patch->node[n]->x[1];
      double z = patch->node[n]->x[2];

      alpha_real->v[n] = 1./sqrt(Pow2(x)+Pow2(y)+Pow2(z)+1);
    }
  }
  
  pr_field_difference(grid,"alpha","alpha_real");
  
  analytic_numeric_convergence_test(grid,"alpha_real","alpha");
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Field_T *f = patch->fields[Ind("alpha_real")];
    remove_field(f);
  }
  
  return EXIT_SUCCESS;
}
