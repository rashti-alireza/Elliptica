/*
// Alireza Rashti
// Feb 2019
*/

#include "convergence_tests.h"

/* calculate L1 norm, L2 norm, L infinit of:
// (analytic field - numeric field)/max(analytic field)
// for convergence test purposes.
// NOTE: it is not completed yet!
// ->return value: EXIT_SUCCESS. */
int analytic_numeric_convergence_test(const Grid_T *const grid,const char *const f_analytic_name,const char *const f_numeric_name)
{
  Patch_T *patch;
  Field_T *f_analytic;
  Field_T *f_numeric;
  //FILE *file_patch;
  double scale;
  double *diff;
  double L1,L2,Linf;
  unsigned nn;
  unsigned p,i;
  
  pr_line_custom('=');
  printf("Convergence test ... \n\n");
  
  FOR_ALL_PATCHES(p,grid)
  {
    patch      = grid->patch[p];
    nn         = patch->nn;
    f_analytic = patch->pool[Ind(f_analytic_name)];
    f_numeric  = patch->pool[Ind(f_numeric_name)];
    scale      = L_inf(nn,f_analytic->v);
    diff       = alloc_double(nn);
    
    if (EQL(scale,0)) scale = 1;
    
    for (i = 0; i < nn; ++i)
      diff[i] = f_analytic->v[i]-f_numeric->v[i];
      
    Linf  = L_inf(nn,diff);
    Linf /= scale;
    L2    = L2_norm(nn,f_analytic->v,f_numeric->v);
    L2   /= scale;
    L1    = L1_norm(nn,f_analytic->v,f_numeric->v);
    L1   /= scale;
    
    printf("%s Linf L2 L1\n",patch->name);
    printf("%u %g %g %g\n",nn,Linf,L2,L1);
    
    free(diff);
  }
  
  printf("\nConvergence test ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
  return EXIT_SUCCESS;
}
