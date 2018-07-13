/*
// Alireza Rashti
// July 2018
*/

#include "tests.h"
#define ArgM(a) a,#a
/* testing:
// ========
//
// given a set of analytic functions, it expands them in these basis and
// take its various derivatives ans finally 
// compares the analytic results with numeric results.
// note: the results will be printed accordingly in 
// "ExpantionTests_Derivatives" folder.
// note: only those patches that use basis will be compared.
// ->return value: EXIT_SUCCESS;
*/
int DerivativeTest(const Grid_T *const grid)
{
  sFunc_Grid2Double_T **func;
  unsigned fi;
  
  init_func_Grid2Double(&func);
  /* below is all the available analytic functions with their derivatives. 
  // one can add more function, with the same style below.
  // note: derivative must be shown by underline '_' and direction of
  // derivative as it has been shown. it is "recommended" to use the same
  // notation for naming of functions.
  // note: functions are defined in Analytic folder in Maths.
  */
  add_func_Grid2Double(ArgM(c_f));
  add_func_Grid2Double(ArgM(x_f));
  add_func_Grid2Double(ArgM(y_f));
  add_func_Grid2Double(ArgM(z_f));
  add_func_Grid2Double(ArgM(poly3_f));
  add_func_Grid2Double(ArgM(r_f,));
  add_func_Grid2Double(ArgM(sinxyz_f,));
  add_func_Grid2Double(ArgM(mix2_f));
  add_func_Grid2Double(ArgM(c_f_x));
  add_func_Grid2Double(ArgM(c_f_y));
  add_func_Grid2Double(ArgM(c_f_z));
  add_func_Grid2Double(ArgM(x_f_x));
  add_func_Grid2Double(ArgM(x_f_y));
  add_func_Grid2Double(ArgM(x_f_z));
  add_func_Grid2Double(ArgM(x_f_xx));
  add_func_Grid2Double(ArgM(y_f_x));
  add_func_Grid2Double(ArgM(y_f_y));
  add_func_Grid2Double(ArgM(y_f_z));
  add_func_Grid2Double(ArgM(y_f_yy));
  add_func_Grid2Double(ArgM(z_f_x));
  add_func_Grid2Double(ArgM(z_f_y));
  add_func_Grid2Double(ArgM(z_f_z));
  add_func_Grid2Double(ArgM(z_f_zz));
  add_func_Grid2Double(ArgM(poly3_f_x));
  add_func_Grid2Double(ArgM(poly3_f_y));
  add_func_Grid2Double(ArgM(poly3_f_z));
  add_func_Grid2Double(ArgM(poly3_f_xx));
  add_func_Grid2Double(ArgM(poly3_f_yy));
  add_func_Grid2Double(ArgM(poly3_f_zz));
  add_func_Grid2Double(ArgM(poly3_f_xy));
  add_func_Grid2Double(ArgM(poly3_f_xz));
  add_func_Grid2Double(ArgM(poly3_f_yz));
  add_func_Grid2Double(ArgM(poly3_f_xyz));
  add_func_Grid2Double(ArgM(r_f_x));
  add_func_Grid2Double(ArgM(r_f_y));
  add_func_Grid2Double(ArgM(r_f_z));
  add_func_Grid2Double(ArgM(r_f_xx));
  add_func_Grid2Double(ArgM(r_f_yy));
  add_func_Grid2Double(ArgM(r_f_zz));
  add_func_Grid2Double(ArgM(r_f_xy));
  add_func_Grid2Double(ArgM(r_f_xz));
  add_func_Grid2Double(ArgM(r_f_yz));
  add_func_Grid2Double(ArgM(r_f_xyz));
  add_func_Grid2Double(ArgM(sinxyz_f_x));
  add_func_Grid2Double(ArgM(sinxyz_f_y));
  add_func_Grid2Double(ArgM(sinxyz_f_z));
  add_func_Grid2Double(ArgM(sinxyz_f_xx));
  add_func_Grid2Double(ArgM(sinxyz_f_yy));
  add_func_Grid2Double(ArgM(sinxyz_f_zz));
  add_func_Grid2Double(ArgM(sinxyz_f_xy));
  add_func_Grid2Double(ArgM(sinxyz_f_xz));
  add_func_Grid2Double(ArgM(sinxyz_f_yz));
  add_func_Grid2Double(ArgM(sinxyz_f_xyz));
  add_func_Grid2Double(ArgM(mix2_f_x));
  add_func_Grid2Double(ArgM(mix2_f_y));
  add_func_Grid2Double(ArgM(mix2_f_z));
  add_func_Grid2Double(ArgM(mix2_f_xx));
  add_func_Grid2Double(ArgM(mix2_f_yy));
  add_func_Grid2Double(ArgM(mix2_f_zz));
  add_func_Grid2Double(ArgM(mix2_f_xy));
  add_func_Grid2Double(ArgM(mix2_f_xz));
  add_func_Grid2Double(ArgM(mix2_f_yz));
  add_func_Grid2Double(ArgM(mix2_f_xyz));
  
  FOR_ALL(fi,func)
  {
    sFunc_Grid2Double_T *F[N_FUNC];
    double *anal[N_FUNC];/* analytic */
    double *numc[N_FUNC];/* numeric */
    enum FUNC_E e;
    
    if (func[fi]->done == 1) continue;
    
    init_anal(anal);/* make them to point to 0 */
    init_numc(numc);/* make them to point to 0 */
    read_F(F,func);/* get a function name and find all of it derivative in data base */
    
    /* if any of these functions below has not been defined,anal pointer,
    // points to NULL; moreover, one can add more types of derivative
    // by the same format, just remember to change the enum consequently.
    */
    
    anal[FUNC]     = F[FUNC]->func(grid);
    anal[FUNC_x]   = F[FUNC_x]->func(grid);
    anal[FUNC_y]   = F[FUNC_y]->func(grid);
    anal[FUNC_z]   = F[FUNC_z]->func(grid);
    anal[FUNC_xx]  = F[FUNC_xx]->func(grid);
    anal[FUNC_yy]  = F[FUNC_yy]->func(grid);
    anal[FUNC_zz]  = F[FUNC_zz]->func(grid);
    anal[FUNC_xy]  = F[FUNC_xy]->func(grid);
    anal[FUNC_xz]  = F[FUNC_xz]->func(grid);
    anal[FUNC_yz]  = F[FUNC_yz]->func(grid);
    anal[FUNC_xyz] = F[FUNC_xyz]->func(grid);
    
    assert(anal[FUNC]);
    
    for (e = FUNC_x; e < N_FUNC; ++e)
    {
      /* if anal is define and so not null */
      if (anal[e])
      {
        numc[e] = derivative(anal[FUNC],e);
        analyze(numc[e],anal[e],grid);
      }
    }
    
    free_anal(anal);
    free_numc(numc);
  }
  
  return EXIT_SUCCESS;
}

/* testing: Chebyshev first kind basis. */
static void ChebyshevFirstKindBasis_DerivativeTest(const Patch_T *const patch)
{
}