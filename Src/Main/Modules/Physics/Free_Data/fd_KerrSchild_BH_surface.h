#include "fd_header.h"
#include "maths_spectral_methods_lib.h"
#include "maths_equation_solvings_lib.h"
#include "fd_KerrSchild_header.h"

#undef x
#undef y
#undef z

#define KS_func_pass_args_sed KS_func_pass_args_macro

#define KS_set_args \
  struct Analytic_Func_Arg_S farg[1];\
  farg->x = x;\
  farg->y = y;\
  farg->z = z;\
  farg->X = fd_ks_X KS_func_pass_args_macro;\
  farg->Y = fd_ks_Y KS_func_pass_args_macro;\
  farg->Z = fd_ks_Z KS_func_pass_args_macro;\
  farg->R = fd_ks_R KS_func_pass_args_macro;



/* struct for BH surface root finder */
struct BH_surface_RootFinder_S
{
  double th;/* theta */
  double ph;/* phi */
  double R;/* R in non-boosted and non-rotated coords */
};


double *fd_find_KerrSchild_BH_surface(Physics_T *const phys);
static double BH_surface_root_finder_eq(void *params,const double *const r);

