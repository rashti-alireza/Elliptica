#include "bbn_headers.h"


#define IJ(i,j) ((j)+Nphi*(i))

/* getting field name (xName) and put the interpolant into ixName */
#define INTERPOLATE_macro(xName) \
  { \
    Interpolation_T *interp_##xName = init_interpolation(); \
    interp_##xName->field = patch->pool[Ind(#xName)]; \
    interp_##xName->XY_dir_flag  = 1; \
    interp_##xName->X            = X[0]; \
    interp_##xName->Y            = X[1];\ 
    interp_##xName->K            = type_flg == 1/* NS? */? patch->n[2]-1 : 0; \
    plan_interpolation(interp_##xName); \
    i##xName = execute_interpolation(interp_##xName); \
    free_interpolation(interp_##xName); \
  }

static void find_XYZ_and_patch_of_theta_phi_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid,const char *const type);
void bbn_populate_2d_induced_metric_S2_theta_phi(double *const h_D0D0,double *const h_D0D1,double *const h_D1D1,const double *const g_D0D0,const double *const g_D0D1,const double *const g_D0D2,const double *const g_D1D1,const double *const g_D1D2,const double *const g_D2D2,const double theta,const double phi);
void bbn_populate_3d_included_metric_S2_theta_phi(const double *const h_D0D0,const double *const h_D0D1,const double *const h_D1D1,double *const g_D0D0, double *const g_D0D1, double *const g_D0D2, double *const g_D1D1, double *const g_D1D2, double *const g_D2D2, double theta,const double phi);



