#include "bbn_headers.h"


#define IJ(i,j) ((j)+Nphi*(i))
/* math */
#define Sqrt(x)     sqrt((x))
#define Cos(x)      cos((x))
#define Sin(x)      sin((x))
#define Power(x,n)  pow((x),(n))

/* getting field name (xName) and put the interpolant into ixName */
#define INTERPOLATE_macro(xName) \
  { \
    Interpolation_T *interp_##xName = init_interpolation(); \
    interp_##xName->field = patch->pool[Ind(#xName)]; \
    interp_##xName->XY_dir_flag  = 1; \
    interp_##xName->X            = X[0]; \
    interp_##xName->Y            = X[1]; \
    interp_##xName->K            = type_flg == 1/* NS? */? patch->n[2]-1 : 0; \
    plan_interpolation(interp_##xName); \
    i##xName = execute_interpolation(interp_##xName); \
    free_interpolation(interp_##xName); \
  }

static void find_XYZ_and_patch_of_theta_phi_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid,const char *const type);
static void find_theta_phi_of_XYZ_CS(double *const theta,double *const phi,const double *const X,const Flag_T side);
void bbn_populate_2d_induced_metric_S2_theta_phi(double *const h_D0D0,double *const h_D0D1,double *const h_D1D1,const double *const g_D0D0,const double *const g_D0D1,const double *const g_D0D2,const double *const g_D1D1,const double *const g_D1D2,const double *const g_D2D2,const double r,const double theta,const double phi);
void bbn_compute_induced_metric_on_S2_CS_Ylm_CTS(Grid_T *const grid,const char *const type,const unsigned lmax,double **const ph_D0D0,double **const ph_D0D1,double **const ph_D1D1);
void bbn_test_induced_metric_algorithm(Grid_T *const grid);

void 
bbn_compute_induced_metric_on_S2_CS_FT_CTS
  (
  Grid_T *const grid,
  const char *const type,/* NS or BH */
  const unsigned Ntheta,/* # of collocation points in theta direction */
  const unsigned Nphi,/* # of collocation points in phi direction */
  double **const ph_D0D0,/* induced h00  pointer */
  double **const ph_D0D1,/* induced h01  pointer */
  double **const ph_D1D1 /* induced h11  pointer */
  );


void
bbn_compute_AKV_from_z
  (
  Grid_T *const grid/* grid */,
  const double *const akv/* akv scalar values */,
  const char *const dakv_D0/* d/dx akv name */,
  const char *const dakv_D1/* d/dy akv name */,
  const char *const dakv_D2/* d/dz akv name */,
  const char *const type/* NS or BH */,
  const unsigned Ntheta/* number of points in theta direction */,
  const unsigned Nphi/* number of points in theta direction */,
  const unsigned lmax/* l max in Ylm, if asked for spherical harmonic */,
  const int interpolation_type/* 1 double fourier, 0: spherical harmonic */
  );





