#include "bbn_headers.h"
#include "utilities_lib.h"
#include "maths_approximation_lib.h"
#include "maths_analytic_lib.h"
#include "physics_observables_lib.h"
#include "maths_equation_solvings_lib.h"
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

#define MAX_ARR   400
#define MAX_ARRx2 2*MAX_ARR
#define MAX_ARRx3 3*MAX_ARR
#define MAX_ARRx4 4*MAX_ARR
#define MAX_ARRx5 5*MAX_ARR
#define Power(a,b) pow(a,b)
#define Sqrt(a) sqrt(a)
#define prep_and_call(x) REALLOC_v_WRITE_v(x)\
                         const double *const other_##x = patchp->pool[LookUpField_E(#x,patchp)]->v;
#define copy_values(x)   x[ijk] = other_##x[ijk];

/* handy macros for extrapolating inside BH */
#define STRING_IT(x)  #x

/* BH-filler */
#define WTGR_EXTRAPOLATE_scalar(x)   \
        /* u */ \
        double x##_onAH       = interpolate_from_patch_prim(STRING_IT(x)        ,X_on_BHsurf,BHsurf_patch); \
        /* du/dx^i */ \
        double d##x##D0_onAH = interpolate_from_patch_prim(STRING_IT(d##x##_D0),X_on_BHsurf,BHsurf_patch); \
        double d##x##D1_onAH = interpolate_from_patch_prim(STRING_IT(d##x##_D1),X_on_BHsurf,BHsurf_patch); \
        double d##x##D2_onAH = interpolate_from_patch_prim(STRING_IT(d##x##_D2),X_on_BHsurf,BHsurf_patch); \
        /* d^2u/d(x^i)^2 */ \
        double dd##x##D0D0_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##_D0D0),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D0D1_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##_D0D1),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D0D2_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##_D0D2),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D1D2_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##_D1D2),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D1D1_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##_D1D1),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D2D2_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##_D2D2),X_on_BHsurf,BHsurf_patch); \
        /* du/dr */ \
        double dur_##x = (N[0]*d##x##D0_onAH+N[1]*d##x##D1_onAH+N[2]*d##x##D2_onAH); \
        /* ddu/dr^2 */ \
        WTGR_ddu_r2(x); \
        double a_##x = (2*ddur2_##x - 6*dur_##x + 11*x##_onAH + 5*u0_##x)/16.; \
        double b_##x = (-2*ddur2_##x + 2*dur_##x + 15*(x##_onAH - u0_##x))/32.; \
        double c_##x = (-2*ddur2_##x + 6*dur_##x - 3*x##_onAH + 3*u0_##x)/16.;\
        double d_##x = (2*ddur2_##x - 2*dur_##x + x##_onAH - u0_##x)/32.; \
        x[ijk]       = a_##x*Cheb_Tn(0,t)+b_##x*Cheb_Tn(1,t)+c_##x*Cheb_Tn(2,t)+d_##x*Cheb_Tn(3,t);


/* BH-filler */
#define WTGR_EXTRAPOLATE_Beta(x)   \
        /* u */ \
        double x##_onAH      = interpolate_from_patch_prim(STRING_IT(x)       ,X_on_BHsurf,BHsurf_patch); \
        /* du/dx^i */ \
        double d##x##D0_onAH = interpolate_from_patch_prim(STRING_IT(d##x##D0),X_on_BHsurf,BHsurf_patch); \
        double d##x##D1_onAH = interpolate_from_patch_prim(STRING_IT(d##x##D1),X_on_BHsurf,BHsurf_patch); \
        double d##x##D2_onAH = interpolate_from_patch_prim(STRING_IT(d##x##D2),X_on_BHsurf,BHsurf_patch); \
        /* d^2u/d(x^i)^2 */ \
        double dd##x##D0D0_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##D0D0),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D0D1_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##D0D1),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D0D2_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##D0D2),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D1D2_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##D1D2),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D1D1_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##D1D1),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D2D2_onAH = interpolate_from_patch_prim(STRING_IT(dd##x##D2D2),X_on_BHsurf,BHsurf_patch); \
        /* du/dr */ \
        double dur_##x = (N[0]*d##x##D0_onAH+N[1]*d##x##D1_onAH+N[2]*d##x##D2_onAH); \
        /* ddu/dr^2 */ \
        WTGR_ddu_r2(x); \
        double a_##x = (2*ddur2_##x - 6*dur_##x + 11*x##_onAH + 5*u0_##x)/16.; \
        double b_##x = (-2*ddur2_##x + 2*dur_##x + 15*(x##_onAH - u0_##x))/32.; \
        double c_##x = (-2*ddur2_##x + 6*dur_##x - 3*x##_onAH + 3*u0_##x)/16.;\
        double d_##x = (2*ddur2_##x - 2*dur_##x + x##_onAH - u0_##x)/32.; \
        x[ijk]       = a_##x*Cheb_Tn(0,t)+b_##x*Cheb_Tn(1,t)+c_##x*Cheb_Tn(2,t)+d_##x*Cheb_Tn(3,t);


/* BH-filler */        
#define WTGR_EXTRAPOLATE_gammabar(x)   \
        /* u */ \
        double x##_onAH      = interpolate_from_patch_prim(STRING_IT(_##x)     ,X_on_BHsurf,BHsurf_patch); \
        /* du/dx^i */ \
        double d##x##D0_onAH = interpolate_from_patch_prim(STRING_IT(_d##x##D0),X_on_BHsurf,BHsurf_patch); \
        double d##x##D1_onAH = interpolate_from_patch_prim(STRING_IT(_d##x##D1),X_on_BHsurf,BHsurf_patch); \
        double d##x##D2_onAH = interpolate_from_patch_prim(STRING_IT(_d##x##D2),X_on_BHsurf,BHsurf_patch); \
        /* d^2u/d(x^i)^2 */ \
        double dd##x##D0D0_onAH = interpolate_from_patch_prim(STRING_IT(_dd##x##D0D0),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D0D1_onAH = interpolate_from_patch_prim(STRING_IT(_dd##x##D0D1),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D0D2_onAH = interpolate_from_patch_prim(STRING_IT(_dd##x##D0D2),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D1D2_onAH = interpolate_from_patch_prim(STRING_IT(_dd##x##D1D2),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D1D1_onAH = interpolate_from_patch_prim(STRING_IT(_dd##x##D1D1),X_on_BHsurf,BHsurf_patch); \
        double dd##x##D2D2_onAH = interpolate_from_patch_prim(STRING_IT(_dd##x##D2D2),X_on_BHsurf,BHsurf_patch); \
        /* du/dr */ \
        double dur_##x = (N[0]*d##x##D0_onAH+N[1]*d##x##D1_onAH+N[2]*d##x##D2_onAH); \
        /* ddu/dr^2 */ \
        WTGR_ddu_r2(x); \
        double a_##x = (2*ddur2_##x - 6*dur_##x + 11*x##_onAH + 5*u0_##x)/16.; \
        double b_##x = (-2*ddur2_##x + 2*dur_##x + 15*(x##_onAH - u0_##x))/32.; \
        double c_##x = (-2*ddur2_##x + 6*dur_##x - 3*x##_onAH + 3*u0_##x)/16.;\
        double d_##x = (2*ddur2_##x - 2*dur_##x + x##_onAH - u0_##x)/32.; \
        _##x[ijk]    = a_##x*Cheb_Tn(0,t)+b_##x*Cheb_Tn(1,t)+c_##x*Cheb_Tn(2,t)+d_##x*Cheb_Tn(3,t);

/* compute ddu/dr^2  */
#define WTGR_ddu_r2(a) \
        double ddur2_##a = 0; \
        {/* local */ \
          double KD[2]  = {0,1}; \
          double _ddur2_##a[3] = {0,0,0}; \
          double _x[3]  = {x,y,z}; \
          double _r3    = pow(r,3); \
          double _du[3] = {d##a##D0_onAH,d##a##D1_onAH,d##a##D2_onAH}; \
          double _ddu[3][3] = { \
            {dd##a##D0D0_onAH,dd##a##D0D1_onAH,dd##a##D0D2_onAH}, \
            {dd##a##D0D1_onAH,dd##a##D1D1_onAH,dd##a##D1D2_onAH}, \
            {dd##a##D0D2_onAH,dd##a##D1D2_onAH,dd##a##D2D2_onAH} \
                              }; \
          unsigned _i,_j; \
          for (_i = 0; _i < 3; ++_i) \
          { \
            for (_j = 0; _j < 3; ++_j) \
            { \
              _ddur2_##a[_i] += (KD[_i==_j]/r - _x[_i]*_x[_j]/_r3)*_du[_j]; \
              _ddur2_##a[_i] += N[_j]*_ddu[_j][_i]; \
            } \
            ddur2_##a += _ddur2_##a[_i]*N[_i]; \
          } \
        }
          

/* _gamma inverse */
#define COMPUTE_gammaI(a00,a01,a02,a10,a11,a12,a20,a21,a22) \
  { \
  _gammaI_U0U0[ijk] = (a11*a22 - a12*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  _gammaI_U0U1[ijk] = (-a01*a22 + a02*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  _gammaI_U0U2[ijk] = (a01*a12 - a02*a11)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  _gammaI_U1U1[ijk] = a00*(a00*a22 - a02*a20)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  _gammaI_U1U2[ijk] =-a00*(a00*a12 - a02*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  _gammaI_U2U2[ijk] = a00*(a00*a11 - a01*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  }

typedef void fAdjustment_t (Grid_T *const grid);

/* root finder struct for NS surface eq */
struct NS_surface_RootFinder_S
{
  Patch_T *patch;
  void *root_finder;
  double x0[3];/* (x,y,z) at the surface */
  double *N;/* the direction of increasing or decreasing of x = x0+N*d */
  //double Euler_C;/* Euler equation const. */
  //double scale;/* to avoid long step in root finder */
  //double maxR;/* max R allowed for NS surrounding */
};

/* this struct is adjustments and requirements for making of a new grid */
struct Grid_Params_S
{
  Grid_T *grid_prev;/* previous grid */
  double Max_R_NS_l;/* max of NS radius */
  double R_BH_r;/* BH radius */
  double a_BH;/* BH spin */
  double BH_r_center;/* center of right BH */
  double NS_l_center;/* center of left NS */
  const char *NS_R_type;/* type of NS which determines how to fill the radius field
                        // PerfectSphere
                        // SphericalHarmonic. 
                        // CubedSpherical. */
  const char *BH_R_type;/* type of BH which determines how to fill the radius field
                        // this gives you the shape of excision region (apparent horizon).
                        // PerfectSphere
                        // Boosted_Kerr-Schild. */
  /* if NS radius Ylm expanded */
  struct
  {
    double *realClm;
    double *imagClm;
    unsigned Lmax;
  }NS_R_Ylm[1];
  
};

/* root finder struct for finding center of mass */
struct CM_RootFinder_S
{
  Grid_T *grid;
  Observable_T *obs;
  double Omega_BHNS;
  double Vr;
  double D;
};

/* root finder structure for NS center */
struct NC_Center_RootFinder_S
{
  Patch_T *patch;
  Root_Finder_T *root_finder;
};

/* root finder struc for Euler eq const */  
struct Euler_eq_const_RootFinder_S
{
  Grid_T *grid;
  double NS_baryonic_mass;
};

/* root finder struct for {Px(y_CM) = 0 && Py(x_CM) = 0} */  
struct PxPy_RootFinder_S
{
  Grid_T *grid;
  Grid_T *freedata_grid;
  double x_CM0;
  double y_CM0;
};


Grid_T *bbn_initialize_next_grid(Grid_T *const grid_prev);
static Grid_T *make_next_grid_using_previous_grid(Grid_T *const grid_prev);
static Grid_T *TOV_KerrSchild_approximation(void);
static Grid_T *creat_bbn_grid_CS(struct Grid_Params_S *const GridParams);
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void init_field_TOV_plus_KerrSchild(Grid_T *const grid,const TOV_T *const tov);
void bbn_make_normal_vector_on_BH_horizon(Grid_T *const grid);
void bbn_add_fields(Grid_T *const grid);
void bbn_partial_derivatives_fields(Grid_T *const grid);
void bbn_populate_free_data(Grid_T *const grid);
void bbn_update_psi10A_UiUj(Patch_T *const patch);
void bbn_update_Aij(Grid_T *const grid);
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev);
static int find_X_and_patch(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch);
static double interpolate_from_patch_prim(const char *const field,const double *const X,Patch_T *const patch);
static void find_Euler_eq_const_TOV_KerrSchild(Grid_T *const grid);
static void find_Euler_eq_const(Grid_T *const grid);
static double Euler_eq_const_rootfinder_eq(void *params,const double *const x);
static void extrapolate_outsideNS_CS_exp_continuity_method(Grid_T *const grid);
static void extrapolate_outsideNS_CS_slop_method(Grid_T *const grid);
static void extrapolate_outsideNS_CS_Ylm_method(Grid_T *const grid,const char *const field_name);
static double interpolate_Clm_r_Ylm_3d(double *const realClm,double *const imagClm,const unsigned lmax,const double r,const double theta,const double phi);
static void extrapolate_fluid_fields_outsideNS(Grid_T *const grid);
static void find_NS_surface_Ylm_bisect_CS(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void find_NS_surface_Ylm_SD_CS(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void find_XYZ_and_patch_of_theta_phi_NS_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid);
static void find_theta_phi_of_XYZ_NS_CS(double *const theta,double *const phi,const double *const X,const Flag_T side);
static void free_Grid_Params_S(struct Grid_Params_S *par);
static struct Grid_Params_S *init_GridParams(void);
static void find_NS_center(Grid_T *const grid);
static void adjust_NS_center_draw_enthalpy(Grid_T *const grid);
static void adjust_NS_center_tune_enthalpy(Grid_T *const grid,const double dhx0,const double dhz0);
static void keep_NS_center_fixed(Grid_T *const grid);
static double dh_dx0_root_finder_eq(void *params,const double *const x);
static double dh_dx1_root_finder_eq(void *params,const double *const x);
static double dh_dx2_root_finder_eq(void *params,const double *const x);
static double bbn_NS_surface_enthalpy_eq(void *params,const double *const x);
static double bbn_NS_surface_denthalpy_dr(void *params,const double *const x,const unsigned dir);
static void extrapolate_fluid_fields_outsideNS(Grid_T *const grid);
static void find_NS_surface(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void update_B1_dB1_Beta_dBete_Aij_dAij(Grid_T *const grid);
void bbn_extrapolate_metric_fields_insideBH(Grid_T *const grid);
static void add_patches_insideBH(Grid_T *const grid);
static void extrapolate_insideBH_CS_linear(Grid_T *const grid);
static void parse_adjust_parameter(const char *const par,char *adjust[3]);
static void P_ADM_control(Grid_T *const grid);
fAdjustment_t *get_func_force_balance_adjustment(const char *const adjust);
fAdjustment_t *get_func_P_ADM_adjustment(const char *const adjust);
static void Pxy_ADM_is0_by_xy_CM_roots(Grid_T *const grid);
static void Px_ADM_is0_by_x_boost(Grid_T *const grid);
static void Py_ADM_is0_by_y_boost(Grid_T *const grid);
static void Pz_ADM_is0_by_z_boost(Grid_T *const grid);
static void Py_ADM_is0_by_x_CM(Grid_T *const grid);
static void Px_ADM_is0_by_y_CM(Grid_T *const grid);
static void Px_ADM_is0_by_BH_center_y(Grid_T *const grid);
static void Py_ADM_is0_by_BH_center_x(Grid_T *const grid);
static double x_CM_root_of_Py(void *params,const double *const x);
static double y_CM_root_of_Px(void *params,const double *const x);
static void force_balance_eq_root_finders(Grid_T *const grid,const int dir, const char *const par);
static void force_balance_eq(Grid_T *const grid);
static void force_balance_ddx_x_CM(Grid_T *const grid);
static void force_balance_ddy_x_CM(Grid_T *const grid);
static void force_balance_ddz_x_CM(Grid_T *const grid);
static void force_balance_ddx_y_CM(Grid_T *const grid);
static void force_balance_ddy_y_CM(Grid_T *const grid);
static void force_balance_ddz_y_CM(Grid_T *const grid);
static void force_balance_ddx_Omega(Grid_T *const grid);
static void force_balance_ddy_Omega(Grid_T *const grid);
static void force_balance_ddz_Omega(Grid_T *const grid);
static void adjust_AH_radius(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void adjust_BH_Omega(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static double AH_surface_function_PerfectSphere_CS(const double a,const double b,const double R,const double *const c,const Flag_T side);
static void move_geometry(Grid_T *const grid_next,Grid_T *const grid_prev);
static void move_solve_man_jacobian(Patch_T *const patch2,Patch_T *const patch1);
void bbn_free_grid_and_its_parameters(Grid_T *grid);
static Grid_T *load_checkpoint_file(void);
static int IsThereAnyUsefulCheckpointFile(void);
static void Pz_ADM_is0_by_BH_Vz(Grid_T *const grid);
static void extrapolate_insideBH_CS_C0_Ylm(Grid_T *const grid,const char *const field_name);
static void extrapolate_insideBH_CS_WTGR(Grid_T *const grid);
static void find_XYZ_and_patch_of_theta_phi_BH_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid);
static void force_balance_ddCM_Omega(Grid_T *const grid);
void bbn_create_grid_prototype_BC(Grid_T *const grid);
void bbn_bam_error(const char *const msg,const char *const file,const int line);




