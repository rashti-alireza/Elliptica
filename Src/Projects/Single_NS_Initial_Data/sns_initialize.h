#include "sns_headers.h"
#include "utilities_lib.h"
#include "maths_approximation_lib.h"
#include "maths_analytic_lib.h"
#include "physics_observables_lib.h"
#include "maths_equation_solvings_lib.h"

#define Power(a,b) pow(a,b)
#define Sqrt(a) sqrt(a)

/* this struct is adjustments and requirements for making of a new grid */
struct Grid_Params_S
{
  double Max_R_NS_l;/* max of NS radius */
  double R_BH_r;/* BH radius */
  double a_BH;/* BH spin */
  double BH_r_center;/* center of right BH */
  double NS_l_center;/* center of left NS */
  const char *NS_R_type;/* type of NS which determines how to fill the radius field
                        // PerfectSphere
                        // SphericalHarmonic. 
                        // CubedSpherical*/
  /* if NS radius Ylm expanded */
  struct
  {
    double *realClm;
    double *imagClm;
    unsigned Lmax;
  }NS_R_Ylm[1];
  
};

/* root finder struc for Euler eq const */  
struct Euler_eq_const_RootFinder_S
{
  Grid_T *grid;
  double NS_baryonic_mass;
};

/* root finder structure for NS center */
struct NC_Center_RootFinder_S
{
  Patch_T *patch;
  Root_Finder_T *root_finder;
};

Grid_T *sns_initialize_next_grid(Grid_T *const grid_prev);
static Grid_T *make_next_grid_using_previous_grid(Grid_T *const grid_prev);
static Grid_T *TOV_approximation(void);
static Grid_T *creat_sns_grid_CS(struct Grid_Params_S *const GridParams);
static void NS_surface_CubedSpherical_grid(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void init_field_TOV(Grid_T *const grid,const TOV_T *const tov);
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev);
static void find_Xp_and_patchp(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch);
static double interpolate_from_patch_prim(const char *const field,const double *const X,Patch_T *const patch);
static void find_Euler_eq_const_TOV(Grid_T *const grid);
static void find_Euler_eq_const(Grid_T *const grid);
static double Euler_eq_const_rootfinder_eq(void *params,const double *const x);
static void extrapolate_fluid_fields_outsideNS_CS(Grid_T *const grid);
static void find_NS_surface_Ylm_method_CS(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void find_XYZ_and_patch_of_theta_phi_NS_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid);
static void find_theta_phi_of_XYZ_NS_CS(double *const theta,double *const phi,const double *const X,const Flag_T side);
static void free_Grid_Params_S(struct Grid_Params_S *par);
static struct Grid_Params_S *init_GridParams(void);
static void sns_update_Aij(Grid_T *const grid);
static void find_NS_center(Grid_T *const grid);
static double dh_dx0_root_finder_eq(void *params,const double *const x);
static double dh_dx1_root_finder_eq(void *params,const double *const x);
static double dh_dx2_root_finder_eq(void *params,const double *const x);
static void adjust_NS_center(Grid_T *const grid);
static double sns_NS_surface_enthalpy_eq(void *params,const double *const x);

