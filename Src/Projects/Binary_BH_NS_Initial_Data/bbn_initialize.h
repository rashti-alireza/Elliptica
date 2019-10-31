#include "bbn_headers.h"
#include "utilities_lib.h"
#include "maths_approximation_lib.h"
#include "maths_analytic_lib.h"
#include "physics_observables_lib.h"
#include "maths_equation_solvings_lib.h"

#define Power(a,b) pow(a,b)
#define Sqrt(a) sqrt(a)

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

Grid_T *bbn_initialize_next_grid(Grid_T *const grid_prev);
static Grid_T *make_next_grid_using_previous_grid(Grid_T *const grid_prev);
static Grid_T *TOV_KerrSchild_approximation(void);
static Grid_T *creat_bbn_grid_CS(struct Grid_Params_S *const GridParams);
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void init_field_TOV_plus_KerrSchild(Grid_T *const grid,const TOV_T *const tov, const double a_BH, const double M_BH);
static void make_normal_vector_on_BH_horizon(Grid_T *const grid);
void bbn_allocate_fields(Grid_T *const grid);
void bbn_partial_derivatives_fields(Grid_T *const grid);
void bbn_populate_free_data(Grid_T *const grid);
void bbn_update_psi10A_UiUj(Patch_T *const patch);
static void bbn_update_Aij(Grid_T *const grid);
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev);
static void find_X_and_patch(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch);
static double interpolate_from_patch_prim(const char *const field,const double *const X,Patch_T *const patch);
static void find_Euler_eq_const_TOV_KerrSchild(Grid_T *const grid);
static void find_Euler_eq_const(Grid_T *const grid);
static double Euler_eq_const_rootfinder_eq(void *params,const double *const x);
static void extrapolate_fluid_fields_outsideNS_CS(Grid_T *const grid);
//static void find_NS_surface_CS_method_CS(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void find_NS_surface_Ylm_method_CS(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void find_XYZ_and_patch_of_theta_phi_NS_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid);
static void find_theta_phi_of_XYZ_NS_CS(double *const theta,double *const phi,const double *const X,const Flag_T side);
static void free_Grid_Params_S(struct Grid_Params_S *par);
static struct Grid_Params_S *init_GridParams(void);
static double CenterOfMass_for_P_ADM_root_finder_eq(void *params,const double *const x);
static void find_center_of_mass(Grid_T *const grid);
static void find_NS_center(Grid_T *const grid);
static double dh_dx0_root_finder_eq(void *params,const double *const x);
static double dh_dx1_root_finder_eq(void *params,const double *const x);
static double dh_dx2_root_finder_eq(void *params,const double *const x);
static void find_BH_NS_Omega_force_balance_eq(Grid_T *const grid);
static void adjust_NS_center(Grid_T *const grid);
static double bbn_NS_surface_enthalpy_eq(void *params,const double *const x);
static double bbn_NS_surface_denthalpy_dr(void *params,const double *const x,const unsigned dir);
static void extrapolate_fluid_fields_outsideNS(Grid_T *const grid);
static void find_NS_surface(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void update_B1_then_Beta_and_Aij(Grid_T *const grid,const double Omega_BHNS,const double Vr,const double y_CM,const double D);
