#include "bbn_headers.h"
#include "utilities_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_approximation_lib.h"
#include "maths_analytic_lib.h"

struct Grid_Params_S
{
  double Max_R_NS_l;/* max of NS radius */
  double R_BH_r;/* BH radius */
  double a_BH;/* BH spin */
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

Grid_T *bbn_initialize_next_grid(Grid_T *const grid_prev);
static Grid_T *make_next_grid_using_previous_grid(Grid_T *const grid_prev);
static Grid_T *TOV_KerrShild_approximation(void);
static Grid_T *creat_grid_CS(struct Grid_Params_S *const GridParams);
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void init_field_TOV_plus_KerrSchild(Grid_T *const grid,const TOV_T *const tov, const double a_BH, const double M_BH);
static void make_normal_vector_on_BH_horizon(Grid_T *const grid);
void bbn_allocate_fields(Grid_T *const grid);
void bbn_partial_derivatives_fields(Grid_T *const grid);
void bbn_populate_free_data(Grid_T *const grid);
void bbn_update_psi10A_UiUj(Patch_T *const patch);
static void bbn_update_Aij(Grid_T *const grid);
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev);
static void find_Xp_and_patchp(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch);
static double interpolate_from_prev_grid(const char *const field,const double *const X,Patch_T *const patch);
static void find_Euler_eq_const_TOV_KerrSchild(Grid_T *const grid);
static void find_Euler_eq_const(Grid_T *const grid);
static double Euler_eq_const_rootfinder_eq(void *params,const double *const x);
static void extrapolate_fluid_fields_outsideNS_CS(Grid_T *const grid);
static void find_NS_surface_CS_method_CS(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void find_NS_surface_Ylm_method_CS(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void find_XYZ_of_theta_phi_NS_CS(double *const X,const double theta,const double phi,Patch_T *const patch);
static double XYZ_of_theta_phi_NS_CS_RT_EQ(void *params,const double *const dr);
static Patch_T *find_patch_of_theta_phi_NS_CS(const double theta,const double phi,Grid_T *const grid);
