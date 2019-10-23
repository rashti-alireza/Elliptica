#include "sbh_headers.h"
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
  double R_BH;/* BH radius */
  double a_BH;/* BH spin */
  double BH_r_center;/* center of right BH */
};

Grid_T *sbh_initialize_next_grid(Grid_T *const grid_prev);
static Grid_T *make_next_grid_using_previous_grid(Grid_T *const grid_prev);
static Grid_T *KerrShild_approximation(void);
static Grid_T *creat_sbh_grid_CS(struct Grid_Params_S *const GridParams);
static void BH_surface_CubedSpherical_grid(Grid_T *const grid,struct Grid_Params_S *const GridParams);
static void init_field_KerrSchild(Grid_T *const grid,const double a_BH, const double M_BH);
static void make_normal_vector_on_BH_horizon(Grid_T *const grid);
void sbh_allocate_fields(Grid_T *const grid);
void sbh_partial_derivatives_fields(Grid_T *const grid);
void sbh_populate_free_data(Grid_T *const grid);
void sbh_update_psi10A_UiUj(Patch_T *const patch);
static void sbh_update_Aij(Grid_T *const grid);
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev);
static void find_X_and_patch(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch);
static double interpolate_from_patch_prim(const char *const field,const double *const X,Patch_T *const patch);
static void free_Grid_Params_S(struct Grid_Params_S *par);
static struct Grid_Params_S *init_GridParams(void);
static void calculate_P_ADMs(Grid_T *const grid);
