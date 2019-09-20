#include "bbn_headers.h"
#include "utilities_lib.h"

Grid_T *bbn_initialize_next_grid(Grid_T *const grid_prev);
static Grid_T *make_next_grid_using_previous_grid(Grid_T *const grid_prev);
static Grid_T *TOV_KerrShild_approximation(void);
static Grid_T *creat_grid_TOV_KerrShild(const double R_NS_l,const double R_BH_r,const double a_BH);
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,const double R_NS_l,const double R_BH_r,const double a_BH);
static void init_field_TOV_plus_KerrSchild(Grid_T *const grid,const TOV_T *const tov, const double a_BH, const double M_BH);
static void make_normal_vector_on_BH_horizon(Grid_T *const grid);
static void find_NS_surface(Grid_T *const grid);
void bbn_allocate_fields(Grid_T *const grid);
void bbn_partial_derivatives_fields(Grid_T *const grid);
void bbn_populate_free_data(Grid_T *const grid);
void bbn_update_psi10A_UiUj(Patch_T *const patch);
static void bbn_update_Aij(Grid_T *const grid);
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev);
static void find_Xp_and_patchp(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch);
static double interpolate_from_prev_grid(const char *const field,const double *const X,Patch_T *const patch);
static void find_Euler_eq_const_TOV_KerrSchild(Grid_T *const grid);


