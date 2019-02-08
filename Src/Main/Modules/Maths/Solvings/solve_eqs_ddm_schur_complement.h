#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "maths_calculus_lib.h"
#include "coordinates_lib.h"
#include "maths_approximation_lib.h"
#include "memory_managing_lib.h"
#include "maths_solvers_lib.h"
#include "maths_general_lib.h"
#include "maths_linear_algebra_lib.h"

#define DDM_SCHUR_COMPLEMENT_OpenMP(x) _Pragma ( #x )

typedef enum DDM_SCHUR_COMPLEMENT_FLAG_T
{
  ITS,
  OTHERS
}DDM_SC_Flag_T;

int ddm_schur_complement(Grid_T *const grid);
fJs_T *get_j_reader(const Matrix_T *const m);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type);
double dinterp_x_df(const Patch_T *const patch,const double *const X,const unsigned df);
double dinterp_y_df(const Patch_T *const patch,const double *const X,const unsigned df);
double dinterp_z_df(const Patch_T *const patch,const double *const X,const unsigned df);
double dinterp_df_1d_index(const Patch_T *const patch,const double *const X,const unsigned df);
double dinterp_df_2d_index(const Patch_T *const patch,const double *const X,const unsigned df);
double dinterp_df_3d_index(const Patch_T *const patch,const double *const X,const unsigned df);
static void preparing_ingredients(Grid_T *const grid);
static void make_map_and_inv(Patch_T *const patch);
static Flag_T check_residual(const Grid_T *const grid,const double res_input);
static void make_f(Patch_T *const patch);
static char **read_fields_in_order(unsigned *const nf);
static int solve_field(Grid_T *const grid);
static void set_cf(Grid_T *const grid,const char *const field_name);
static void f_in_equation_part(Patch_T *const patch);
static void f_in_outerboundary_part(Patch_T *const patch);
static void make_partial_g(Patch_T *const patch);
static void pg_collocation(Patch_T *const patch, Pair_T *const pair);
static void pg_interpolation(Patch_T *const patch, Pair_T *const pair);
static void make_pg(Patch_T *const patch, Pair_T *const pair);
static void make_g(Grid_T *const grid);
static void populate_sewing(Patch_T *const patch);
static void make_others_sewing(const Patch_T *const patch,const Patch_T *const patch2,Sewing_T **const sewing);
static void make_its_sewing(const Patch_T *const patch,Sewing_T **const sewing);
static void populate_pair(Sewing_T *const sewing,SubFace_T *const subface,const DDM_SC_Flag_T flag);
static Pair_T *find_pair_in_sewing(const Sewing_T *const sewing,const SubFace_T *const subface);
static void mirror_pairs(Patch_T *const patch);
static unsigned const_index_of_face(Patch_T *const patch,const SubFace_T *const sf);
static void fill_interpolation_flags(Interpolation_T *const it,Patch_T *const patch,const SubFace_T *const sf);
static unsigned OnFace(const unsigned *const n, const unsigned p);
static void making_B_and_E(Patch_T *const patch);
static void making_E_prime_and_f_prime(Patch_T *const patch);
static void making_F_and_C(Patch_T *const patch);
static void populate_F_and_C(Patch_T *const patch, Pair_T *const pair);
static void fill_C_F_collocation(Patch_T *const patch, Pair_T *const pair);
static void fill_C_F_interpolation(Patch_T *const patch, Pair_T *const pair);
static void others_in_sewing(Patch_T *const patch);
static double J_dInterp_df(const Patch_T *const patch,const SubFace_T *const subface,const double *const X,const unsigned df,const char *const dir);





