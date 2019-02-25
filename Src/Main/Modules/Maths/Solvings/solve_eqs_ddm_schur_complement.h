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
#include "prints_lib.h"

#define DDM_SCHUR_COMPLEMENT_OpenMP(x) _Pragma ( #x )

typedef enum DDM_SCHUR_COMPLEMENT_FLAG_T
{
  ITS,
  OTHERS
}DDM_SC_Flag_T;

int ddm_schur_complement(Grid_T *const grid);
void test_solve_ddm_schur_complement(Grid_T *const grid);
fJs_T *get_j_reader(const Matrix_T *const m);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type);
fdInterp_dfs_T *get_dInterp_df(const Patch_T *const patch,const SubFace_T *const sf,const char *const dir);
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
static void making_F_by_f_prime(Patch_T *const patch);
static void making_F_by_E_prime(Patch_T *const patch);
static double *compute_g_prime(Grid_T *const grid);
static Matrix_T *compute_S(Grid_T *const grid);
static void compute_x(Patch_T *const patch);
static void solve_Sy_g_prime(Matrix_T *const S,double *const g_prime,Grid_T *const grid);
static void populate_F_and_C(Patch_T *const patch, Pair_T *const pair);
static void fill_C_F_collocation(Patch_T *const patch, Pair_T *const pair);
static void fill_C_F_interpolation(Patch_T *const patch, Pair_T *const pair);
static void miscellany_in_sewing(Patch_T *const patch);
static void set_NSs_NIs(Patch_T *const patch);
static void free_E_Trans_prime(Patch_T *const patch);
static void update_field(Patch_T *const patch);
static void free_x(Patch_T *const patch);
static void free_y(Grid_T *const grid);
static void making_B_single_patch(Patch_T *const patch);
static void solve_Bx_f(Patch_T *const patch);
static void update_field_single_patch(Patch_T *const patch);
static Flag_T check_residual_single_patch(const Patch_T *const patch,const double res_input);
static void checks_and_constraints(const Grid_T *const grid);
static Matrix_T *making_J_Schur_Method(Grid_T *const grid);
static Matrix_T *making_J_Old_Fashion(Grid_T *const grid);
static int solve_field_test(Grid_T *const grid);
static double *make_col_F(Grid_T *const grid);
static Matrix_T *making_J_Old_Fashion(Grid_T *const grid);
static int compare_Js(Grid_T *const grid,const Matrix_T *const J_Reg,const Matrix_T *const J_Schur);
static void free_schur_f_g(Grid_T *const grid);







