#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "maths_calculus_lib.h"
#include "manifold_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_matrix_solvers_lib.h"
#include "maths_general_lib.h"
#include "maths_linear_algebra_lib.h"
#include "maths_equation_solvings_lib.h"
#include "fields_lib.h"

#define DDM_SCHUR_COMPLEMENT_OpenMP(x) _Pragma ( #x )

typedef enum DDM_SCHUR_COMPLEMENT_FLAG_T
{
  ITS,
  OTHERS
}DDM_SC_Flag_T;

int ddm_schur_complement(Solve_Equations_T *const SolveEqs);
void test_Jacobian_of_equations(Solve_Equations_T *const SolveEqs);
fJs_T *get_j_reader(const Matrix_T *const m);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type);
fdInterp_dfs_T *get_dInterp_df(const Patch_T *const patch,const SubFace_T *const sf,const char *const dir);
static void preparing_ingredients(Solve_Equations_T *const SolveEqs);
static void make_map_and_inv(Patch_T *const patch);
static void calculate_residual(Grid_T *const grid);
static void calculate_residual_single_patch(Patch_T *const patch);
static void set_solving_man_settings_Frms_i(Grid_T *const grid);
static void set_solving_man_settings_Frms_i_single_patch(Patch_T *const patch);
static void set_solving_man_settings_solver_step(Grid_T *const grid,const int current_step);
static void free_solving_man_settings(Grid_T *const grid);
static void make_f(Patch_T *const patch);
char **get_solving_field_name(const char *const solving_order,unsigned *const nf);
static int solve_field(Solve_Equations_T *const SolveEqs);
static void set_solving_man_settings(Solve_Equations_T *const SolveEqs);
static void set_solving_man_cf(Solve_Equations_T *const SolveEqs);
static void f_in_equation_part(Patch_T *const patch);
static void f_in_boundary_part(Patch_T *const patch);
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
static char *making_B_and_E(Patch_T *const patch);
static char *making_E_prime_and_f_prime(Patch_T *const patch);
static char *making_F_and_C(Patch_T *const patch);
static char *making_F_by_f_prime(Patch_T *const patch);
static char *making_F_by_E_prime(Patch_T *const patch);
static double *compute_g_prime(Grid_T *const grid);
static Matrix_T *compute_S(Grid_T *const grid);
static Matrix_T *compute_S_CCS_long(Grid_T *const grid);
static Matrix_T *compute_S_CCS(Grid_T *const grid);
static void compute_x(Patch_T *const patch);
static char *solve_Sy_g_prime(Matrix_T *const S,double *const g_prime,Grid_T *const grid);
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
static void checks_and_constraints(const Grid_T *const grid);
static Matrix_T *making_J_Schur_Method(Solve_Equations_T *const SolveEqs);
static Matrix_T *making_J_Old_Fashion(Solve_Equations_T *const SolveEqs);
static int Jwritten_vs_Jequation(Solve_Equations_T *const SolveEqs);
static double *make_col_F(Grid_T *const grid);
static int compare_Js(Grid_T *const grid,const Matrix_T *const J_Reg,const Matrix_T *const J_Schur);
static void free_schur_f_g(Grid_T *const grid);
static void making_F_and_C_Regular(Patch_T *const patch);
static void pr_intro_ddm_schur_complement(void);
void calculate_equation_residual(Solve_Equations_T *const SolveEqs);
void sync_patch_pools(const Grid_T*const latest_grid,Solve_Equations_T *const solve);
void free_patch_SolMan_method_Schur(Patch_T *const patch);
Sewing_T *alloc_sewing(void);






