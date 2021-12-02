#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_calculus_lib.h"
#include "macros_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"


#define DO 1
#define DO_NOT 0


void test_dfs_df_spectral_vs_FiniteDiff(Grid_T *const grid);
void test_dfs_df_Spectral_vs_analytic(Grid_T *const grid);
void test_dInterp_a_df(Grid_T *const grid);
fdInterp_dfs_T *get_dInterp_df(const Patch_T *const patch,const SubFace_T *const sf,const char *const dir);
static void test_dInterp_x_df_YZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_y_df_YZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_z_df_YZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_df_YZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_x_df_XZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_y_df_XZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_z_df_XZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_df_XZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_x_df_XY_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_y_df_XY_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_z_df_XY_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_df_XY_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_x_df_XYZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_y_df_XYZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_z_df_XYZ_Tn_Ex(Field_T *const phi_field);
static void test_dInterp_df_XYZ_Tn_Ex(Field_T *const phi_field);
void test_root_finders(Grid_T *const grid);
static int root_finder_SteepestDescent(Grid_T *const grid);
static double root_finder_df0_dx_eq(void *params,const double *const x,const Uint dir);
static double root_finder_df1_dx_eq(void *params,const double *const x,const Uint dir);
static double root_finder_df2_dx_eq(void *params,const double *const x,const Uint dir);
static double root_finder_f0_eq(void *params,const double *const x);
static double root_finder_f1_eq(void *params,const double *const x);
static double root_finder_f2_eq(void *params,const double *const x);
void test_dfs_df_Spectral_vs_Spectral(Grid_T *const grid);

