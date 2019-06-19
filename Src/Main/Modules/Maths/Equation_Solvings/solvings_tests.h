#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_approximation_lib.h"
#include "maths_calculus_lib.h"
#include "macros_lib.h"
#include "coordinates_lib.h"


#define DO 1
#define NOT_DO 0


void test_dfs_df_values(Grid_T *const grid);
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


