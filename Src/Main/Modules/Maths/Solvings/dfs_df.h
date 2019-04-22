#include "core_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "coordinates_lib.h"
#include "error_handling_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "memory_managing_lib.h"
#include "maths_linear_algebra_lib.h"

#define _MAX_STR_ 400

/* Jacobian type */
typedef enum JTYPE_E
{
  T_x/* df_x/df */,
  T_xx/* df_xx/df */,
  T_y/* df_y/df */,
  T_yy/* df_yy/df */,
  T_z/* df_z/df */,
  T_zz/* df_zz/df */,
  T_UNDEF
}JType_E;

/* Jacobain for equation */
typedef void Js_Jacobian_eq_F(double **const J,Patch_T *const patch,JType_E jt_e);

void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
void make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void obsolete_fill_jacobian_spectral_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type);
double read_matrix_entry_ccs(Matrix_T *const m, const long r,const long c);
fdInterp_dfs_T *get_dInterp_df(const Patch_T *const patch,const SubFace_T *const sf,const char *const dir);
static JType_E str2JType_E(const char *const str);
static void JType_E2str(const JType_E e,char *const str);
static void make_jacobian_spectral_method(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_spectral_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_spectral_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir);
static void make_jacobian_direct_method(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_direct_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_direct_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir);
static double ChebExtrema_1point(const unsigned n, const unsigned p);
static double dc_df(const unsigned n,const unsigned i,const unsigned l);
static double dT_dx(const int n,const double x);
static void read_1st_and_2nd_deriv(const JType_E deriv_dir,JType_E *const deriv_1st,JType_E *const deriv_2nd);
static void JType_E2Dd_T(const JType_E jt_e, Dd_T *const q_dir);
static void write_J_in_disk_ccs(void);
static double J_sizeMb_ccs(const Matrix_T *const m);
fJs_T *get_j_reader(const Matrix_T *const m);
static double dInterp_x_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_y_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_z_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_x_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_y_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_z_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_x_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_y_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_z_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_x_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_y_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_z_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
static double dInterp_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const unsigned df);
