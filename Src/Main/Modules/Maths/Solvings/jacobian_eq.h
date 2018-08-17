#include "core_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "coordinates_lib.h"
#include "error_handling_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "memory_managing_lib.h"

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
typedef void Jacobian_eq_F(double **const J,Patch_T *const patch,JType_E jt_e);

static double SIGN[2] = {1.0,-1.0};

void make_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_make_jacobian_eq(Grid_T *const grid, const char * const* types);
static JType_E str2JType_E(const char *const str);
static void JType_E2str(const JType_E e,char *const str);
static void make_jacobian_spectral_method(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_spectral_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_spectral_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir);
static void make_jacobian_direct_method(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_direct_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_direct_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir);
static double ChebExtrema_1point(const unsigned n, const unsigned p);
static double dc0_df(const unsigned n0,const unsigned i,const unsigned l);
static double dc1_df(const unsigned n1,const unsigned j,const unsigned m);
static double dc2_df(const unsigned n2,const unsigned k,const unsigned n);
static double dT_dx(const int n,const double x);
static void read_1st_and_2nd_deriv(const JType_E deriv_dir,JType_E *const deriv_1st,JType_E *const deriv_2nd);
static void JType_E2Dd_T(const JType_E jt_e, Dd_T *const q_dir);

