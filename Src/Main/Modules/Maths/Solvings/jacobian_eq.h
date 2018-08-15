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
  T_x,
  T_xx,
  T_y,
  T_yy,
  T_z,
  T_zz
}JType_E;

static double SIGN[2] = {1.0,-1.0};

static void fill_jacobian_spectral_method_x(double **const j,Patch_T *const patch);
static void fill_jacobian_spectral_method_xx(double **const j,Patch_T *const patch);
static void fill_jacobian_spectral_method_y(double **const j,Patch_T *const patch);
static void fill_jacobian_spectral_method_yy(double **const j,Patch_T *const patch);
static void fill_jacobian_spectral_method_z(double **const j,Patch_T *const patch);
static void fill_jacobian_spectral_method_zz(double **const j,Patch_T *const patch);
static double ChebExtrema_1point(const unsigned n, const unsigned p);
static double dc0_df(const unsigned n0,const unsigned i,const unsigned l);
static double dc1_df(const unsigned n1,const unsigned j,const unsigned m);
static double dc2_df(const unsigned n2,const unsigned k,const unsigned n);
static double dT_dx(const int n,const double x);


