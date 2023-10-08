#include "core_lib.h"
#include "manifold_lib.h"
#include "macros_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_special_functions_lib.h"
#include "maths_calculus_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_complex_lib.h"
#include "fields_lib.h"

#define DO (1)
#define DO_NOT (0)
#define IJ(i,j,n)  ((j)+(i)*(n))
#define Cos(a) cos(a)
#define Sin(a) sin(a)
#define Power(a,b) pow(a,b)

/* types of derivatives; new one "must" be added to one before the last */
enum FUNC_E
{
  FUNC = 0,
  FUNC_x = 1,/* this number must to be changed */
  FUNC_y,
  FUNC_z,
  FUNC_xx,
  FUNC_yy,
  FUNC_zz,
  FUNC_xy,
  FUNC_xz,
  FUNC_yz,
  FUNC_xyz,
  N_FUNC
};

struct Error_S
{
  double E_an;/* analytic */
  double E_nu;/* numeric */
};

int derivative_tests(Grid_T *const grid);
int interpolation_tests(Grid_T *const grid);
Interpolation_T *init_interpolation(void);
double execute_interpolation(Interpolation_T *const interp_s);
void plan_interpolation(Interpolation_T *const interp_s);
static Flag_T read_F(sFunc_Patch2Pdouble_T **const F,sFunc_Patch2Pdouble_T **const DataBase_func,const enum FUNC_E fn);
static void enum2strcat(enum FUNC_E e,char *const f_derivative);
static void enum2str(enum FUNC_E e,char *const str);
static Flag_T compare_derivative(const char *const name,const double *const numc,const double *const anac,const Field_T *const func,const enum FUNC_E fn,const Patch_T *const patch,const char *const path,struct Error_S *const error);
static double *make_random_number(const Uint N,const double min,const double max);
static int interpolation_tests_X(Field_T *const field,const double *const X,const Uint N);
static int interpolation_tests_Y(Field_T *const field,const double *const Y,const Uint N);
static int interpolation_tests_Z(Field_T *const field,const double *const Z,const Uint N);
static int interpolation_tests_XY(Field_T *const field,const double *const X,const double *const Y,const Uint Nx,const Uint Ny);
static int interpolation_tests_XZ(Field_T *const field,const double *const X,const double *const Z,const Uint Nx,const Uint Nz);
static int interpolation_tests_YZ(Field_T *const field,const double *const Y,const double *const Z,const Uint Ny,const Uint Nz);
static int interpolation_tests_XYZ(Field_T *const field,const double *const X,const double *const Y,const double *const Z,const Uint Nx,const Uint Ny,const Uint Nz);
static int interpolation_tests_Neville_1d(void);
static int interpolation_tests_N_cubic_spline_1d(void);
static Uint order_of_derivative(const enum FUNC_E fn);
static double calculate_expected_precision_for_derivative(const Field_T *const func,const enum FUNC_E fn,const Patch_T *const patch);
int fourier_transformation_tests(Grid_T *const grid);
static int cft_c2r_r2c_1d_EquiSpaced_test(Grid_T *const grid);
static int Ylm_trans_test(Grid_T *const grid);
int Ylm_transformation_tests(Grid_T *const grid);
static int Ylm_derivatives_test(Grid_T *const grid);
static void free_func_Patch2Pdouble(sFunc_Patch2Pdouble_T **func);
static int r2cft_2d_EquiSpaced_test(Grid_T *const grid);
static int r2cft_2d_EquiSpaced_S2_test(Grid_T *const grid);
static int interpolation_tests_Hermite_1d(void);
static double f_poly_3deg1(const double x) __attribute__((unused));
static double df_poly_3deg1(const double x) __attribute__((unused));
static double ddf_poly_3deg1(const double x) __attribute__((unused));
static void derivative_tests_spectral_method(Grid_T *const grid);
static void derivative_tests_finite_diff_method(Grid_T *const grid);




