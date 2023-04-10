#include "core_lib.h"
#include "manifold_lib.h"
#include "macros_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_special_functions_lib.h"
#include "fields_lib.h"
#include "maths_linear_algebra_lib.h"
#include "maths_spectral_methods_lib.h"

#define Tx(i,x) T(n[0],i,x)
#define Ty(j,y) T(n[1],j,y)
#define Tz(k,z) T(n[2],k,z)

typedef fInterpolation_T *fPick_Func_T(Interpolation_T *const interp_s);

Interpolation_T *init_interpolation(void);
double execute_interpolation(Interpolation_T *const interp_struct);
void plan_interpolation(Interpolation_T *const interp_s);
static double T(const Uint n,const Uint i,const double x);
static fInterpolation_T *interpolation_Chebyshev_Tn(Interpolation_T *const interp_s);
static double interpolation_Neville_1d(Interpolation_T *const interp_s);
static double interpolation_natural_cubic_spline_1d(Interpolation_T *const interp_s);
static double interpolation_NCS_derivative(Interpolation_T *const interp_s);
static void find_coeffs_natural_cubic_spline_1d(Interpolation_T *const interp_s);
static void order_arrays_natural_cubic_spline_1d(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_X(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_Y(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_Z(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_XY(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_XZ(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_YZ(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_XYZ(Interpolation_T *const interp_s);
void free_interpolation(Interpolation_T *interp_s);

/////////////////////////////Extra interpolation methods.
void assign_interpolation_ptrs(Interpolation_T *const interp_s);
static void set_interp_order_flag(Interpolation_T *const interp_s, Uint flag);
static void set_interp_alloc_mem_flag(Interpolation_T *const interp_s, Uint flag);
void set_interp_warn_flag(Interpolation_T *const interp_s, Uint flag);
//static Uint get_order_flag(Interpolation_T *const interp_s);
//static Uint get_alloc_mem_flag(Interpolation_T *const interp_s);
//static Uint get_warn_flag(Interpolation_T *const interp_s);
static void order_arrays_spline_1d(Interpolation_T *const interp_s);
static double interpolation_finite_difference(Interpolation_T *const interp_s);
static Uint FDM_min(Uint n, Uint M);
static void find_coeffs_Hermite_cubic_spline(Interpolation_T *const interp_s);
static double interpolation_Hermite_cubic_spline(Interpolation_T *const interp_s);
static double interpolation_HCS_derivative(Interpolation_T *const interp_s);
static void find_coeffs_clamped_cubic_spline_1d(Interpolation_T *const interp_s);
static double interpolation_clamped_cubic_spline_1d(Interpolation_T *const interp_s);
static double interpolation_CCS_derivative(Interpolation_T *const interp_s);
//static double interpolation_log_linear(Interpolation_T *const interp_s);
//static double interpolation_log_derivative(Interpolation_T *const interp_s);
//static void prepare_log_interpolation(Interpolation_T *const interp_s);









