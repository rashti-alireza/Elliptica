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











