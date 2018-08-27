#include "core_lib.h"
//#include "coordinates_lib.h"
#include "macros_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
//#include "maths_analytic_lib.h"
//#include "maths_calculus_lib.h" 

typedef double Interpolation_Func_T(Interpolation_T *const interp_s);

Interpolation_T *init_interpolation(void);
double interpolation(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn(Interpolation_T *const interp_s);
static Interpolation_Func_T *ChooseInterpolationFunction(const Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_X(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_Y(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_Z(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_XY(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_XZ(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_YZ(Interpolation_T *const interp_s);
static double interpolation_Chebyshev_Tn_XYZ(Interpolation_T *const interp_s);











