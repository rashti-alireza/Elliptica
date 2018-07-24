#include "core_lib.h"
//#include "coordinates_lib.h"
#include "macros_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
//#include "maths_general_lib.h"
//#include "maths_analytic_lib.h"
//#include "maths_calculus_lib.h" 

typedef double Interpolation_Func_T(const Field_T *const f);

double interpolation(const Field_T *const f);
static double spectral_interpolation(const Field_T *const f);
static Interpolation_Func_T *WhichSpectralInterpolation(const Field_T *const f);

