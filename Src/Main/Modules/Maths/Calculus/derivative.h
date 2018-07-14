#include "core_lib.h"
#include "maths_approximation_lib.h"
#include "memory_managing_lib.h"
#include "macros_lib.h"
#include "error_handling_lib.h"
//#include "coordinates_lib.h"
//#include "utilities_lib.h"
//#include "maths_general_lib.h"
//#include "maths_analytic_lib.h"

typedef enum DERIVATIVE_T
{
  UNDEFINED_DERIVATIVE = 0,
  SPECTRAL,
  FINITE_DIFF
}Derivative_T;

double *Df(const Field_T *const f,const char *task);
static Derivative_T derivative_type(const char *const par,const char *const task);
static Derivative_T str2enum(const char *const str);