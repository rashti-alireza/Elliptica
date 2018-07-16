#include "core_lib.h"
#include "maths_approximation_lib.h"
#include "memory_managing_lib.h"
#include "macros_lib.h"
#include "error_handling_lib.h"
#include "coordinates_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
//#include "maths_analytic_lib.h"

/* methods of derivative */
typedef enum METHOD_T
{
  UNDEFINED_METHOD = 0,
  SPECTRAL,
  FINITE_DIFF
}Method_T;

/* directions of derivative */
typedef enum DIRECTION_T
{
  UNDEFINED_DIR = 0,
  x_DIR,
  y_DIR,
  z_DIR,
  a_DIR,
  b_DIR,
  c_DIR,
  N_DIR/* total number of direction */
}Direction_T;

double *Df(Field_T *const f,const char *task);
static Method_T derivative_method(const char *const par,const char *const task);
static Method_T str2enum_method(const char *const str);
static Direction_T *derivative_direction(const char *const task,unsigned *const n);
static Direction_T str2enum_direction(const char *const str);
static double *take_spectral_derivative(Field_T *const f,const Direction_T  *const dir_e,const unsigned Ndir);
static double *spectral_derivative_1d(Field_T *const f,const Direction_T dir_e);
static double *derivative_Chebyshev_Tn_1d(Field_T *const f,const Patch_T *const patch,const unsigned dir);