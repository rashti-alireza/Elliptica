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

/* function for taking spectral derivative */
typedef double *SpecDerivative_Func_T(Field_T *const f,const Dd_T dir);

double *Partial_Derivative(Field_T *const f,const char *task);
static Method_T derivative_method(const char *const par,const char *const task);
static Method_T str2enum_method(const char *const str);
static Dd_T *derivative_direction(const char *const task,unsigned *const n);
static Dd_T str2enum_direction(const char *const str);
static double *take_spectral_derivative(Field_T *const f,const Dd_T  *const dir_e,const unsigned Ndir);
static double *spectral_derivative_1stOrder(Field_T *const f,const Dd_T dir_e);
static double *spectral_derivative_2ndOrder(Field_T *const f,const Dd_T dir_e);
static double *derivative_Chebyshev_Tn_1stOrder(Field_T *const f,const Dd_T dir);
static double *derivative_Chebyshev_Tn_2ndOrder(Field_T *const f,const Dd_T dir);
static double *make_1Dcollocation_ChebExtrema(const unsigned N);
static void get_dp_1stOrder(const Patch_T *const patch,SpecDerivative_Func_T **func,const Dd_T dir,Dd_T *dp);
static void get_dp_2ndOrder(const Patch_T *const patch,SpecDerivative_Func_T **func,const Dd_T dir,Dd_T *dp);
static void get_dependency(const Patch_T *const patch,const Dd_T dir, unsigned *dep);
static void get_SpecDerivative_func_1stOrder(const Patch_T *const patch,SpecDerivative_Func_T **func);
static void get_SpecDerivative_func_2ndOrder(const Patch_T *const patch,SpecDerivative_Func_T **func);
static unsigned IsSecondOrderFormula(Field_T *const f,const Dd_T *const dir_e,const unsigned Ndir);
static unsigned JacobianFormat_2ndOrder(const Patch_T *const patch,const Dd_T dir,Dd_T dp);
