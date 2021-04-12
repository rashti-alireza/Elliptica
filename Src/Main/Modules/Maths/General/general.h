#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_special_functions_lib.h"

/* these following are used for compatibility with CForm in mathematica */
#define E M_E
#define Cos(a) cos(a)
#define Sin(a) sin(a)
#define Cosh(a) cosh(a)
#define Sinh(a) sinh(a)
#define Log(a) log(a)
#define Power(a,b) pow(a,b)
#define Sqrt(a) sqrt(a)
#define Csc(a) (1./sin(a))
#define Cot(a) (1./tan(a))

double root_square(const Uint n, const double *const v2,const double *const v1);
long double root_square_long(const long Uint n, const double *const v2, const double *const v1);
double dot(const Uint n, const double *const v2,const double *const v1);
double sum_1_N_cos_ia(const Uint N, const double a);
double sum_0_N_dCi_dfj_by_Ti_q(const Uint N,const Uint j,const double q);
double sum_0_N_dCi_dfj_by_dTi_dq(const Uint N,const Uint j,const double q);
double d_dq_sum_1_N_cos_ixb_cos_ixa(const int N, const double b,const double a);
double MaxMag_d(const double a,const double b);
double L_inf(const Uint N,const double *const v);
double L1_norm(const Uint n, const double *const v2,const double *const v1);
double L2_norm(const Uint n, const double *const v2,const double *const v1);
double arctan(const double y,const double x);
void arctan_argument_signum(double *const y_sign,double *const x_sign,const double phi);

