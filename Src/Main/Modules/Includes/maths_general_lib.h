#ifndef maths_general_LIB_H
#define maths_general_LIB_H
#include "elliptica_system_lib.h"

/* x^2 */
#define Pow2(x) ((x)*(x))
/* x^3 */
#define Pow3(x) ((x)*(x)*(x))
/* x^4 */
#define Pow4(x) ((x)*(x)*(x)*(x))

/* absolute value for double data type */
#define ABSd(x) ((x) > 0. ? (x) : -(x))

#ifndef M_PI
#define M_PI (3.141592653589793)
#endif

#ifndef M_E
#define M_E  (2.718281828459045)
#endif

long double rms_l(const long Uint n, const double *v2, const double *v1);
double root_square(const Uint n, const double *const v2, const double *const v1);
double dot(const Uint n, const double *const v2, const double *const v1);
double sum_1_N_cos_ia(const Uint N, const double a);
double sum_0_N_dCi_dfj_by_Ti_q(const Uint N,const Uint j,const double q);
double sum_0_N_dCi_dfj_by_dTi_dq(const Uint N,const Uint j,const double q);
double d_dq_sum_1_N_cos_ixb_cos_ixa(const int N, const double b,const double a);
void summation_tests(void);
double MaxMag_d(const double a,const double b);
double L_inf(const Uint n,const double *const v);
double L1_norm(const Uint n, const double *const v2,const double *const v1);
double L2_norm(const Uint n, const double *const v2,const double *const v1);
double arctan(const double y,const double x);
void arctan_argument_signum(double *const y_sign,double *const x_sign,const double phi);


#endif


