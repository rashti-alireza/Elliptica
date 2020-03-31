#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"

/* these following are used for compatibility with CForm in mathematica */
#define E M_E
#define Cos(a) cos(a)
#define Sin(a) sin(a)
#define Cosh(a) cosh(a)
#define Sinh(a) sinh(a)
#define Log(a) log(a)
#define Power(a,b) pow(a,b)
#define Sqrt(a) sqrt(a)
#define Csc(a) 1/sin(a)
#define Cot(a) 1/tan(a)

double root_square(const unsigned n, const double *const v2,const double *const v1);
long double rmsL(const long unsigned n, const double *const v2, const double *const v1);
double dot(const unsigned n, const double *const v2,const double *const v1);
double ABS(const double v);
double Cheb_U(const int n, const double x);
double Cheb_Tn(const int n, const double x);
double d2T_dx2(const int n, const double x);
double sum_1_N_cos_ia(const unsigned N, const double a);
double sum_0_N_dCi_dfj_by_Ti_q(const unsigned N,const unsigned j,const double q);
double sum_0_N_dCi_dfj_by_dTi_dq(const unsigned N,const unsigned j,const double q);
double d_dq_sum_1_N_cos_ixb_cos_ixa(const int N, const double b,const double a);
double MaxMag_d(const double a,const double b);
double L_inf(const unsigned N,const double *const v);
double L1_norm(const unsigned n, const double *const v2,const double *const v1);
double L2_norm(const unsigned n, const double *const v2,const double *const v1);
double arctan(const double y,const double x);
void arctan_argument_signum(double *const y_sign,double *const x_sign,const double phi);

