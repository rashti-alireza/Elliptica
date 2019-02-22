#define SQR(x) ((x)*(x))

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#ifndef M_E
#define M_E 2.718281828459045
#endif


long double rms_l(const long unsigned n, const double *v2, const double *v1);
double rms(const unsigned n, const double *const v2, const double *const v1);
double dot(const unsigned n, const double *const v2, const double *const v1);
double ABS(const double v);
double Cheb_Tn(const int n, const double x);
double Cheb_Un(const int n, const double x);
double d2T_dx2(const int n, const double x);
double sum_1_N_cos_ia(const unsigned N, const double a);
double sum_0_N_dCi_dfj_by_Ti_q(const unsigned N,const unsigned j,const double q);
double sum_0_N_dCi_dfj_by_dTi_dq(const unsigned N,const unsigned j,const double q);
double d_dq_sum_1_N_cos_ixb_cos_ixa(const int N, const double b,const double a);
void summation_tests(void);
double MaxMag_d(const double a,const double b);
double L_inf(const double *const stack,const double N);
double L1_norm(const unsigned n, const double *const v2,const double *const v1);
double L2_norm(const unsigned n, const double *const v2,const double *const v1);

