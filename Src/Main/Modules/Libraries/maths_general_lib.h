#define SQR(x) ((x)*(x))

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#ifndef M_E
#define M_E 2.718281828459045
#endif


long double rms_l(const unsigned long n, const double *v2, const double *v1);
double rms(const unsigned n, const double *const v2, const double *const v1);
double dot(const unsigned n, const double *const v2, const double *const v1);
double ABS(const double v);
double Cheb_Tn(const int n, const double x);
double Cheb_Un(const int n, const double x);
double d2T_dx2(const int n, const double x);