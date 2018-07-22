#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"

double rms(const unsigned n, const double *const v2,const double *const v1);
long double rmsL(const unsigned long n, const double *const v2, const double *const v1);
double dot(const unsigned n, const double *const v2,const double *const v1);
double ABS(const double v);
double Cheb_U(const int n, const double x);
double Cheb_Tn(const int n, const double x);
double d2T_dx2(const int n, const double x);
