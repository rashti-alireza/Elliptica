#include <complex.h>
#include <tgmath.h>//for the type generate macros.

double complex Ylm(const unsigned l, int m, const double theta, const double phi);
double complex dYlm_dphi(const unsigned l, const int m, const double theta, const double phi);
double complex dYlm_dtheta(const unsigned l, const int m, const double theta, const double phi);


