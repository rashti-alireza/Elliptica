#include "core_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_analytic_lib.h"
#include "memory_managing_lib.h"
#include "maths_calculus_lib.h"
#include "maths_approximation_lib.h"
#include "maths_complex_numbers_lib.h"

void get_Ylm_coeffs(double *const realClm,double *const imagClm,const double *const f,const unsigned Ntheta,const unsigned Nphi,const unsigned Lmax);
static void get_Ylm_coeffs_GaussLegendre_EquiSpaced(double *const realClm,double *const imagClm,const double *const f,const unsigned Ntheta,const unsigned Nphi,const unsigned Lmax);
static double complex integrate_expImphi(const double *const f, const unsigned n,const int m);
void n2lm_Ylm(const int n, int *const l, int *const m,const int lmax);
int lm2n_Ylm(const int l,const int m, const int lmax);
unsigned lm2n(const unsigned l, const unsigned m);
double *alloc_ClmYlm(unsigned Lmax);
double interpolation_Ylm(const double *const realClm,const double *const imagClm,const unsigned Lmax, const double theta, const double phi);
double *df_dphi_Ylm(const double *const realClm,const double *const imagClm,const unsigned Ntheta, const unsigned Nphi,const unsigned Lmax);
double *df_dtheta_Ylm(const double *const realClm,const double *const imagClm,const unsigned Ntheta, const unsigned Nphi,const unsigned Lmax);

