#include "core_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_special_functions_lib.h"
#include "maths_calculus_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_complex_lib.h"

void get_Ylm_coeffs(double *const realClm,double *const imagClm,const double *const f,const Uint Ntheta,const Uint Nphi,const Uint Lmax);
static void get_Ylm_coeffs_GaussLegendre_EquiSpaced(double *const realClm,double *const imagClm,const double *const f,const Uint Ntheta,const Uint Nphi,const Uint Lmax);
static double complex integrate_expImphi(const double *const f, const Uint n,const int m);
void n2lm_Ylm(const int n, int *const l, int *const m,const int lmax);
int lm2n_Ylm(const int l,const int m, const int lmax);
Uint lm2n(const Uint l, const Uint m);
double *alloc_ClmYlm(Uint Lmax);
double interpolation_Ylm(const double *const realClm,const double *const imagClm,const Uint Lmax, const double theta, const double phi);
double *df_dphi_Ylm(const double *const realClm,const double *const imagClm,const Uint Ntheta, const Uint Nphi,const Uint Lmax);
double *df_dtheta_Ylm(const double *const realClm,const double *const imagClm,const Uint Ntheta, const Uint Nphi,const Uint Lmax);

