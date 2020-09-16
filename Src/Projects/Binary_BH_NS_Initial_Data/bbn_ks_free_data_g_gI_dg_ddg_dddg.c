#include "bbn_headers.h"
#include "bbn_ks_free_date_analytic.h"
void bbn_free_data_g_gI_analytic(
        Patch_T *const patch,
        double *(*get_v)(const char *const fname,void *params),
        void *params);
void bbn_free_data_dg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params);
void bbn_free_data_ddg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params);
void bbn_free_data_dddg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params);
void bbn_free_data_g_gI_analytic(
 Patch_T *const patch,
 double *(*get_v)(const char *const fname,void *params),
 void *params)
{
  const double r0          = Pgetd("BH_KerrSchild_RollOff");
  const double BH_center_x = Pgetd("BH_center_x");
  const double BH_center_y = Pgetd("BH_center_y");
  const double BH_center_z = Pgetd("BH_center_z");
  const double M_BH        = Pgetd("BH_irreducible_mass");
  const double a_BH        = Pgetd("BH_net_spin");
  const double chi_U0   = Pgetd("BH_chi_U0");
  const double chi_U1   = Pgetd("BH_chi_U1");
  const double chi_U2   = Pgetd("BH_chi_U2");
  const double y_CM = Pgetd("y_CM");
  const double x_CM = Pgetd("x_CM");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double chi = sqrt(Pow2(chi_U0)+Pow2(chi_U1)+Pow2(chi_U2));
  const unsigned nn = patch->nn;
  double phiy = 0,phiz = 0;
  double Bx,By,Bz,B2;/* B = v/c */
  unsigned ijk;

  assert(LSSEQL(chi,1));

  /* boost */
  Bx = -Omega_BHNS*(BH_center_y-y_CM);
  By =  Omega_BHNS*(BH_center_x-x_CM);
  Bz = Pgetd("BH_Vz");
  B2 = Pow2(Bx)+Pow2(By)+Pow2(Bz);

  /* rotation */
  if (!EQL(chi,0))/* otherwise tR is 0 */
  {
    phiz = -arctan(chi_U1,chi_U0);
    phiy = -acos(chi_U2/chi);
    assert(isfinite(phiy));
  }
    double *const _gamma_D2D2 = get_v("_gamma_D2D2",params);
    double *const _gamma_D0D2 = get_v("_gamma_D0D2",params);
    double *const _gamma_D0D0 = get_v("_gamma_D0D0",params);
    double *const _gamma_D0D1 = get_v("_gamma_D0D1",params);
    double *const _gamma_D1D2 = get_v("_gamma_D1D2",params);
    double *const _gamma_D1D1 = get_v("_gamma_D1D1",params);
    double *const _gammaI_U0U2 = get_v("_gammaI_U0U2",params);
    double *const _gammaI_U0U0 = get_v("_gammaI_U0U0",params);
    double *const _gammaI_U0U1 = get_v("_gammaI_U0U1",params);
    double *const _gammaI_U1U2 = get_v("_gammaI_U1U2",params);
    double *const _gammaI_U1U1 = get_v("_gammaI_U1U1",params);
    double *const _gammaI_U2U2 = get_v("_gammaI_U2U2",params);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x,y,z,r,H,k0,k1,k2,kt;
      x = patch->node[ijk]->x[0]-BH_center_x;
      y = patch->node[ijk]->x[1]-BH_center_y;
      z = patch->node[ijk]->x[2]-BH_center_z;
      r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
      _gamma_D0D0[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
bbn_ks_c(x, y, z)*pow(bbn_ks_k0(x, y, z), 2) + 1.0;
      _gamma_D0D1[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z);
      _gamma_D0D2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z);
      _gamma_D1D1[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k1
bbn_ks_c(x, y, z)*pow(bbn_ks_k1(x, y, z), 2) + 1.0;
      _gamma_D1D2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z);
      _gamma_D2D2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k2
bbn_ks_c(x, y, z)*pow(bbn_ks_k2(x, y, z), 2) + 1.0;
      _gammaI_U0U0[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
((pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*bbn_ks_c(x, y, z) + 1)/
((pow(bbn_ks_k0(x, y, z), 2) + pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*
bbn_ks_c(x, y, z) + 1);
      _gammaI_U0U1[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
-bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)/((pow(bbn_ks_k0(x, y, z), 2) + 
pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*bbn_ks_c(x, y, z) + 
1);
      _gammaI_U0U2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
-bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)/((pow(bbn_ks_k0(x, y, z), 2) + 
pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*bbn_ks_c(x, y, z) + 
1);
      _gammaI_U1U1[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
((pow(bbn_ks_k0(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*bbn_ks_c(x, y, z) + 1)/
((pow(bbn_ks_k0(x, y, z), 2) + pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*
bbn_ks_c(x, y, z) + 1);
      _gammaI_U1U2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
-bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)/((pow(bbn_ks_k0(x, y, z), 2) + 
pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*bbn_ks_c(x, y, z) + 
1);
      _gammaI_U2U2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
((pow(bbn_ks_k0(x, y, z), 2) + pow(bbn_ks_k1(x, y, z), 2))*bbn_ks_c(x, y, z) + 1)/
((pow(bbn_ks_k0(x, y, z), 2) + pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*
bbn_ks_c(x, y, z) + 1);
      /* quick test check _gamma * _gammaI = delta */
      if (0)
      {
          double delta_U0D0 = 
        _gammaI_U0U0[ijk]*_gamma_D0D0[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D0D2[ijk];
          double delta_U0D1 = 
        _gammaI_U0U0[ijk]*_gamma_D0D1[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D1D2[ijk];
          double delta_U0D2 = 
        _gammaI_U0U0[ijk]*_gamma_D0D2[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U0U2[ijk]*_gamma_D2D2[ijk];
          double delta_U1D2 = 
        _gammaI_U0U1[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U1U2[ijk]*_gamma_D2D2[ijk];
          double delta_U1D0 = 
        _gammaI_U0U1[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D0D2[ijk];
         double delta_U1D1 = 
        _gammaI_U0U1[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D1D2[ijk];
        
          double delta_U2D2 = 
        _gammaI_U0U2[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U2U2[ijk]*_gamma_D2D2[ijk];
          double delta_U2D0 = 
        _gammaI_U0U2[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D0D2[ijk];
          double delta_U2D1 = 
        _gammaI_U0U2[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D1D2[ijk];
        if(!EQL(delta_U1D1,1)||!isfinite(delta_U1D1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D1,0)||!isfinite(delta_U0D1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D2,0)||!isfinite(delta_U0D2))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U1D2,0)||!isfinite(delta_U1D2))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D0,1)||!isfinite(delta_U0D0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D1,0)||!isfinite(delta_U2D1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D2,1)||!isfinite(delta_U2D2))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D0,0)||!isfinite(delta_U2D0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U1D0,0)||!isfinite(delta_U1D0))  Error0("_gammaI is not correct!\n");
      }
    }
}
void bbn_free_data_dg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params)
{
  const double r0          = Pgetd("BH_KerrSchild_RollOff");
  const double BH_center_x = Pgetd("BH_center_x");
  const double BH_center_y = Pgetd("BH_center_y");
  const double BH_center_z = Pgetd("BH_center_z");
  const double M_BH        = Pgetd("BH_irreducible_mass");
  const double a_BH        = Pgetd("BH_net_spin");
  const double chi_U0   = Pgetd("BH_chi_U0");
  const double chi_U1   = Pgetd("BH_chi_U1");
  const double chi_U2   = Pgetd("BH_chi_U2");
  const double y_CM = Pgetd("y_CM");
  const double x_CM = Pgetd("x_CM");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double chi = sqrt(Pow2(chi_U0)+Pow2(chi_U1)+Pow2(chi_U2));
  const unsigned nn = patch->nn;
  double phiy = 0,phiz = 0;
  double Bx,By,Bz,B2;/* B = v/c */
  unsigned ijk;

  assert(LSSEQL(chi,1));

  /* boost */
  Bx = -Omega_BHNS*(BH_center_y-y_CM);
  By =  Omega_BHNS*(BH_center_x-x_CM);
  Bz = Pgetd("BH_Vz");
  B2 = Pow2(Bx)+Pow2(By)+Pow2(Bz);

  /* rotation */
  if (!EQL(chi,0))/* otherwise tR is 0 */
  {
    phiz = -arctan(chi_U1,chi_U0);
    phiy = -acos(chi_U2/chi);
    assert(isfinite(phiy));
  }
  
double *const _dgamma_D1D2D2 = get_v("_dgamma_D1D2D2",params);
double *const _dgamma_D0D0D1 = get_v("_dgamma_D0D0D1",params);
double *const _dgamma_D0D2D1 = get_v("_dgamma_D0D2D1",params);
double *const _dgamma_D0D1D0 = get_v("_dgamma_D0D1D0",params);
double *const _dgamma_D1D2D1 = get_v("_dgamma_D1D2D1",params);
double *const _dgamma_D2D2D0 = get_v("_dgamma_D2D2D0",params);
double *const _dgamma_D0D0D0 = get_v("_dgamma_D0D0D0",params);
double *const _dgamma_D0D0D2 = get_v("_dgamma_D0D0D2",params);
double *const _dgamma_D0D2D2 = get_v("_dgamma_D0D2D2",params);
double *const _dgamma_D2D2D1 = get_v("_dgamma_D2D2D1",params);
double *const _dgamma_D0D1D1 = get_v("_dgamma_D0D1D1",params);
double *const _dgamma_D0D2D0 = get_v("_dgamma_D0D2D0",params);
double *const _dgamma_D1D2D0 = get_v("_dgamma_D1D2D0",params);
double *const _dgamma_D1D1D1 = get_v("_dgamma_D1D1D1",params);
double *const _dgamma_D0D1D2 = get_v("_dgamma_D0D1D2",params);
double *const _dgamma_D1D1D0 = get_v("_dgamma_D1D1D0",params);
double *const _dgamma_D1D1D2 = get_v("_dgamma_D1D1D2",params);
double *const _dgamma_D2D2D2 = get_v("_dgamma_D2D2D2",params);

  for (ijk = 0; ijk < nn; ++ijk) 
  {
    double x,y,z;
    x = patch->node[ijk]->x[0]-BH_center_x;
    y = patch->node[ijk]->x[1]-BH_center_y;
    z = patch->node[ijk]->x[2]-BH_center_z;

    _dgamma_D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z);

    _dgamma_D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k0
(2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y))*bbn_ks_k0(x, y, z);

    _dgamma_D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y) + bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y);

    _dgamma_D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x) + bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x);

    _dgamma_D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y);

    _dgamma_D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k2
(2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x))*bbn_ks_k2(x, y, z);

    _dgamma_D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k0
(2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x))*bbn_ks_k0(x, y, z);

    _dgamma_D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k0
(2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z))*bbn_ks_k0(x, y, z);

    _dgamma_D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z);

    _dgamma_D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k2
(2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y))*bbn_ks_k2(x, y, z);

    _dgamma_D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y) + bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y);

    _dgamma_D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x) + bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x);

    _dgamma_D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x);

    _dgamma_D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k1
(2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y))*bbn_ks_k1(x, y, z);

    _dgamma_D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z);

    _dgamma_D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k1
(2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x))*bbn_ks_k1(x, y, z);

    _dgamma_D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k1
(2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z))*bbn_ks_k1(x, y, z);

    _dgamma_D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k2
(2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z))*bbn_ks_k2(x, y, z);

  }
}

void bbn_free_data_ddg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params)
{

  const double r0          = Pgetd("BH_KerrSchild_RollOff");
  const double BH_center_x = Pgetd("BH_center_x");
  const double BH_center_y = Pgetd("BH_center_y");
  const double BH_center_z = Pgetd("BH_center_z");
  const double M_BH        = Pgetd("BH_irreducible_mass");
  const double a_BH        = Pgetd("BH_net_spin");
  const double chi_U0   = Pgetd("BH_chi_U0");
  const double chi_U1   = Pgetd("BH_chi_U1");
  const double chi_U2   = Pgetd("BH_chi_U2");
  const double y_CM = Pgetd("y_CM");
  const double x_CM = Pgetd("x_CM");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double chi = sqrt(Pow2(chi_U0)+Pow2(chi_U1)+Pow2(chi_U2));
  const unsigned nn = patch->nn;
  double phiy = 0,phiz = 0;
  double Bx,By,Bz,B2;/* B = v/c */
  unsigned ijk;

  assert(LSSEQL(chi,1));

  /* boost */
  Bx = -Omega_BHNS*(BH_center_y-y_CM);
  By =  Omega_BHNS*(BH_center_x-x_CM);
  Bz = Pgetd("BH_Vz");
  B2 = Pow2(Bx)+Pow2(By)+Pow2(Bz);

  /* rotation */
  if (!EQL(chi,0))/* otherwise tR is 0 */
  {
    phiz = -arctan(chi_U1,chi_U0);
    phiy = -acos(chi_U2/chi);
    assert(isfinite(phiy));
  }

double *const _ddgamma_D2D2D0D2 = get_v("_ddgamma_D2D2D0D2",params);
double *const _ddgamma_D0D1D0D1 = get_v("_ddgamma_D0D1D0D1",params);
double *const _ddgamma_D1D1D0D1 = get_v("_ddgamma_D1D1D0D1",params);
double *const _ddgamma_D1D2D0D1 = get_v("_ddgamma_D1D2D0D1",params);
double *const _ddgamma_D0D0D1D0 = get_v("_ddgamma_D0D0D1D0",params);
double *const _ddgamma_D1D1D2D0 = get_v("_ddgamma_D1D1D2D0",params);
double *const _ddgamma_D0D0D0D1 = get_v("_ddgamma_D0D0D0D1",params);
double *const _ddgamma_D1D1D0D0 = get_v("_ddgamma_D1D1D0D0",params);
double *const _ddgamma_D0D0D2D2 = get_v("_ddgamma_D0D0D2D2",params);
double *const _ddgamma_D2D2D1D1 = get_v("_ddgamma_D2D2D1D1",params);
double *const _ddgamma_D2D2D2D2 = get_v("_ddgamma_D2D2D2D2",params);
double *const _ddgamma_D0D0D2D0 = get_v("_ddgamma_D0D0D2D0",params);
double *const _ddgamma_D0D2D2D1 = get_v("_ddgamma_D0D2D2D1",params);
double *const _ddgamma_D1D2D0D2 = get_v("_ddgamma_D1D2D0D2",params);
double *const _ddgamma_D1D2D1D0 = get_v("_ddgamma_D1D2D1D0",params);
double *const _ddgamma_D2D2D2D0 = get_v("_ddgamma_D2D2D2D0",params);
double *const _ddgamma_D0D0D0D0 = get_v("_ddgamma_D0D0D0D0",params);
double *const _ddgamma_D0D0D1D1 = get_v("_ddgamma_D0D0D1D1",params);
double *const _ddgamma_D0D1D0D2 = get_v("_ddgamma_D0D1D0D2",params);
double *const _ddgamma_D1D2D1D1 = get_v("_ddgamma_D1D2D1D1",params);
double *const _ddgamma_D0D2D1D1 = get_v("_ddgamma_D0D2D1D1",params);
double *const _ddgamma_D1D2D2D2 = get_v("_ddgamma_D1D2D2D2",params);
double *const _ddgamma_D2D2D1D2 = get_v("_ddgamma_D2D2D1D2",params);
double *const _ddgamma_D0D1D2D0 = get_v("_ddgamma_D0D1D2D0",params);
double *const _ddgamma_D0D2D2D0 = get_v("_ddgamma_D0D2D2D0",params);
double *const _ddgamma_D0D1D1D0 = get_v("_ddgamma_D0D1D1D0",params);
double *const _ddgamma_D1D2D2D0 = get_v("_ddgamma_D1D2D2D0",params);
double *const _ddgamma_D1D1D1D0 = get_v("_ddgamma_D1D1D1D0",params);
double *const _ddgamma_D1D2D1D2 = get_v("_ddgamma_D1D2D1D2",params);
double *const _ddgamma_D1D2D0D0 = get_v("_ddgamma_D1D2D0D0",params);
double *const _ddgamma_D0D0D0D2 = get_v("_ddgamma_D0D0D0D2",params);
double *const _ddgamma_D2D2D2D1 = get_v("_ddgamma_D2D2D2D1",params);
double *const _ddgamma_D0D2D0D0 = get_v("_ddgamma_D0D2D0D0",params);
double *const _ddgamma_D0D2D0D1 = get_v("_ddgamma_D0D2D0D1",params);
double *const _ddgamma_D1D1D1D1 = get_v("_ddgamma_D1D1D1D1",params);
double *const _ddgamma_D0D1D0D0 = get_v("_ddgamma_D0D1D0D0",params);
double *const _ddgamma_D1D1D2D2 = get_v("_ddgamma_D1D1D2D2",params);
double *const _ddgamma_D0D2D1D2 = get_v("_ddgamma_D0D2D1D2",params);
double *const _ddgamma_D0D1D2D1 = get_v("_ddgamma_D0D1D2D1",params);
double *const _ddgamma_D0D1D1D2 = get_v("_ddgamma_D0D1D1D2",params);
double *const _ddgamma_D1D1D2D1 = get_v("_ddgamma_D1D1D2D1",params);
double *const _ddgamma_D0D1D1D1 = get_v("_ddgamma_D0D1D1D1",params);
double *const _ddgamma_D0D2D1D0 = get_v("_ddgamma_D0D2D1D0",params);
double *const _ddgamma_D2D2D1D0 = get_v("_ddgamma_D2D2D1D0",params);
double *const _ddgamma_D0D1D2D2 = get_v("_ddgamma_D0D1D2D2",params);
double *const _ddgamma_D2D2D0D0 = get_v("_ddgamma_D2D2D0D0",params);
double *const _ddgamma_D0D0D2D1 = get_v("_ddgamma_D0D0D2D1",params);
double *const _ddgamma_D0D2D2D2 = get_v("_ddgamma_D0D2D2D2",params);
double *const _ddgamma_D1D1D1D2 = get_v("_ddgamma_D1D1D1D2",params);
double *const _ddgamma_D1D2D2D1 = get_v("_ddgamma_D1D2D2D1",params);
double *const _ddgamma_D0D0D1D2 = get_v("_ddgamma_D0D0D1D2",params);
double *const _ddgamma_D0D2D0D2 = get_v("_ddgamma_D0D2D0D2",params);
double *const _ddgamma_D1D1D0D2 = get_v("_ddgamma_D1D1D0D2",params);
double *const _ddgamma_D2D2D0D1 = get_v("_ddgamma_D2D2D0D1",params);

  for (ijk = 0; ijk < nn; ++ijk) 
  {
    double x,y,z;
    x = patch->node[ijk]->x[0]-BH_center_x;
    y = patch->node[ijk]->x[1]-BH_center_y;
    z = patch->node[ijk]->x[2]-BH_center_z;

    _ddgamma_D2D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x);

    _ddgamma_D0D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D1D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x);

    _ddgamma_D1D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x);

    _ddgamma_D0D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D1D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x);

    _ddgamma_D0D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D1D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 2*
bbn_ks_c(x, y, z)*pow(Derivative(bbn_ks_k1(x, y, z), x), 2) + pow(bbn_ks_k1(x, y, z), 2)*
Derivative(bbn_ks_c(x, y, z), (x, 2)) + 4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), x);

    _ddgamma_D0D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*pow(Derivative(bbn_ks_k0(x, y, z), z), 2) + pow(bbn_ks_k0(x, y, z), 2)*
Derivative(bbn_ks_c(x, y, z), (z, 2)) + 4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), z);

    _ddgamma_D2D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*pow(Derivative(bbn_ks_k2(x, y, z), y), 2) + pow(bbn_ks_k2(x, y, z), 2)*
Derivative(bbn_ks_c(x, y, z), (y, 2)) + 4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y);

    _ddgamma_D2D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*pow(Derivative(bbn_ks_k2(x, y, z), z), 2) + pow(bbn_ks_k2(x, y, z), 2)*
Derivative(bbn_ks_c(x, y, z), (z, 2)) + 4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z);

    _ddgamma_D0D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D0D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y);

    _ddgamma_D1D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x);

    _ddgamma_D1D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x);

    _ddgamma_D2D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x);

    _ddgamma_D0D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 2*
bbn_ks_c(x, y, z)*pow(Derivative(bbn_ks_k0(x, y, z), x), 2) + pow(bbn_ks_k0(x, y, z), 2)*
Derivative(bbn_ks_c(x, y, z), (x, 2)) + 4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D0D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*pow(Derivative(bbn_ks_k0(x, y, z), y), 2) + pow(bbn_ks_k0(x, y, z), 2)*
Derivative(bbn_ks_c(x, y, z), (y, 2)) + 4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), y);

    _ddgamma_D0D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D1D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y);

    _ddgamma_D0D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y);

    _ddgamma_D1D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), z);

    _ddgamma_D2D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y);

    _ddgamma_D0D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D0D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D0D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D1D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x);

    _ddgamma_D1D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x);

    _ddgamma_D1D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y);

    _ddgamma_D1D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x);

    _ddgamma_D0D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D2D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y);

    _ddgamma_D0D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D0D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D1D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*pow(Derivative(bbn_ks_k1(x, y, z), y), 2) + pow(bbn_ks_k1(x, y, z), 2)*
Derivative(bbn_ks_c(x, y, z), (y, 2)) + 4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), y);

    _ddgamma_D0D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D1D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*pow(Derivative(bbn_ks_k1(x, y, z), z), 2) + pow(bbn_ks_k1(x, y, z), 2)*
Derivative(bbn_ks_c(x, y, z), (z, 2)) + 4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), z);

    _ddgamma_D0D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y);

    _ddgamma_D0D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y);

    _ddgamma_D0D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y);

    _ddgamma_D1D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y);

    _ddgamma_D0D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y);

    _ddgamma_D0D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D2D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _ddgamma_D0D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z);

    _ddgamma_D2D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 2*
bbn_ks_c(x, y, z)*pow(Derivative(bbn_ks_k2(x, y, z), x), 2) + pow(bbn_ks_k2(x, y, z), 2)*
Derivative(bbn_ks_c(x, y, z), (x, 2)) + 4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _ddgamma_D0D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y);

    _ddgamma_D0D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z);

    _ddgamma_D1D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y);

    _ddgamma_D1D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y);

    _ddgamma_D0D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y);

    _ddgamma_D0D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x);

    _ddgamma_D1D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x);

    _ddgamma_D2D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);
  }	
}

void bbn_free_data_dddg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params)
{
  const double r0          = Pgetd("BH_KerrSchild_RollOff");
  const double BH_center_x = Pgetd("BH_center_x");
  const double BH_center_y = Pgetd("BH_center_y");
  const double BH_center_z = Pgetd("BH_center_z");
  const double M_BH        = Pgetd("BH_irreducible_mass");
  const double a_BH        = Pgetd("BH_net_spin");
  const double chi_U0   = Pgetd("BH_chi_U0");
  const double chi_U1   = Pgetd("BH_chi_U1");
  const double chi_U2   = Pgetd("BH_chi_U2");
  const double y_CM = Pgetd("y_CM");
  const double x_CM = Pgetd("x_CM");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double chi = sqrt(Pow2(chi_U0)+Pow2(chi_U1)+Pow2(chi_U2));
  const unsigned nn = patch->nn;
  double phiy = 0,phiz = 0;
  double Bx,By,Bz,B2;/* B = v/c */
  unsigned ijk;

  assert(LSSEQL(chi,1));

  /* boost */
  Bx = -Omega_BHNS*(BH_center_y-y_CM);
  By =  Omega_BHNS*(BH_center_x-x_CM);
  Bz = Pgetd("BH_Vz");
  B2 = Pow2(Bx)+Pow2(By)+Pow2(Bz);

  /* rotation */
  if (!EQL(chi,0))/* otherwise tR is 0 */
  {
    phiz = -arctan(chi_U1,chi_U0);
    phiy = -acos(chi_U2/chi);
    assert(isfinite(phiy));
  }
double *const _dddgamma_D0D1D2D2D1 = get_v("_dddgamma_D0D1D2D2D1",params);
double *const _dddgamma_D0D1D2D0D1 = get_v("_dddgamma_D0D1D2D0D1",params);
double *const _dddgamma_D1D2D0D2D0 = get_v("_dddgamma_D1D2D0D2D0",params);
double *const _dddgamma_D1D1D1D2D1 = get_v("_dddgamma_D1D1D1D2D1",params);
double *const _dddgamma_D0D2D0D0D1 = get_v("_dddgamma_D0D2D0D0D1",params);
double *const _dddgamma_D0D2D0D2D2 = get_v("_dddgamma_D0D2D0D2D2",params);
double *const _dddgamma_D2D2D2D1D2 = get_v("_dddgamma_D2D2D2D1D2",params);
double *const _dddgamma_D0D1D0D2D0 = get_v("_dddgamma_D0D1D0D2D0",params);
double *const _dddgamma_D1D1D2D2D0 = get_v("_dddgamma_D1D1D2D2D0",params);
double *const _dddgamma_D0D2D0D2D1 = get_v("_dddgamma_D0D2D0D2D1",params);
double *const _dddgamma_D1D1D2D1D0 = get_v("_dddgamma_D1D1D2D1D0",params);
double *const _dddgamma_D1D1D0D1D2 = get_v("_dddgamma_D1D1D0D1D2",params);
double *const _dddgamma_D0D2D1D1D0 = get_v("_dddgamma_D0D2D1D1D0",params);
double *const _dddgamma_D1D1D1D0D1 = get_v("_dddgamma_D1D1D1D0D1",params);
double *const _dddgamma_D0D0D2D1D1 = get_v("_dddgamma_D0D0D2D1D1",params);
double *const _dddgamma_D0D1D2D1D2 = get_v("_dddgamma_D0D1D2D1D2",params);
double *const _dddgamma_D0D0D0D1D1 = get_v("_dddgamma_D0D0D0D1D1",params);
double *const _dddgamma_D2D2D1D2D1 = get_v("_dddgamma_D2D2D1D2D1",params);
double *const _dddgamma_D1D1D1D2D2 = get_v("_dddgamma_D1D1D1D2D2",params);
double *const _dddgamma_D2D2D2D2D1 = get_v("_dddgamma_D2D2D2D2D1",params);
double *const _dddgamma_D2D2D1D0D0 = get_v("_dddgamma_D2D2D1D0D0",params);
double *const _dddgamma_D0D1D2D1D1 = get_v("_dddgamma_D0D1D2D1D1",params);
double *const _dddgamma_D1D1D2D2D2 = get_v("_dddgamma_D1D1D2D2D2",params);
double *const _dddgamma_D1D2D1D0D1 = get_v("_dddgamma_D1D2D1D0D1",params);
double *const _dddgamma_D0D1D1D1D0 = get_v("_dddgamma_D0D1D1D1D0",params);
double *const _dddgamma_D2D2D2D0D1 = get_v("_dddgamma_D2D2D2D0D1",params);
double *const _dddgamma_D1D2D2D2D2 = get_v("_dddgamma_D1D2D2D2D2",params);
double *const _dddgamma_D1D2D1D2D1 = get_v("_dddgamma_D1D2D1D2D1",params);
double *const _dddgamma_D2D2D0D2D1 = get_v("_dddgamma_D2D2D0D2D1",params);
double *const _dddgamma_D1D1D2D1D2 = get_v("_dddgamma_D1D1D2D1D2",params);
double *const _dddgamma_D1D2D2D1D0 = get_v("_dddgamma_D1D2D2D1D0",params);
double *const _dddgamma_D1D2D0D0D1 = get_v("_dddgamma_D1D2D0D0D1",params);
double *const _dddgamma_D0D2D0D1D1 = get_v("_dddgamma_D0D2D0D1D1",params);
double *const _dddgamma_D0D2D1D0D1 = get_v("_dddgamma_D0D2D1D0D1",params);
double *const _dddgamma_D0D0D0D1D2 = get_v("_dddgamma_D0D0D0D1D2",params);
double *const _dddgamma_D0D0D2D1D2 = get_v("_dddgamma_D0D0D2D1D2",params);
double *const _dddgamma_D0D1D0D1D0 = get_v("_dddgamma_D0D1D0D1D0",params);
double *const _dddgamma_D1D2D1D0D0 = get_v("_dddgamma_D1D2D1D0D0",params);
double *const _dddgamma_D0D2D1D0D0 = get_v("_dddgamma_D0D2D1D0D0",params);
double *const _dddgamma_D1D2D0D1D1 = get_v("_dddgamma_D1D2D0D1D1",params);
double *const _dddgamma_D0D1D2D0D2 = get_v("_dddgamma_D0D1D2D0D2",params);
double *const _dddgamma_D1D1D2D0D0 = get_v("_dddgamma_D1D1D2D0D0",params);
double *const _dddgamma_D0D1D2D2D2 = get_v("_dddgamma_D0D1D2D2D2",params);
double *const _dddgamma_D0D0D2D0D2 = get_v("_dddgamma_D0D0D2D0D2",params);
double *const _dddgamma_D0D0D0D0D2 = get_v("_dddgamma_D0D0D0D0D2",params);
double *const _dddgamma_D1D1D0D1D1 = get_v("_dddgamma_D1D1D0D1D1",params);
double *const _dddgamma_D2D2D0D0D1 = get_v("_dddgamma_D2D2D0D0D1",params);
double *const _dddgamma_D0D0D1D2D2 = get_v("_dddgamma_D0D0D1D2D2",params);
double *const _dddgamma_D1D1D0D2D0 = get_v("_dddgamma_D1D1D0D2D0",params);
double *const _dddgamma_D0D1D1D2D0 = get_v("_dddgamma_D0D1D1D2D0",params);
double *const _dddgamma_D1D1D0D0D1 = get_v("_dddgamma_D1D1D0D0D1",params);
double *const _dddgamma_D1D1D0D0D2 = get_v("_dddgamma_D1D1D0D0D2",params);
double *const _dddgamma_D1D2D0D1D2 = get_v("_dddgamma_D1D2D0D1D2",params);
double *const _dddgamma_D0D0D1D0D2 = get_v("_dddgamma_D0D0D1D0D2",params);
double *const _dddgamma_D0D0D0D2D0 = get_v("_dddgamma_D0D0D0D2D0",params);
double *const _dddgamma_D2D2D2D1D1 = get_v("_dddgamma_D2D2D2D1D1",params);
double *const _dddgamma_D0D2D0D0D2 = get_v("_dddgamma_D0D2D0D0D2",params);
double *const _dddgamma_D2D2D2D2D0 = get_v("_dddgamma_D2D2D2D2D0",params);
double *const _dddgamma_D0D2D0D0D0 = get_v("_dddgamma_D0D2D0D0D0",params);
double *const _dddgamma_D2D2D1D2D0 = get_v("_dddgamma_D2D2D1D2D0",params);
double *const _dddgamma_D1D1D0D1D0 = get_v("_dddgamma_D1D1D0D1D0",params);
double *const _dddgamma_D1D2D2D0D2 = get_v("_dddgamma_D1D2D2D0D2",params);
double *const _dddgamma_D0D1D0D1D1 = get_v("_dddgamma_D0D1D0D1D1",params);
double *const _dddgamma_D2D2D0D2D2 = get_v("_dddgamma_D2D2D0D2D2",params);
double *const _dddgamma_D1D2D2D1D2 = get_v("_dddgamma_D1D2D2D1D2",params);
double *const _dddgamma_D2D2D1D1D1 = get_v("_dddgamma_D2D2D1D1D1",params);
double *const _dddgamma_D0D0D0D0D0 = get_v("_dddgamma_D0D0D0D0D0",params);
double *const _dddgamma_D0D0D1D2D1 = get_v("_dddgamma_D0D0D1D2D1",params);
double *const _dddgamma_D1D1D1D1D0 = get_v("_dddgamma_D1D1D1D1D0",params);
double *const _dddgamma_D1D1D1D1D1 = get_v("_dddgamma_D1D1D1D1D1",params);
double *const _dddgamma_D1D1D1D2D0 = get_v("_dddgamma_D1D1D1D2D0",params);
double *const _dddgamma_D0D1D0D2D1 = get_v("_dddgamma_D0D1D0D2D1",params);
double *const _dddgamma_D2D2D1D1D0 = get_v("_dddgamma_D2D2D1D1D0",params);
double *const _dddgamma_D1D2D2D2D0 = get_v("_dddgamma_D1D2D2D2D0",params);
double *const _dddgamma_D0D0D1D1D2 = get_v("_dddgamma_D0D0D1D1D2",params);
double *const _dddgamma_D0D2D0D2D0 = get_v("_dddgamma_D0D2D0D2D0",params);
double *const _dddgamma_D0D0D1D0D0 = get_v("_dddgamma_D0D0D1D0D0",params);
double *const _dddgamma_D0D2D2D1D0 = get_v("_dddgamma_D0D2D2D1D0",params);
double *const _dddgamma_D0D2D2D0D0 = get_v("_dddgamma_D0D2D2D0D0",params);
double *const _dddgamma_D1D1D1D1D2 = get_v("_dddgamma_D1D1D1D1D2",params);
double *const _dddgamma_D0D0D2D2D1 = get_v("_dddgamma_D0D0D2D2D1",params);
double *const _dddgamma_D2D2D0D1D0 = get_v("_dddgamma_D2D2D0D1D0",params);
double *const _dddgamma_D1D1D2D0D2 = get_v("_dddgamma_D1D1D2D0D2",params);
double *const _dddgamma_D0D1D1D1D1 = get_v("_dddgamma_D0D1D1D1D1",params);
double *const _dddgamma_D1D2D0D0D0 = get_v("_dddgamma_D1D2D0D0D0",params);
double *const _dddgamma_D1D2D0D0D2 = get_v("_dddgamma_D1D2D0D0D2",params);
double *const _dddgamma_D2D2D0D1D2 = get_v("_dddgamma_D2D2D0D1D2",params);
double *const _dddgamma_D0D0D1D1D1 = get_v("_dddgamma_D0D0D1D1D1",params);
double *const _dddgamma_D0D2D0D1D0 = get_v("_dddgamma_D0D2D0D1D0",params);
double *const _dddgamma_D0D2D1D2D1 = get_v("_dddgamma_D0D2D1D2D1",params);
double *const _dddgamma_D1D2D1D2D2 = get_v("_dddgamma_D1D2D1D2D2",params);
double *const _dddgamma_D0D1D1D0D0 = get_v("_dddgamma_D0D1D1D0D0",params);
double *const _dddgamma_D0D0D1D2D0 = get_v("_dddgamma_D0D0D1D2D0",params);
double *const _dddgamma_D1D2D1D1D0 = get_v("_dddgamma_D1D2D1D1D0",params);
double *const _dddgamma_D1D2D2D2D1 = get_v("_dddgamma_D1D2D2D2D1",params);
double *const _dddgamma_D0D1D0D0D0 = get_v("_dddgamma_D0D1D0D0D0",params);
double *const _dddgamma_D1D2D2D1D1 = get_v("_dddgamma_D1D2D2D1D1",params);
double *const _dddgamma_D0D0D0D1D0 = get_v("_dddgamma_D0D0D0D1D0",params);
double *const _dddgamma_D0D0D0D0D1 = get_v("_dddgamma_D0D0D0D0D1",params);
double *const _dddgamma_D0D0D0D2D2 = get_v("_dddgamma_D0D0D0D2D2",params);
double *const _dddgamma_D1D2D0D2D1 = get_v("_dddgamma_D1D2D0D2D1",params);
double *const _dddgamma_D2D2D1D1D2 = get_v("_dddgamma_D2D2D1D1D2",params);
double *const _dddgamma_D1D2D1D1D1 = get_v("_dddgamma_D1D2D1D1D1",params);
double *const _dddgamma_D0D0D1D0D1 = get_v("_dddgamma_D0D0D1D0D1",params);
double *const _dddgamma_D0D1D0D1D2 = get_v("_dddgamma_D0D1D0D1D2",params);
double *const _dddgamma_D0D2D1D2D2 = get_v("_dddgamma_D0D2D1D2D2",params);
double *const _dddgamma_D1D1D0D2D2 = get_v("_dddgamma_D1D1D0D2D2",params);
double *const _dddgamma_D0D1D0D0D2 = get_v("_dddgamma_D0D1D0D0D2",params);
double *const _dddgamma_D0D1D0D2D2 = get_v("_dddgamma_D0D1D0D2D2",params);
double *const _dddgamma_D0D2D1D2D0 = get_v("_dddgamma_D0D2D1D2D0",params);
double *const _dddgamma_D0D1D0D0D1 = get_v("_dddgamma_D0D1D0D0D1",params);
double *const _dddgamma_D1D2D1D1D2 = get_v("_dddgamma_D1D2D1D1D2",params);
double *const _dddgamma_D2D2D0D1D1 = get_v("_dddgamma_D2D2D0D1D1",params);
double *const _dddgamma_D1D1D2D2D1 = get_v("_dddgamma_D1D1D2D2D1",params);
double *const _dddgamma_D0D1D2D1D0 = get_v("_dddgamma_D0D1D2D1D0",params);
double *const _dddgamma_D0D2D1D0D2 = get_v("_dddgamma_D0D2D1D0D2",params);
double *const _dddgamma_D2D2D2D1D0 = get_v("_dddgamma_D2D2D2D1D0",params);
double *const _dddgamma_D0D1D1D1D2 = get_v("_dddgamma_D0D1D1D1D2",params);
double *const _dddgamma_D1D2D1D0D2 = get_v("_dddgamma_D1D2D1D0D2",params);
double *const _dddgamma_D0D0D2D2D0 = get_v("_dddgamma_D0D0D2D2D0",params);
double *const _dddgamma_D0D1D1D0D1 = get_v("_dddgamma_D0D1D1D0D1",params);
double *const _dddgamma_D0D1D2D2D0 = get_v("_dddgamma_D0D1D2D2D0",params);
double *const _dddgamma_D0D2D0D1D2 = get_v("_dddgamma_D0D2D0D1D2",params);
double *const _dddgamma_D0D2D2D2D1 = get_v("_dddgamma_D0D2D2D2D1",params);
double *const _dddgamma_D1D1D0D2D1 = get_v("_dddgamma_D1D1D0D2D1",params);
double *const _dddgamma_D0D0D2D1D0 = get_v("_dddgamma_D0D0D2D1D0",params);
double *const _dddgamma_D1D2D0D1D0 = get_v("_dddgamma_D1D2D0D1D0",params);
double *const _dddgamma_D0D2D2D1D2 = get_v("_dddgamma_D0D2D2D1D2",params);
double *const _dddgamma_D0D1D1D2D1 = get_v("_dddgamma_D0D1D1D2D1",params);
double *const _dddgamma_D1D2D2D0D0 = get_v("_dddgamma_D1D2D2D0D0",params);
double *const _dddgamma_D0D1D2D0D0 = get_v("_dddgamma_D0D1D2D0D0",params);
double *const _dddgamma_D1D2D1D2D0 = get_v("_dddgamma_D1D2D1D2D0",params);
double *const _dddgamma_D2D2D1D0D1 = get_v("_dddgamma_D2D2D1D0D1",params);
double *const _dddgamma_D2D2D1D2D2 = get_v("_dddgamma_D2D2D1D2D2",params);
double *const _dddgamma_D1D1D1D0D0 = get_v("_dddgamma_D1D1D1D0D0",params);
double *const _dddgamma_D0D1D1D0D2 = get_v("_dddgamma_D0D1D1D0D2",params);
double *const _dddgamma_D0D2D2D2D0 = get_v("_dddgamma_D0D2D2D2D0",params);
double *const _dddgamma_D2D2D2D2D2 = get_v("_dddgamma_D2D2D2D2D2",params);
double *const _dddgamma_D0D2D2D0D1 = get_v("_dddgamma_D0D2D2D0D1",params);
double *const _dddgamma_D0D2D1D1D2 = get_v("_dddgamma_D0D2D1D1D2",params);
double *const _dddgamma_D2D2D0D0D2 = get_v("_dddgamma_D2D2D0D0D2",params);
double *const _dddgamma_D2D2D2D0D0 = get_v("_dddgamma_D2D2D2D0D0",params);
double *const _dddgamma_D2D2D0D2D0 = get_v("_dddgamma_D2D2D0D2D0",params);
double *const _dddgamma_D1D1D2D0D1 = get_v("_dddgamma_D1D1D2D0D1",params);
double *const _dddgamma_D2D2D0D0D0 = get_v("_dddgamma_D2D2D0D0D0",params);
double *const _dddgamma_D2D2D2D0D2 = get_v("_dddgamma_D2D2D2D0D2",params);
double *const _dddgamma_D0D2D2D0D2 = get_v("_dddgamma_D0D2D2D0D2",params);
double *const _dddgamma_D0D2D2D2D2 = get_v("_dddgamma_D0D2D2D2D2",params);
double *const _dddgamma_D1D1D0D0D0 = get_v("_dddgamma_D1D1D0D0D0",params);
double *const _dddgamma_D1D2D2D0D1 = get_v("_dddgamma_D1D2D2D0D1",params);
double *const _dddgamma_D0D0D1D1D0 = get_v("_dddgamma_D0D0D1D1D0",params);
double *const _dddgamma_D1D1D1D0D2 = get_v("_dddgamma_D1D1D1D0D2",params);
double *const _dddgamma_D2D2D1D0D2 = get_v("_dddgamma_D2D2D1D0D2",params);
double *const _dddgamma_D0D2D1D1D1 = get_v("_dddgamma_D0D2D1D1D1",params);
double *const _dddgamma_D0D0D2D0D0 = get_v("_dddgamma_D0D0D2D0D0",params);
double *const _dddgamma_D0D1D1D2D2 = get_v("_dddgamma_D0D1D1D2D2",params);
double *const _dddgamma_D0D0D2D0D1 = get_v("_dddgamma_D0D0D2D0D1",params);
double *const _dddgamma_D0D0D0D2D1 = get_v("_dddgamma_D0D0D0D2D1",params);
double *const _dddgamma_D0D2D2D1D1 = get_v("_dddgamma_D0D2D2D1D1",params);
double *const _dddgamma_D0D0D2D2D2 = get_v("_dddgamma_D0D0D2D2D2",params);
double *const _dddgamma_D1D1D2D1D1 = get_v("_dddgamma_D1D1D2D1D1",params);
double *const _dddgamma_D1D2D0D2D2 = get_v("_dddgamma_D1D2D0D2D2",params);


  for (ijk = 0; ijk < nn; ++ijk) 
  {
    double x,y,z;
    x = patch->node[ijk]->x[0]-BH_center_x;
    y = patch->node[ijk]->x[1]-BH_center_y;
    z = patch->node[ijk]->x[2]-BH_center_z;

    _dddgamma_D0D1D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D0D1D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D1D2D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D1D1D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 4*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
4*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k1(x, y, z), y), 2);

    _dddgamma_D0D2D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D2D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D2D2D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 4*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k2(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z);

    _dddgamma_D0D1D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D1D1D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 4*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k1(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z);

    _dddgamma_D0D2D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D1D1D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D1D1D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D0D2D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D1D1D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 4*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k1(x, y, z), y), 2) + 4*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D0D0D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 4*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
4*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k0(x, y, z), y), 2);

    _dddgamma_D0D1D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D0D0D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 4*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k0(x, y, z), y), 2) + 4*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y);

    _dddgamma_D2D2D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 4*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
4*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k2(x, y, z), y), 2);

    _dddgamma_D1D1D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 4*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k1(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z);

    _dddgamma_D2D2D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 4*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k2(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z);

    _dddgamma_D2D2D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), y) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 4*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
pow(Derivative(bbn_ks_k2(x, y, z), x), 2);

    _dddgamma_D0D1D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D1D1D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 3)) + 6*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (z, 3)) + 6*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 6*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
6*Derivative(bbn_ks_c(x, y, z), z)*pow(Derivative(bbn_ks_k1(x, y, z), z), 2);

    _dddgamma_D1D2D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D1D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D2D2D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D1D2D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (z, 3)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 3)) + 3*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 3*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 3)) + 3*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 3*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
6*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z);

    _dddgamma_D1D2D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D2D2D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D1D1D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 4*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k1(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z);

    _dddgamma_D1D2D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D1D2D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D2D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D2D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D0D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y);

    _dddgamma_D0D0D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 4*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k0(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z);

    _dddgamma_D0D1D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), y) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D1D2D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D2D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D1D2D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D1D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D1D1D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 4*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k1(x, y, z), x), 2);

    _dddgamma_D0D1D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 3)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 3)) + 3*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 3*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 3)) + 3*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 3*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
3*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
3*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
6*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), z);

    _dddgamma_D0D0D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 4*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k0(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z);

    _dddgamma_D0D0D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 4*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k0(x, y, z), x), 2);

    _dddgamma_D1D1D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 4*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k1(x, y, z), y), 2) + 4*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D2D2D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), y) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 4*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
pow(Derivative(bbn_ks_k2(x, y, z), x), 2);

    _dddgamma_D0D0D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 4*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k0(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z);

    _dddgamma_D1D1D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 4*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k1(x, y, z), x), 2);

    _dddgamma_D0D1D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D1D1D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), y) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 4*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
pow(Derivative(bbn_ks_k1(x, y, z), x), 2);

    _dddgamma_D1D1D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 4*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k1(x, y, z), x), 2);

    _dddgamma_D1D2D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D0D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y);

    _dddgamma_D0D0D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 4*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k0(x, y, z), x), 2);

    _dddgamma_D2D2D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 4*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
4*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k2(x, y, z), y), 2);

    _dddgamma_D0D2D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D2D2D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 4*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k2(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z);

    _dddgamma_D0D2D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 3)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 3)) + 3*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 3*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 3)) + 3*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 3*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
6*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D2D2D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D1D1D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), y) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 4*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
pow(Derivative(bbn_ks_k1(x, y, z), x), 2);

    _dddgamma_D1D2D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D1D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D2D2D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 4*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k2(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z);

    _dddgamma_D1D2D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D2D2D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 3)) + 6*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 3)) + 6*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 6*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
6*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k2(x, y, z), y), 2);

    _dddgamma_D0D0D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 3)) + 6*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 3)) + 6*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 6*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
6*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k0(x, y, z), x), 2);

    _dddgamma_D0D0D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 4*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
4*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k0(x, y, z), y), 2);

    _dddgamma_D1D1D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 4*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k1(x, y, z), y), 2) + 4*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D1D1D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 3)) + 6*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 3)) + 6*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 6*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
6*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k1(x, y, z), y), 2);

    _dddgamma_D1D1D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D0D1D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D2D2D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 4*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k2(x, y, z), y), 2) + 4*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D1D2D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D0D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 4*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
4*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k0(x, y, z), y), 2);

    _dddgamma_D0D2D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D0D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), y) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 4*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
pow(Derivative(bbn_ks_k0(x, y, z), x), 2);

    _dddgamma_D0D2D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D2D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D1D1D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 4*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
4*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k1(x, y, z), y), 2);

    _dddgamma_D0D0D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 4*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k0(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z);

    _dddgamma_D2D2D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), y) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 4*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
pow(Derivative(bbn_ks_k2(x, y, z), x), 2);

    _dddgamma_D1D1D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 4*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k1(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z);

    _dddgamma_D0D1D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 3)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 3)) + 3*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 3*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 3)) + 3*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 3*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
3*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
3*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
6*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D1D2D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 3)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 3)) + 3*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 3*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 3)) + 3*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 3*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
6*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D1D2D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D2D2D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D0D0D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 3)) + 6*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 3)) + 6*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 6*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
6*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k0(x, y, z), y), 2);

    _dddgamma_D0D2D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D2D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D1D2D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D0D1D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), y) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D0D0D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y);

    _dddgamma_D1D2D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D1D2D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D0D1D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 3)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 3)) + 3*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 3*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 3)) + 3*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 3*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
3*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
3*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
6*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D1D2D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D0D0D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), y) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 4*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
pow(Derivative(bbn_ks_k0(x, y, z), x), 2);

    _dddgamma_D0D0D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), y) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 4*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
pow(Derivative(bbn_ks_k0(x, y, z), x), 2);

    _dddgamma_D0D0D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 4*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k0(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z);

    _dddgamma_D1D2D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D2D2D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 4*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
4*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k2(x, y, z), y), 2);

    _dddgamma_D1D2D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 3)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 3)) + 3*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 3*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 3)) + 3*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 3*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
6*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D0D0D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 4*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k0(x, y, z), y), 2) + 4*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y);

    _dddgamma_D0D1D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D0D2D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D1D1D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 4*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k1(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z);

    _dddgamma_D0D1D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D0D1D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D0D2D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D1D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), y) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D1D2D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D2D2D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 4*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k2(x, y, z), y), 2) + 4*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D1D1D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 4*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k1(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z);

    _dddgamma_D0D1D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D0D2D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D2D2D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D0D1D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D1D2D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D0D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 4*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k0(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z);

    _dddgamma_D0D1D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (y, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D0D1D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D0D2D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D2D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D1D1D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D0D0D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y);

    _dddgamma_D1D2D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), y) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), y) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D2D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D0D1D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D1D2D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D1D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D1D2D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D2D2D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 4*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k2(x, y, z), y), 2) + 4*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D2D2D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 4*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*pow(Derivative(bbn_ks_k2(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z);

    _dddgamma_D1D1D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2), y) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), y) + 4*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x, y) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, y) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + 2*Derivative(bbn_ks_c(x, y, z), y)*
pow(Derivative(bbn_ks_k1(x, y, z), x), 2);

    _dddgamma_D0D1D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x);

    _dddgamma_D0D2D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D2D2D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (z, 3)) + 6*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (z, 3)) + 6*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 6*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
6*Derivative(bbn_ks_c(x, y, z), z)*pow(Derivative(bbn_ks_k2(x, y, z), z), 2);

    _dddgamma_D0D2D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D2D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D2D2D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 4*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k2(x, y, z), x), 2);

    _dddgamma_D2D2D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 4*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k2(x, y, z), x), 2);

    _dddgamma_D2D2D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 4*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k2(x, y, z), x), 2);

    _dddgamma_D1D1D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D2D2D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (x, 3)) + 6*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 3)) + 6*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (x, 2)) + 6*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
6*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k2(x, y, z), x), 2);

    _dddgamma_D2D2D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 4*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
4*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k2(x, y, z), z), 2) + 4*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z);

    _dddgamma_D0D2D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D2D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (z, 3)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 3)) + 3*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 3*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 3)) + 3*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 3*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
6*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z);

    _dddgamma_D1D1D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (x, 3)) + 6*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 3)) + 6*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (x, 2)) + 6*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
6*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k1(x, y, z), x), 2);

    _dddgamma_D1D2D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), z) + Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x) + 
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y) + Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x);

    _dddgamma_D0D0D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, (y, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
4*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, (y, 2)) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 4*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), x) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*pow(Derivative(bbn_ks_k0(x, y, z), y), 2) + 4*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y);

    _dddgamma_D1D1D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D2D2D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
pow(bbn_ks_k2(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k2(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y, z) + 2*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, y) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D0D2D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 3)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 3)) + 3*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 3*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 3)) + 3*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 3*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), y) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
3*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
6*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D0D0D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (x, 2), z) + 4*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (x, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (x, 2)) + 
4*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), x, z) + 
4*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k0(x, y, z), x), 2);

    _dddgamma_D0D1D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y, (z, 2)) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), y) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), y);

    _dddgamma_D0D0D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y);

    _dddgamma_D0D0D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x, y, z) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), x, y, z) + 2*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y, z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x, y) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, y) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), x)*Derivative(bbn_ks_k0(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), x)*
Derivative(bbn_ks_k0(x, y, z), y);

    _dddgamma_D0D2D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), (y, 2), z) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (y, 2), z) + 2*bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k0(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 2*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k2(x, y, z), y, z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k2(x, y, z), z) + 
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), (y, 2)) + 
2*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k2(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (y, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k0(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
2*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), y)*
Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), y) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), y)*
Derivative(bbn_ks_k2(x, y, z), y);

    _dddgamma_D0D0D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*Derivative(bbn_ks_k0(x, y, z), (z, 3)) + 6*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k0(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 
pow(bbn_ks_k0(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (z, 3)) + 6*bbn_ks_k0(x, y, z)*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k0(x, y, z), (z, 2)) + 6*
bbn_ks_k0(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k0(x, y, z), z) + 
6*Derivative(bbn_ks_c(x, y, z), z)*pow(Derivative(bbn_ks_k0(x, y, z), z), 2);

    _dddgamma_D1D1D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2), z) + 4*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
pow(bbn_ks_k1(x, y, z), 2)*Derivative(bbn_ks_c(x, y, z), (y, 2), z) + 4*bbn_ks_k1(x, y, z)*
Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y, z) + 2*
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (y, 2))*Derivative(bbn_ks_k1(x, y, z), z) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), (y, 2)) + 
4*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k1(x, y, z), y)*Derivative(bbn_ks_c(x, y, z), y, z) + 
4*Derivative(bbn_ks_c(x, y, z), y)*Derivative(bbn_ks_k1(x, y, z), y)*
Derivative(bbn_ks_k1(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
pow(Derivative(bbn_ks_k1(x, y, z), y), 2);

    _dddgamma_D1D2D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), x, (z, 2)) + bbn_ks_c(x, y, z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 2*
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_c(x, y, z)*Derivative(bbn_ks_k1(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_c(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x, (z, 2)) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), (z, 2)) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k2(x, y, z), x, z) + 
bbn_ks_k1(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k2(x, y, z), x) + 
2*bbn_ks_k1(x, y, z)*Derivative(bbn_ks_k2(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), (z, 2)) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), x, z) + 
bbn_ks_k2(x, y, z)*Derivative(bbn_ks_c(x, y, z), (z, 2))*Derivative(bbn_ks_k1(x, y, z), x) + 
2*bbn_ks_k2(x, y, z)*Derivative(bbn_ks_k1(x, y, z), z)*Derivative(bbn_ks_c(x, y, z), x, z) + 
2*Derivative(bbn_ks_c(x, y, z), x)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), z) + 2*Derivative(bbn_ks_c(x, y, z), z)*
Derivative(bbn_ks_k1(x, y, z), x)*Derivative(bbn_ks_k2(x, y, z), z) + 2*
Derivative(bbn_ks_c(x, y, z), z)*Derivative(bbn_ks_k1(x, y, z), z)*
Derivative(bbn_ks_k2(x, y, z), x);
  }
}


