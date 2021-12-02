#ifndef maths_special_function_LIB_H
#define maths_special_function_LIB_H
#include "elliptica_system_lib.h"

#include <complex.h>
#undef I
#define imagI _Complex_I

/* Chebyshev polynomial of first kind Tn(x). x MUST be normalized value.
// ->return value: Tn(x)
*/
INLINE double Cheb_Tn(const int n, const double x)
{
  double t = DBL_MAX;
  
  if (n == 0)
    t = 1;
  else if (EQL(x,1))
    t = 1;
  else if (EQL(x,-1))
  {
    if (n%2)
      t = -1;
    else
      t = 1;
  }
  else
  {
    double th = acos(x);
    t = cos(n*th);
  }
  
  return t;
}

/* Chebyshev polynomial of second kind Un(x). x MUST be normalized value.
// ->return value: Un(x)
*/
INLINE double Cheb_Un(const int n, const double x)
{
  double u = DBL_MAX;
  
  if (n == 0) 
    u = 1;
  else if (EQL(x,1))
    u = n+1;
  else if (EQL(x,-1)) 
  {
    if (n%2)
      u = -n-1;
    else
      u = n+1;
  }  
  else
  {
    double th = acos(x);
    u = sin((n+1)*th)/sin(th);
  }
  
  return u;
}



double *c_f(Patch_T *const patch);
double *x_f(Patch_T *const patch);
double *y_f(Patch_T *const patch);
double *z_f(Patch_T *const patch);
double *poly5_f(Patch_T *const patch);
double *r_f(Patch_T *const patch);
double *inv_rP1_f(Patch_T *const patch);
double *cosx_f(Patch_T *const patch);
double *cosy_f(Patch_T *const patch);
double *cosz_f(Patch_T *const patch);
double *sinx_f(Patch_T *const patch);
double *siny_f(Patch_T *const patch);
double *sinz_f(Patch_T *const patch);
double *cosxyz_f(Patch_T *const patch);
double *cos4xyz_f(Patch_T *const patch);
double *cos5xyz_f(Patch_T *const patch);
double *sinxyz_f(Patch_T *const patch);
double *sin3xyz_f(Patch_T *const patch);
double *coshxyz_f(Patch_T *const patch);
double *sinhxyz_f(Patch_T *const patch);
double *tanhxyz_f(Patch_T *const patch);
double *logxyz_f(Patch_T *const patch);
double *mix1_f(Patch_T *const patch);
double *mix2_f(Patch_T *const patch);
double *c_f_x(Patch_T *const patch);
double *c_f_y(Patch_T *const patch);
double *c_f_z(Patch_T *const patch);
double *x_f_x(Patch_T *const patch);
double *x_f_y(Patch_T *const patch);
double *x_f_z(Patch_T *const patch);
double *x_f_xx(Patch_T *const patch);
double *y_f_x(Patch_T *const patch);
double *y_f_y(Patch_T *const patch);
double *y_f_z(Patch_T *const patch);
double *y_f_yy(Patch_T *const patch);
double *z_f_x(Patch_T *const patch);
double *z_f_y(Patch_T *const patch);
double *z_f_z(Patch_T *const patch);
double *z_f_zz(Patch_T *const patch);
double *r_f_x(Patch_T *const patch);
double *r_f_y(Patch_T *const patch);
double *r_f_z(Patch_T *const patch);
double *r_f_xx(Patch_T *const patch);
double *r_f_yy(Patch_T *const patch);
double *r_f_zz(Patch_T *const patch);
double *r_f_xy(Patch_T *const patch);
double *r_f_xz(Patch_T *const patch);
double *r_f_yz(Patch_T *const patch);
double *r_f_xyz(Patch_T *const patch);
double *inv_rP1_f_x(Patch_T *const patch);
double *inv_rP1_f_y(Patch_T *const patch);
double *inv_rP1_f_z(Patch_T *const patch);
double *inv_rP1_f_xx(Patch_T *const patch);
double *inv_rP1_f_yy(Patch_T *const patch);
double *inv_rP1_f_zz(Patch_T *const patch);
double *inv_rP1_f_xy(Patch_T *const patch);
double *inv_rP1_f_xz(Patch_T *const patch);
double *inv_rP1_f_yz(Patch_T *const patch);
double *inv_rP1_f_xyz(Patch_T *const patch);
double *sinxyz_f_x(Patch_T *const patch);
double *sinxyz_f_y(Patch_T *const patch);
double *sinxyz_f_z(Patch_T *const patch);
double *sinxyz_f_xx(Patch_T *const patch);
double *sinxyz_f_yy(Patch_T *const patch);
double *sinxyz_f_zz(Patch_T *const patch);
double *sinxyz_f_xy(Patch_T *const patch);
double *sinxyz_f_xz(Patch_T *const patch);
double *sinxyz_f_yz(Patch_T *const patch);
double *sinxyz_f_xyz(Patch_T *const patch);
double *poly5_f_x(Patch_T *const patch);
double *poly5_f_y(Patch_T *const patch);
double *poly5_f_z(Patch_T *const patch);
double *poly5_f_xx(Patch_T *const patch);
double *poly5_f_yy(Patch_T *const patch);
double *poly5_f_zz(Patch_T *const patch);
double *poly5_f_xy(Patch_T *const patch);
double *poly5_f_xz(Patch_T *const patch);
double *poly5_f_yz(Patch_T *const patch);
double *poly5_f_xyz(Patch_T *const patch);
double *mix2_f_x(Patch_T *const patch);
double *mix2_f_y(Patch_T *const patch);
double *mix2_f_z(Patch_T *const patch);
double *mix2_f_xx(Patch_T *const patch);
double *mix2_f_yy(Patch_T *const patch);
double *mix2_f_zz(Patch_T *const patch);
double *mix2_f_xy(Patch_T *const patch);
double *mix2_f_xz(Patch_T *const patch);
double *mix2_f_yz(Patch_T *const patch);
double *mix2_f_xyz(Patch_T *const patch);
double *sinx_f_x(Patch_T *const patch);
double *sinx_f_y(Patch_T *const patch);
double *sinx_f_z(Patch_T *const patch);
double *sinx_f_xx(Patch_T *const patch);
double *sinx_f_yy(Patch_T *const patch);
double *sinx_f_zz(Patch_T *const patch);
double *sinx_f_xy(Patch_T *const patch);
double *sinx_f_xz(Patch_T *const patch);
double *sinx_f_yz(Patch_T *const patch);
double *sinx_f_xyz(Patch_T *const patch);
double poly5_f_point(const double x,const double y,const double z);
int Factorial(const int n);
INLINE double Cheb_Tn(const int n, const double x) INLINE_WARN_UNUSED_FUNC;
INLINE double Cheb_Un(const int n, const double x) INLINE_WARN_UNUSED_FUNC;
double d2Cheb_Tn_dx2(const int n, const double x);
double dCheb_Tn_dx(const int n, const double x);
double d2Cheb_Tn_dx2(const int n, const double x);

double 
d_associated_legendre_dtheta
  (
  const int l/* l in P^l_m(x) */,
  const int m/* m in P^l_m(x), 0 <= m <= l*/,
  const double x/* x in P^l_m(x), -1 <= x=cos(theta) <= 1 */
  );

double 
associated_legendre
  (
  const int l/* l in P^l_m(x) */,
  const int m/* m in P^l_m(x), 0 <= m <= l*/,
  const double x/* x in P^l_m(x), -1 <= x=cos(theta) <= 1 */
  );

double complex 
Ylm
  (
  const int l/* l in Y_l^m(theta,phi) */, 
  const int m/* m in Y_l^m(theta,phi), -l <= m <= l */, 
  const double theta/* theta in Y_l^m(theta,phi) */,
  const double phi/*  phi in Y_l^m(theta,phi) */
  );

double complex 
dYlm_dphi
  (
  const int l/* l in Y_l^m(theta,phi) */, 
  const int m/* m in Y_l^m(theta,phi), -l <= m <= l */, 
  const double theta/* theta in Y_l^m(theta,phi) */,
  const double phi/*  phi in Y_l^m(theta,phi) */
  );

double complex 
dYlm_dtheta
  (
  const int l/* l in Y_l^m(theta,phi) */, 
  const int m/* m in Y_l^m(theta,phi), -l <= m <= l */, 
  const double theta/* theta in Y_l^m(theta,phi) */,
  const double phi/*  phi in Y_l^m(theta,phi) */
  );


#endif



