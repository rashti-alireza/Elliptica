#ifndef maths_analytic_LIB_H
#define maths_analytic_LIB_H


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
double associated_legendre(const int l, const int m, const double x);
void init_associated_legendre(void);
void init_Ylm(void);
void init_dYlm_dphi(void);
void init_dYlm_dtheta(void);

#endif


