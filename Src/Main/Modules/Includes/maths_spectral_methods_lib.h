#ifndef maths_spectral_methods_LIB_H
#define maths_spectral_methods_LIB_H
#include "elliptica_system_lib.h"

/* 1d array index for Ylm */
#define IJ_Ylm(i,j,Nphi) ((j)+(i)*(Nphi))
/* Nphi at Ylm */
#define Nphi_Ylm(lmax) (2*(lmax)+1)
/* Ntheta at Ylm */
#define Ntheta_Ylm(lmax) (2*(lmax)+1)
/* Ntheta*Nphi at Ylm */
#define Ntotal_Ylm(lmax) (Ntheta_Ylm(lmax)*Nphi_Ylm(lmax))
/* number of coeffs */
#define Ncoeffs_Ylm(lmax) (((lmax)+1)*(lmax)/2 + (lmax)+1)

struct INTERPOLATION_T;
struct FIELD_T;
struct PATCH_T;

/* interpolation function typedef */
typedef double fInterpolation_T(struct INTERPOLATION_T *const interp_s);

/* interpolation struct used in interpolation function */
typedef struct INTERPOLATION_T
{
  const char *method;
  fInterpolation_T *interpolation_func;/* interpolation function (interpolant) */
  fInterpolation_T *interpolation_1st_deriv;/* first derivative of the interpolant(if available) */

  //////////////
  // spectral //
  //////////////
  struct FIELD_T *field;/* interesting field for interpolation */
  double X,Y,Z;/* where interpolant calculated. 
               // MUST be provided in coords sys. used by patch.
               */
  Uint X_dir_flag   : 1;/* 1-D interpolation in X direction */
  Uint Y_dir_flag   : 1;/* 1-D interpolation in Y direction */
  Uint Z_dir_flag   : 1;/* 1-D interpolation in Z direction */
  Uint XY_dir_flag  : 1;/* 2-D interpolation in X&Y direction */
  Uint XZ_dir_flag  : 1;/* 2-D interpolation in X&Z direction */
  Uint YZ_dir_flag  : 1;/* 2-D interpolation in Y&Z direction */
  Uint XYZ_dir_flag : 1;/* 3-D interpolation in X&Y&Z direction */
  Uint I;/* the index held constant in case of interpolation in 1-D and 2-D */
  Uint J;/* the index held constant in case of interpolation in 1-D and 2-D */
  Uint K;/* the index held constant in case of interpolation in 1-D and 2-D */
  
  /////////////
  // Neville //
  /////////////
  struct
  {
   const double *f;/* f(xi)'s */
   double *x;/* xi's */
   double h;/* desired point to interpolate f */
   Uint N;/* number of xi's */
   Uint max;/* desired number of xi to be used, if 0 then the value N is picked */
  }Neville_1d[1];/* the method is Neville's iterated interpolation */
  
  //////////////////
  // cubic spline //
  //////////////////
  struct
  {
   double *f;/* f(xi)'s */
   double *x;/* xi's, note: it must be x0 < x1 < ...< xN */
   double h;/* desired point to interpolate f */
   Uint N;/* number of xi's */
   double *a,*b,*c,*d;/* coefficents in s(h) = a+b(h-xi)+c(h-xi)^2+d(h-xi)^3 */
   Uint Order: 1;/* if xi's in the order 1, otherwise 0 */
   Uint Alloc_Mem: 1;/* if it allocates memory for x and f */
   Uint No_Warn: 1;/* if 1 it prints NO warning in case of an error */
  }N_cubic_spline_1d[1];/* natural cubic spline 1d */
  
  ////////////////////
  // Hermite spline //
  ////////////////////
  struct
  { 
   double *f;// f(xi)
   double *fp;// df(xi)/dx
   double *x;// coordinate grid
   double h;// point to interpolate
   Uint N;// number of grid points
   Uint fd_accuracy_order;// order of finite difference accuracy
   Uint fd_derivative_order;// order of derivative required from finite difference method
   Uint num_points;// the number of points being used for the interpolant, 
                   // namely, polynomial degree = (2*num_points-1)
   Uint Order: 1;// 1 if x array in order
   Uint Alloc_fx: 1;// if f and x ordered 1
   Uint Alloc_fp: 1;// if fp allocated 1
   Uint No_Warn: 1;// emit warning if 0
  }Hermite_1d[1];
}Interpolation_T;

void rft_1d_ChebyshevExtrema_coeffs(double *const values ,double *const coeffs,const Uint n);
void rft_1d_ChebyshevNodes_coeffs(double *const values ,double *const coeffs,const Uint n);
void *r2cft_1d_EquiSpaced_coeffs(const double *const value,const Uint n);
double *c2rft_1d_EquiSpaced_values(void *const coeffs,const Uint N);
int derivative_tests(Grid_T *const grid);
int interpolation_tests(Grid_T *const grid);
int fourier_transformation_tests(Grid_T *const grid);
int Ylm_transformation_tests(Grid_T *const grid);
Interpolation_T *init_interpolation(void);
double execute_interpolation(Interpolation_T *const interp_s);
double execute_1st_deriv_interpolation(Interpolation_T *const interp_s);
void plan_interpolation(Interpolation_T *const interp_s);
void get_Ylm_coeffs(double *const realClm,double *const imagClm,const double *const f,const Uint Ntheta,const Uint Nphi,const Uint Lmax);
double interpolation_Ylm(const double *const realClm,const double *const imagClm,const Uint Lmax, const double theta, const double phi);
double *df_dphi_Ylm(const double *const realClm,const double *const imagClm,const Uint Ntheta, const Uint Nphi,const Uint Lmax);
double *df_dtheta_Ylm(const double *const realClm,const double *const imagClm,const Uint Ntheta, const Uint Nphi,const Uint Lmax);
Uint lm2n(const Uint l, const Uint m);
double *alloc_ClmYlm(Uint Lmax);
void free_interpolation(Interpolation_T *interp_s);

double *
r2cft_2d_df_dphi1
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Nphi0/* number of point in phi0 direction */,
  const Uint Nphi1/* number of point in phi1 direction */
);

double 
r2cft_2d_interpolation
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Nphi0/* number of point in phi0 direction */,
  const Uint Nphi1/* number of point in phi1 direction */,
  const double phi0/* point of interest at phi0 dir */,
  const double phi1/* point of interest at phi0 dir */
);

double *
r2cft_2d_df_dphi0
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Nphi0/* number of point in phi0 direction */,
  const Uint Nphi1/* number of point in phi1 direction */
);

void
r2cft_2d_coeffs
(
  const double *const f/* field values */,
  const Uint Nphi0/* number of point in phi0 direction */, 
  const Uint Nphi1/* number of point in phi1 direction */,
  double **const realC/* real part of coeffs, allocates memory */,
  double **const imagC/* imag part of coeffs, allocates memory*/
);

double
r2cft_2d_coeffs_S2
(
  const double *const f/* field values given on theta and phi coords. */,
  Uint Ntheta/* number of point in theta direction */, 
  const Uint Nphi/* number of point in phi direction */,
  double **const realC/* real part of coeffs, allocates memory */,
  double **const imagC/* imag part of coeffs, allocates memory*/,
  const int improve/* if 1, it tries to improve the expansion, otherwise no. */
);

double 
r2cft_2d_interpolation_S2
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Ntheta/* number of point in theta direction */,
  const Uint Nphi/* number of point in phi direction */,
  const double theta/* point of interest at theta dir */,
  const double phi/* point of interest at phi dir */
);

double *
r2cft_2d_df_dphi_S2
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Ntheta/* number of point in theta direction */,
  const Uint Nphi/* number of point in phi direction */
);

double *
r2cft_2d_df_dtheta_S2
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Ntheta/* number of point in theta direction */,
  const Uint Nphi/* number of point in phi direction */
);

// encapsulate filter args.
typedef struct SPECTRAL_FILTER_T
{
 struct PATCH_T *patch;// patch that has the field
 const char *field;// field name, e.g., "rho0".
 const char *filter;// name of the filter, e.g., erfclog
 int erfclog_p;// p arg for erfclog filter.
}spectral_filter_T;

void spectral_filter(const spectral_filter_T *const args);

#endif



