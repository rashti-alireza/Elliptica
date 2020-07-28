#ifndef maths_approximation_LIB_H
#define maths_approximation_LIB_H


struct INTERPOLATION_T;
struct FIELD_T;


/* interpolation function typedef */
typedef double fInterpolation_T(struct INTERPOLATION_T *const interp_s);

/* interpolation struct used in interpolation function */
typedef struct INTERPOLATION_T
{
  const char *method;
  struct FIELD_T *field;/* interesting field for interpolation */
  fInterpolation_T *interpolation_func;/* interpolation function */
  double X,Y,Z;/* where interpolant calculated. 
               // MUST be provided in coords sys. used by patch.
               */
  unsigned X_dir_flag   : 1;/* 1-D interpolation in X direction */
  unsigned Y_dir_flag   : 1;/* 1-D interpolation in Y direction */
  unsigned Z_dir_flag   : 1;/* 1-D interpolation in Z direction */
  unsigned XY_dir_flag  : 1;/* 2-D interpolation in X&Y direction */
  unsigned XZ_dir_flag  : 1;/* 2-D interpolation in X&Z direction */
  unsigned YZ_dir_flag  : 1;/* 2-D interpolation in Y&Z direction */
  unsigned XYZ_dir_flag : 1;/* 3-D interpolation in X&Y&Z direction */
  unsigned I;/* the index held constant in case of interpolation in 1-D and 2-D */
  unsigned J;/* the index held constant in case of interpolation in 1-D and 2-D */
  unsigned K;/* the index held constant in case of interpolation in 1-D and 2-D */
  struct
  {
   const double *f;/* f(xi)'s */
   double *x;/* xi's */
   double h;/* desired point to interpolate f */
   unsigned N;/* number of xi's */
   unsigned max;/* desired number of xi to be used, if 0 then the value N is picked */
  }Neville_1d[1];/* the method is Neville's iterated interpolation */
  struct
  {
   double *f;/* f(xi)'s */
   double *x;/* xi's, note: it must be x0 < x1 < ...< xN */
   double h;/* desired point to interpolate f */
   unsigned N;/* number of xi's */
   double *a,*b,*c,*d;/* coefficents in s(h) = a+b(h-xi)+c(h-xi)^2+d(h-xi)^3 */
   unsigned Order: 1;/* if xi's in the order 1, otherwise 0 */
   unsigned Alloc_Mem: 1;/* if it allocates memory for x and f */
  }N_cubic_spline_1d[1];/* natural cubic spline 1d */
}Interpolation_T;

void rft_1d_ChebyshevExtrema_coeffs(double *const values ,double *const coeffs,const unsigned n);
void rft_1d_ChebyshevNodes_coeffs(double *const values ,double *const coeffs,const unsigned n);
void *r2cft_1d_EquiSpaced_coeffs(const double *const value,const unsigned n);
double *c2rft_1d_EquiSpaced_values(void *const coeffs,const unsigned N);
int derivative_tests(Grid_T *const grid);
int interpolation_tests(Grid_T *const grid);
int fourier_transformation_tests(Grid_T *const grid);
int Ylm_transformation_tests(Grid_T *const grid);
Interpolation_T *init_interpolation(void);
double execute_interpolation(Interpolation_T *const interp_s);
void plan_interpolation(Interpolation_T *const interp_s);
void get_Ylm_coeffs(double *const realClm,double *const imagClm,const double *const f,const unsigned Ntheta,const unsigned Nphi,const unsigned Lmax);
double interpolation_Ylm(const double *const realClm,const double *const imagClm,const unsigned Lmax, const double theta, const double phi);
double *df_dphi_Ylm(const double *const realClm,const double *const imagClm,const unsigned Ntheta, const unsigned Nphi,const unsigned Lmax);
double *df_dtheta_Ylm(const double *const realClm,const double *const imagClm,const unsigned Ntheta, const unsigned Nphi,const unsigned Lmax);
unsigned lm2n(const unsigned l, const unsigned m);
double *alloc_ClmYlm(unsigned Lmax);
void free_interpolation(Interpolation_T *interp_s);

double *
r2cft_2d_df_dphi1
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const unsigned Nphi0/* number of point in phi0 direction */,
  const unsigned Nphi1/* number of point in phi1 direction */
);

double 
r2cft_2d_interpolation
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const unsigned Nphi0/* number of point in phi0 direction */,
  const unsigned Nphi1/* number of point in phi1 direction */,
  const double phi0/* point of interest at phi0 dir */,
  const double phi1/* point of interest at phi0 dir */
);

double *
r2cft_2d_df_dphi0
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const unsigned Nphi0/* number of point in phi0 direction */,
  const unsigned Nphi1/* number of point in phi1 direction */
);

void
r2cft_2d_coeffs
(
  const double *const f/* field values */,
  const unsigned Nphi0/* number of point in phi0 direction */, 
  const unsigned Nphi1/* number of point in phi1 direction */,
  double **const realC/* real part of coeffs, allocates memory */,
  double **const imagC/* imag part of coeffs, allocates memory*/
);

void
r2cft_2d_coeffs_S2
(
  const double *const f/* field values given on theta and phi coords. */,
  const unsigned Ntheta/* number of point in theta direction */, 
  const unsigned Nphi/* number of point in phi direction */,
  double **const realC/* real part of coeffs, allocates memory */,
  double **const imagC/* imag part of coeffs, allocates memory*/
);

double 
r2cft_2d_interpolation_S2
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const unsigned Ntheta/* number of point in theta direction */,
  const unsigned Nphi/* number of point in phi direction */,
  const double theta/* point of interest at theta dir */,
  const double phi/* point of interest at phi dir */
);

double *
r2cft_2d_df_dphi_S2
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const unsigned Ntheta/* number of point in theta direction */,
  const unsigned Nphi/* number of point in phi direction */
);

double *
r2cft_2d_df_dtheta_S2
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const unsigned Ntheta/* number of point in theta direction */,
  const unsigned Nphi/* number of point in phi direction */
);


#endif



