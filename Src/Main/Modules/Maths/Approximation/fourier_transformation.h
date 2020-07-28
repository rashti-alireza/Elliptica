#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_complex_numbers_lib.h"

#define IJ(i,j,n)  ((j)+(i)*(n))

void rft_1d_ChebyshevExtrema_coeffs(double *const values ,double *const coeffs,const unsigned n);
void rft_1d_ChebyshevNodes_coeffs(double *const values ,double *const coeffs,const unsigned n);
void *r2cft_1d_EquiSpaced_coeffs(const double *const value,const unsigned n);
double *c2rft_1d_EquiSpaced_values(void *const coeffs,const unsigned N);

void
r2cft_2d_coeffs
(
  const double *const f/* field values */,
  const unsigned Nphi0/* number of point in phi0 direction */, 
  const unsigned Nphi1/* number of point in phi1 direction */,
  double **const realC/* real part of coeffs, allocates memory */,
  double **const imagC/* imag part of coeffs, allocates memory*/
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

double *
r2cft_2d_df_dphi1
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const unsigned Nphi0/* number of point in phi0 direction */,
  const unsigned Nphi1/* number of point in phi1 direction */
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



