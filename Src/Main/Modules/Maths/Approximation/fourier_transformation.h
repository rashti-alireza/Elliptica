#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_complex_numbers_lib.h"

#define IJ(i,j,n)  ((j)+(i)*(n))

void rft_1d_ChebyshevExtrema_coeffs(double *const values ,double *const coeffs,const unsigned n);
void rft_1d_ChebyshevNodes_coeffs(double *const values ,double *const coeffs,const unsigned n);
void *r2cft_1d_EquiSpaced_coeffs(const double *const value,const unsigned n);
double *c2rft_1d_EquiSpaced_values(void *const coeffs,const unsigned N);
void *r2cft_2d_Coeffs(const double *const f,const unsigned Nphi0, const unsigned Nphi1);
double *r2cft_2d_realCs(void *C,const unsigned Nphi0, const unsigned Nphi1);
double *r2cft_2d_ImagCs(void *C,const unsigned Nphi0, const unsigned Nphi1);
double r2cft_2d_interpolation(void *C,const unsigned Nphi0, const unsigned Nphi1,const double phi0,const double phi1);
double *r2cft_2d_df_dphi0(void *C,const unsigned Nphi0, const unsigned Nphi1);
double *r2cft_2d_df_dphi1(void *C,const unsigned Nphi0, const unsigned Nphi1);



