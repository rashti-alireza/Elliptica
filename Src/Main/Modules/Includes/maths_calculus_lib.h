#ifndef maths_calculus_LIB_H
#define maths_calculus_LIB_H
#include "elliptica_system_lib.h"


struct FIELD_T;

/* struct fot integration */
typedef struct INTEGRATION_T
{
 const char *type;
 double err;/* expected error for this method */
 /* metrices components used in surface and volume integrals.
 // note: the metric is assumed to be symmetric */
 const double *g00;
 const double *g01;
 const double *g02;
 const double *g11;
 const double *g12;
 const double *g22;
 struct
 {
  const struct FIELD_T *f;/* integrand */
  Uint X_surface : 1;/* integration on the hyper-surface X = const */
  Uint Y_surface : 1;/* integration on the hyper-surface Y = const */
  Uint Z_surface : 1;/* integration on the hyper-surface Z = const */
  Uint I;/* the index showing the hyper-surface X = const */
  Uint J;/* the index showing the hyper-surface Y = const */
  Uint K;/* the index showing the hyper-surface Z = const */
  Uint Ii,If;/* i in [Ii,If] for X-direction(inclusive), used in vol integral */
  Uint Ji,Jf;/* j in [Ji,Jf] for Y-direction(inclusive), used in vol integral */
  Uint Ki,Kf;/* k in [Ki,Kf] for K-direction(inclusive), used in vol integral */
 }Spectral[1];
 struct
 {
  double a,b;/* integral limits */
  const double *f;/* integrand */
  Uint n;/* odd positive integer  */
 }Composite_Simpson_1D[1];
 struct
 {
   const double *f;/* integrand in \integral_{-1}^{1} f(x)/sqrt(1-x^2) dx */
   Uint n;/* number of collocation points for function f */
 }GQ_ChebyshevExtrema[1];/* Gauss Quadrature using chebyshev Extrema collocation */
 struct
 {
   const double *f;/* integrand in \integral_{-1}^{1} f(x) dx */
   Uint n;/* number of collocation points for function f */
 }GQ_Lobatto[1];/* Gauss Quadrature using Lobatto's collocation */
 struct
 {
   const double *f;/* integrand in \integral_{-1}^{1} f(x) dx */
   Uint n;/* number of collocation points for function f */
 }GQ_Legendre[1];/* Gauss Quadrature using Legendre's collocation */
 double (*integration_func)(struct INTEGRATION_T *const I);/* function that integrates */
}Integration_T;

double *Partial_Derivative(struct FIELD_T *const f,const char *task);
double *partial_derivative(struct FIELD_T *const dfield);
void free_integration(Integration_T *I);
void plan_integration(Integration_T *const I);
double execute_integration(Integration_T *const I);
Integration_T *init_integration(void);
int integration_tests(Grid_T *const grid);
double Lobatto_weight_function(const double x, const Uint N);
double Lobatto_root_function(const Uint rootN, const Uint N);
void init_Lobatto_root_function(void);
double Legendre_root_function(const Uint rootN, const Uint N);
double Legendre_weight_function(const double x, const Uint N);
double dLegendre_dx(const Uint n, const double x);
void init_Legendre_root_function(void);
void init_dLegendre_dx(void);
void partial_derivative_regex(Patch_T *const patch,
                                   const char *const regex_list);

#endif


