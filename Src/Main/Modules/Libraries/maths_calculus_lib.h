#ifndef maths_calculus_LIB_H
#define maths_calculus_LIB_H


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
  unsigned X_surface : 1;/* integration on the hyper-surface X = const */
  unsigned Y_surface : 1;/* integration on the hyper-surface Y = const */
  unsigned Z_surface : 1;/* integration on the hyper-surface Z = const */
  unsigned I;/* the index showing the hyper-surface X = const */
  unsigned J;/* the index showing the hyper-surface Y = const */
  unsigned K;/* the index showing the hyper-surface Z = const */
 }Spectral[1];
 struct
 {
  double a,b;/* integral limits */
  const double *f;/* integrand */
  unsigned n;/* odd positive integer  */
 }Composite_Simpson_1D[1];
 struct
 {
   const double *f;/* integrand in \integral_{-1}^{1} f(x)/sqrt(1-x^2) dx */
   unsigned n;/* number of collocation points for function f */
 }GQ_ChebyshevExtrema[1];/* Gauss Quadrature using chebyshev Extrema collocation */
 struct
 {
   const double *f;/* integrand in \integral_{-1}^{1} f(x) dx */
   unsigned n;/* number of collocation points for function f */
 }GQ_Lobatto[1];/* Gauss Quadrature using Lobatto's collocation */
 struct
 {
   const double *f;/* integrand in \integral_{-1}^{1} f(x) dx */
   unsigned n;/* number of collocation points for function f */
 }GQ_Legendre[1];/* Gauss Quadrature using Legendre's collocation */
 double (*integration_func)(struct INTEGRATION_T *const I);/* function that integrates */
}Integration_T;

double *Partial_Derivative(struct FIELD_T *const f,const char *task);
void free_integration(Integration_T *I);
void plan_integration(Integration_T *const I);
double execute_integration(Integration_T *const I);
Integration_T *init_integration(void);
int integration_tests(Grid_T *const grid);
double Lobatto_weight_function(const double x, const unsigned N);
double Lobatto_root_function(const unsigned rootN, const unsigned N);
void init_Lobatto_root_function(void);
double Legendre_root_function(const unsigned rootN, const unsigned N);
double Legendre_weight_function(const double x, const unsigned N);
double dLegendre_dx(const unsigned n, const double x);
void init_Legendre_root_function(void);
void init_dLegendre_dx(void);
double Integrate_ChebTn(const unsigned n,const double xi,const double xf);


#endif


