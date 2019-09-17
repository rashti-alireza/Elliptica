/* struct fot integration */
typedef struct INTEGRATION_T
{
 const char *type;
 double err;/* expected error for this method */
 struct
 {
  Field_T *f;/* integrand */
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

double *Partial_Derivative(Field_T *const f,const char *task);
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
