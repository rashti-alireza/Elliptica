#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "maths_special_functions_lib.h"
#include "utilities_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"

#define DO 1
#define DO_NOT 0

#define E M_E
#define Cos(a) cos(a)
#define Sin(a) sin(a)
#define Cosh(a) cosh(a)
#define Sinh(a) sinh(a)
#define Log(a) log(a)
#define Power(a,b) pow(a,b)
#define Sqrt(a) sqrt(a)

int integration_tests(Grid_T *const grid);
void free_integration(Integration_T *I);
void plan_integration(Integration_T *const I);
double execute_integration(Integration_T *const I);
Integration_T *init_integration(void);
static int csr_1d(Grid_T *const grid);
static int GQ_ChebExtrema(Grid_T *const grid);
static int GQ_Lobatto(Grid_T *const grid);
static int GQ_Legendre(Grid_T *const grid);
static int fdV_spectral(Grid_T *const grid);
static int fdS_spectral(Grid_T *const grid);

