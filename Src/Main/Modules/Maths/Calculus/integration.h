#include "core_lib.h"
#include "maths_calculus_lib.h"
#include "maths_general_lib.h"
#include "maths_analytic_lib.h"
#include "coordinates_lib.h"
#include "utilities_lib.h"

void free_integration(Integration_T *I);
void plan_integration(Integration_T *const I);
double execute_integration(Integration_T *const I);
Integration_T *init_integration(void);
static double Composite_Simpson_1D(Integration_T *const I);
static double GaussQuadrature_ChebyshevExtrema(Integration_T *const I);
static double GaussQuadrature_Lobatto(Integration_T *const I);
static double GaussQuadrature_Legendre(Integration_T *const I);
static double J_xyzN0N1N2(Patch_T *const patch,const unsigned ijk);
static double Int_ChebTn(const unsigned n,const unsigned N);
static double f_xyz_dV_Cheb_Ext_Spec(Integration_T *const I);






