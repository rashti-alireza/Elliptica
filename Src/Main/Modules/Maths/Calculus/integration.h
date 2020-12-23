#include "core_lib.h"
#include "maths_calculus_lib.h"
#include "maths_general_lib.h"
#include "maths_special_functions_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "fields_lib.h"

void free_integration(Integration_T *I);
void plan_integration(Integration_T *const I);
double execute_integration(Integration_T *const I);
Integration_T *init_integration(void);
double det_h_xyzN0N1_Cheb_Ext(Patch_T *const patch,const Integration_T *const I,const Uint ijk);
double det_h_xyzN0N2_Cheb_Ext(Patch_T *const patch,const Integration_T *const I,const Uint ijk);
double det_h_xyzN1N2_Cheb_Ext(Patch_T *const patch,const Integration_T *const I,const Uint ijk);
static double Composite_Simpson_1D(Integration_T *const I);
static double GaussQuadrature_ChebyshevExtrema(Integration_T *const I);
static double GaussQuadrature_Lobatto(Integration_T *const I);
static double GaussQuadrature_Legendre(Integration_T *const I);
static double J_xyzN0N1N2(Patch_T *const patch,const Uint ijk);
static double f_xyz_whole_dV_Cheb_Ext_Spec(Integration_T *const I);
static double f_xyz_interval_dV_Cheb_Ext_Spec(Integration_T *const I);
static double f_xyz_dS_Cheb_Ext_Spec(Integration_T *const I);
static double Int_ChebTn_OPTM(const Uint n,const Uint N);
static double integral_conventional_ChebTn(const Uint N,
                                           const Uint n,
                                           const double th_i,
                                           const double th_f);
