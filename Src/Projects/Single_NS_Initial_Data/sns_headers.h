#include "core_lib.h"
#include "TOV_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"
#include "memory_managing_lib.h"
#include "utilities_lib.h"
#include "maths_approximation_lib.h"
#include "maths_calculus_lib.h"
#include "physics_EoS_lib.h"
#include "physics_StressEnergyTensor_lib.h"

/* root finder struct for NS surface eq */
struct NS_surface_RootFinder_S
{
  Patch_T *patch;
  void *root_finder;
  double x0[3];/* (x,y,z) at the surface */
  double *N;/* the direction of increasing or decreasing of x = x0+N*d */
  double Euler_C;/* Euler equation const. */
  double scale;/* to avoid long step in root finder */
  double maxR;/* max R allowed for NS surrounding */
};


void sns_study_initial_data(Grid_T *const grid);
void sns_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder);
void sns_update_psi10A_UiUj(Patch_T *const patch);
void sns_populate_free_data(Grid_T *const grid);
void sns_free_data_gammas(Grid_T *const grid);
void sns_free_conformal_metric_derivatives(Patch_T *const patch);
void sns_preparing_conformal_metric_derivatives(Patch_T *const patch);
void sns_free_data_Gamma(Grid_T *const grid);
void sns_free_data_Ricci(Grid_T *const grid);
void sns_free_data_KS_trKij(Patch_T *const patch);
void sns_update_derivative_phi(Patch_T *const patch);
void sns_update_derivative_enthalpy(Patch_T *const patch);
void sns_update_derivative_phi(Patch_T *const patch);
void sns_update_rho0(Patch_T *const patch);
void sns_update_derivative_rho0(Patch_T *const patch);
void sns_update_derivative_u0(Patch_T *const patch);
void sns_update_derivative_HS(Patch_T *const patch);
void sns_update_derivative_K(Patch_T *const patch);
void sns_update_derivative_psi(Patch_T *const patch);
void sns_update_derivative_eta(Patch_T *const patch);
void sns_update_derivative_B0_U0(Patch_T *const patch);
void sns_update_derivative_B0_U1(Patch_T *const patch);
void sns_update_derivative_B0_U2(Patch_T *const patch);
void sns_update_derivative_Beta_U0(Patch_T *const patch);
void sns_update_derivative_Beta_U1(Patch_T *const patch);
void sns_update_derivative_Beta_U2(Patch_T *const patch);
void sns_update_Beta_U0(Patch_T *const patch);
void sns_update_Beta_U1(Patch_T *const patch);
void sns_update_Beta_U2(Patch_T *const patch);
double sns_NS_surface_enthalpy_eq(void *params,const double *const x);
void sns_update_enthalpy_and_denthalpy(Grid_T *const grid);
void sns_update_matter_fields(Grid_T *const grid);
