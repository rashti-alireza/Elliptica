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
#include "physics_transformation_lib.h"

/* root finder struc for force balance equation */  
struct Force_Balance_RootFinder_S
{
  Patch_T *patch;
  //Root_Finder_T *root_finder;
  double dLnGamma;
  double y_CM;
  double x_CM;
  double Vr;
  double D;
  double Omega_BHNS;
  const double *X;
  unsigned find_y_CM: 1;
  unsigned find_x_CM: 1;
  unsigned find_Omega: 1;
  int dir;/* direction of derivative */
};

void bbn_study_initial_data(Grid_T *const grid);
void bbn_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder);
void bbn_update_psi10A_UiUj(Patch_T *const patch);
void bbn_populate_free_data(Grid_T *const grid);
void bbn_free_data_gammas(Grid_T *const grid);
void bbn_free_conformal_metric_derivatives(Patch_T *const patch);
void bbn_preparing_conformal_metric_derivatives(Patch_T *const patch);
void bbn_free_data_Gamma(Grid_T *const grid);
void bbn_free_data_Ricci(Grid_T *const grid);
void bbn_free_data_KS_trKij(Patch_T *const patch);
void bbn_update_derivative_phi(Patch_T *const patch);
void bbn_update_derivative_enthalpy(Patch_T *const patch);
void bbn_update_derivative_phi(Patch_T *const patch);
void bbn_update_rho0(Patch_T *const patch);
void bbn_update_derivative_rho0(Patch_T *const patch);
void bbn_update_derivative_u0(Patch_T *const patch);
void bbn_update_derivative_HS(Patch_T *const patch);
void bbn_update_derivative_K(Patch_T *const patch);
void bbn_update_derivative_psi(Patch_T *const patch);
void bbn_update_derivative_eta(Patch_T *const patch);
void bbn_update_derivative_B0_U0(Patch_T *const patch);
void bbn_update_derivative_B0_U1(Patch_T *const patch);
void bbn_update_derivative_B0_U2(Patch_T *const patch);
void bbn_update_derivative_B1_U0(Patch_T *const patch);
void bbn_update_derivative_B1_U1(Patch_T *const patch);
void bbn_update_derivative_B1_U2(Patch_T *const patch);
void bbn_update_derivative_Beta_U0(Patch_T *const patch);
void bbn_update_derivative_Beta_U1(Patch_T *const patch);
void bbn_update_derivative_Beta_U2(Patch_T *const patch);
void bbn_update_Beta_U0(Patch_T *const patch);
void bbn_update_Beta_U1(Patch_T *const patch);
void bbn_update_Beta_U2(Patch_T *const patch);
void bbn_update_B1_U012(Patch_T *const patch);
double bbn_NS_baryonic_mass(Grid_T *const grid,const double Euler_C);
double force_balance_root_finder_eq(void *params,const double *const x);
double dLnGamma_in_force_balance_eq(Patch_T *const patch,const double *const NS_centerX,const int dir);
void bbn_calculate_constraints(Grid_T *const grid);
void bbn_update_enthalpy_and_denthalpy(Grid_T *const grid);
void bbn_update_stress_energy_tensor(Grid_T *const grid,const int flag);
double bbn_KerrSchild_H(const double M_BH,const double rbar,const double a,const double z);
double bbn_KerrShcild_r(const double x,const double y,const double z,const double a);
void bbn_extrapolate_metric_fields_insideBH(Grid_T *const grid);
double bbn_BH_Kommar_mass(Grid_T *const grid);
