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
void bbn_update_derivative_rho0(Patch_T *const patch);
void bbn_update_derivative_u0(Patch_T *const patch);
void bbn_update_derivative_HS(Patch_T *const patch);
void bbn_update_derivative_K(Patch_T *const patch);
void bbn_update_derivative_psi(Patch_T *const patch);
void bbn_update_derivative_eta(Patch_T *const patch);
void bbn_update_derivative_B0_U0(Patch_T *const patch);
void bbn_update_derivative_B0_U1(Patch_T *const patch);
void bbn_update_derivative_B0_U2(Patch_T *const patch);
void bbn_update_derivative_Beta_U0(Patch_T *const patch);
void bbn_update_derivative_Beta_U1(Patch_T *const patch);
void bbn_update_derivative_Beta_U2(Patch_T *const patch);
void bbn_update_Beta_U0(Patch_T *const patch);
void bbn_update_Beta_U1(Patch_T *const patch);
void bbn_update_Beta_U2(Patch_T *const patch);

