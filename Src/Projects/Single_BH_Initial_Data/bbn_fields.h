#include "bbn_headers.h"

void bbn_allocate_fields(Grid_T *const grid);
void bbn_partial_derivatives_fields(Grid_T *const grid);
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
void bbn_update_derivative_Beta_U0(Patch_T *const patch);
void bbn_update_derivative_Beta_U1(Patch_T *const patch);
void bbn_update_derivative_Beta_U2(Patch_T *const patch);
void bbn_update_Beta_U0(Patch_T *const patch);
void bbn_update_Beta_U1(Patch_T *const patch);
void bbn_update_Beta_U2(Patch_T *const patch);
void bbn_update_enthalpy_and_denthalpy(Grid_T *const grid);
void bbn_update_stress_energy_tensor(Grid_T *const grid);
static void cleaning_enthalpy(Patch_T *const patch);

