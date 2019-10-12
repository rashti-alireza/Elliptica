#include "sns_headers.h"

void sns_allocate_fields(Grid_T *const grid);
void sns_partial_derivatives_fields(Grid_T *const grid);
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
void sns_update_enthalpy_and_denthalpy(Grid_T *const grid);
void sns_update_matter_fields(Grid_T *const grid);
