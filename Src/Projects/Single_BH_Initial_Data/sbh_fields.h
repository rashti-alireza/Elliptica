#include "sbh_headers.h"

void sbh_allocate_fields(Grid_T *const grid);
void sbh_partial_derivatives_fields(Grid_T *const grid);
void sbh_update_derivative_HS(Patch_T *const patch);
void sbh_update_derivative_K(Patch_T *const patch);
void sbh_update_derivative_psi(Patch_T *const patch);
void sbh_update_derivative_eta(Patch_T *const patch);
void sbh_update_derivative_B0_U0(Patch_T *const patch);
void sbh_update_derivative_B0_U1(Patch_T *const patch);
void sbh_update_derivative_B0_U2(Patch_T *const patch);
void sbh_update_derivative_Beta_U0(Patch_T *const patch);
void sbh_update_derivative_Beta_U1(Patch_T *const patch);
void sbh_update_derivative_Beta_U2(Patch_T *const patch);
void sbh_update_Beta_U0(Patch_T *const patch);
void sbh_update_Beta_U1(Patch_T *const patch);
void sbh_update_Beta_U2(Patch_T *const patch);

