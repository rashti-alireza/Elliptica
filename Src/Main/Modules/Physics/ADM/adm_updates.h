#include "adm_header.h"

/* set given field x to 0. */
#define set_field_to_zero(x) REALLOC_v_WRITE_v(x); UNUSED(x);

extern fFunc_adm_update_AConfIJ_T (*adm_update_AConfIJ_patch);
extern fFunc_adm_update_adm_Kij_T (*adm_update_adm_Kij_patch);
extern fFunc_adm_update_adm_KIJ_T (*adm_update_adm_KIJ_patch);
extern fFunc_adm_update_B1I_T (*adm_update_B1I_patch);


void adm_compute_constraints(Physics_T *const phys,
                                  const char *const region,
                                  const char *const method,
                                  const char *const ham,
                                  const char *const mom);
void adm_update_AConfIJ(Physics_T *const phys,const char *const region);
void adm_update_adm_Kij(Physics_T *const phys,const char *const region);
void adm_update_adm_KIJ(Physics_T *const phys,const char *const region);
void adm_update_adm_gij(Physics_T *const phys,const char *const region);
void adm_update_adm_B1I(Physics_T *const phys,const char *const region);
void adm_update_B1I_inspiral(Patch_T *const patch,void *params);
void adm_update_B1I_zero(Patch_T *const patch,void *params);
void adm_update_beta(Physics_T *const phys,const char *const region);
void adm_update_beta_U0(Patch_T *const patch);
void adm_update_beta_U1(Patch_T *const patch);
void adm_update_beta_U2(Patch_T *const patch);


