#include "adm_header.h"

extern fFunc_adm_update_AConfIJ_T (*adm_update_AConfIJ_patch);
extern fFunc_adm_update_adm_Kij_T (*adm_update_adm_Kij_patch);
extern fFunc_adm_update_adm_KIJ_T (*adm_update_adm_KIJ_patch);

void adm_compute_constraints(Physics_T *const phys,
                                  const char *const region,
                                  const char *const method,
                                  const char *const ham,
                                  const char *const mom);
void adm_update_AConfIJ(Physics_T *const phys,const char *const region);
void adm_update_adm_Kij(Physics_T *const phys,const char *const region);
void adm_update_adm_KIJ(Physics_T *const phys,const char *const region);
void adm_update_adm_gij(Physics_T *const phys,const char *const region);


