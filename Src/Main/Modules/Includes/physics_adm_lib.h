#ifndef physics_adm_LIB_H
#define physics_adm_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PHYSICS_T;
struct PATCH_T;

int adm_main(struct PHYSICS_T *const phys);
void adm_update_beta_U0(struct PATCH_T *const patch);
void adm_update_beta_U1(struct PATCH_T *const patch);
void adm_update_beta_U2(struct PATCH_T *const patch);


#endif


