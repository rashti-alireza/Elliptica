#ifndef adm_LIB_H
#define adm_LIB_H

#include "core_lib.h"
#include "fields_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "maths_diff_geom_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_lib.h"

/* parameter prefix */
#define P_ "adm_"

/* function typedefs */
typedef void fFunc_adm_update_AConfIJ_T(Patch_T *const patch);
typedef void fFunc_adm_update_adm_Kij_T(Patch_T *const patch);
typedef void fFunc_adm_update_adm_KIJ_T(Patch_T *const patch);


int adm_main(Physics_T *const phys);
void adm_add_3plus1_fields(Grid_T *const grid);
void adm_update_adm_Kij_useAIJ(Patch_T *const patch);
void adm_update_adm_KIJ_useAIJ(Patch_T *const patch);
void adm_update_adm_KIJ_useKij(Patch_T *const patch);
void adm_update_adm_g(Patch_T *const patch);
void adm_update_AConfIJ_XCTS_MConfIJ0(Patch_T *const patch);
void adm_update_adm_Kij(Physics_T *const phys,const char *const region);
void adm_update_adm_KIJ(Physics_T *const phys,const char *const region);
void adm_ham_and_mom_from_identities(Patch_T *const patch,
              const char *const Ham,const char *const Mom);
void adm_ham_and_mom_from_scratch(Patch_T *const patch,
              const char *const Ham,const char *const Mom);
              
void adm_compute_constraints(Physics_T *const phys,
                                  const char *const region,
                                  const char *const method,
                                  const char *const ham,
                                  const char *const mom);

void adm_update_AConfIJ(Physics_T *const phys,const char *const region);

#endif

