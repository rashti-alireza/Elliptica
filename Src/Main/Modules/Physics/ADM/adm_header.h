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

int adm_main(Physics_T *const phys);
void adm_add_3plus1_fields(Grid_T *const grid);
void adm_update_adm_Kij(Patch_T *const patch);
void adm_update_adm_KIJ(Patch_T *const patch);
void adm_update_adm_g(Patch_T *const patch);
void adm_ham_and_mom_from_identities(Patch_T *const patch,
              const char *const Ham,const char *const Mom);
              
void adm_compute_constraints(Physics_T *const phys,
                                  const char *const region,
                                  const char *const method,
                                  const char *const ham,
                                  const char *const mom);
#endif

