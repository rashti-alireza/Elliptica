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

/* parameter prefix
// PLEASE keep it capitalized. */
#define P_ "ADM_"

/* general struct for argument pass */
struct General_Arg_S
{
  /* for B1I */
  double omega;
  double CM[3];
  double Vr;
  double D;
};

/* function typedefs */
typedef void fFunc_adm_update_AConfIJ_T(Patch_T *const patch);
typedef void fFunc_adm_update_adm_Kij_T(Patch_T *const patch);
typedef void fFunc_adm_update_adm_KIJ_T(Patch_T *const patch);
typedef void fFunc_adm_update_B1I_T(Patch_T *const patch,void *params);


int adm_main(Physics_T *const phys);
void adm_update_adm_B1I(Physics_T *const phys,const char *const region);
void adm_update_B1I_inspiral(Patch_T *const patch,void *params);
void adm_update_B1I_zero(Patch_T *const patch,void *params);
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
void adm_update_adm_g_patch(Patch_T *const patch);
void adm_update_adm_gij(Physics_T *const phys,const char *const region);
void adm_update_beta(Physics_T *const phys,const char *const region);
void adm_doctest_AConfIJ(Physics_T *const phys);
void adm_update_AConfIJ_SCTT_MConfIJ0(Patch_T *const patch);


#endif

