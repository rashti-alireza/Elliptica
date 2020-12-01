#include "adm_header.h"


fFunc_adm_update_AConfIJ_T (*adm_update_AConfIJ_patch);
fFunc_adm_update_adm_Kij_T (*adm_update_adm_Kij_patch);
fFunc_adm_update_adm_KIJ_T (*adm_update_adm_KIJ_patch);

int adm_main(Physics_T *const phys);
static int add_adm_params(Physics_T *const phys);
static int add_adm_fields(Physics_T *const phys);
static int compute_ham_and_mom_constraints(Physics_T *const phys);
static int compute_AConfIJ(Physics_T *const phys);
static int compute_adm_Kij(Physics_T *const phys);
static int compute_adm_KIJ(Physics_T *const phys);
static int compute_adm_gij(Physics_T *const phys);


