#include "bh_header.h" 
#include "physics_freedata_lib.h"


void bh_tune_BH_radius_irreducible_mass_perfect_s2(Physics_T *const phys);
void bh_find_bh_surface_perfect_s2(Physics_T *const phys);
void bh_start_off_KerrSchild_perfect_s2(Physics_T *const phys);
void bh_start_off_IsoSchild_perfect_s2(Physics_T *const phys);
void bh_start_off_PGSchild_perfect_s2(Physics_T *const phys);
void bh_start_off_KerrSchild_general_s2(Physics_T *const phys);
void bh_find_bh_surface_KerrSchild_s2(Physics_T *const phys);
void bh_update_sConf_dsConf(Physics_T *const phys);
void bh_update_inner_BC(Physics_T *const phys);
void bh_start_off_CloseKerrSchild_perfect_s2(Physics_T *const phys);


static void 
set_beta_inner_bc_alpha_omegaXr
 (
 Physics_T *const phys,
 const char *const region,
 const char *const ib_Beta
 );


void bh_tune_BH_chi_simple(Physics_T *const phys);



