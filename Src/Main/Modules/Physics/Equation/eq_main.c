/*
// Alireza Rashti
// December 2020
*/

/* equations manager, to set, to solve etc. */

#include "eq_main.h"

/* main function to issue commands */
int eq_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case EQ_SET_PARAMS:
      ret = set_equation_params(phys);
    break;
    
    case EQ_ADD_FIELDS:
      ret = add_equation_fields(phys);
    break;
    
    case EQ_SOLVE:
      ret = solve_equation(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* set default parameters. 
// the notation might be sophisticated. sorry! */
static int set_equation_params(Physics_T *const phys)
{
  FUNC_TIC
  
  /* prefix of system paramter (no leading '_') */
  Pset_default(P_"system_prefix",phys->ssys);
  
  /* equation type
  //
  // options:
  // ========
  // elliptic: elliptic PDE. */
  Pset_default(P_"type","elliptic");
  
  /* equation test
  //
  // options:
  // ========
  // elliptic_test: test and analyze elliptic eqs after each iteration. */
  Pset_default(P_"elliptic_test","no");
  
  /* field update method:
  //
  // options:
  // ========
  // relaxed_scheme: update fields using weights => new = W1*solved+W2*old */
  Pset_default(P_"update_method","relaxed_scheme");
  
  /* update_weight (soft paramter):
  // one should add required update weight as the following example:
  // Eq_update_weight_phi1 = 0.2. */
  
  /* equation name and region to be solved
  // format: par_name = equation_name,region.
  // for instance, Eq_psi=XCTS_curve_excision_Type1_DDM, .*
  // which shows Eq_psi is "XCTS_curve_excision_Type1_DDM" and
  // it is supposed to be solved everywhere. */
  
  /* params:
  // ==============
  // phi,phi1,ph2: irrotaional piece of NS fluid
  // psi         : conformal factor
  // alphaPsi:   : lapse * psi
  // beta_?      : shifts
  //
  // options:(explained)
  // XCTS         : extended conformal thin sandwich method
  // curve        : metric is not flat
  // excision     : BH is excised so we have B.C on AH
  // DDM          : DDM_Schur_Complement method for solve.
  // Type1(for phi eq)    : arXiv:1910.09690 [gr-qc]+arXiv:0804.3787v3[gr-qc], 
  //                        modify eq for rho0 + exact helical for DiBi.
  // Type2(for phi eq)    : arXiv:1910.09690 [gr-qc], modify eq for rho0.
  // Type1(for metric eqs): arXiv:1607.07962 [gr-qc] on AH.
  // Type2(for alphaPsi)  : set dalpha/dr = 0 on AH.
  // .*           : everywhere
  // NS?          : only for patches covering NS? */
  
  Pset_default(P_"phi","XCTS_curve_Type1_DDM,NS");
  Pset_default(P_"phi1","XCTS_curve_Type1_DDM,NS1");
  Pset_default(P_"phi2","XCTS_curve_Type1_DDM,NS2");
  Pset_default(P_"psi","XCTS_curve_excision_Type1_DDM,.*");
  Pset_default(P_"alphaPsi","XCTS_curve_excision_Type1_DDM,.*");
  Pset_default(P_"B0_U0","XCTS_curve_excision_Type1_DDM,.*");
  Pset_default(P_"B0_U1","XCTS_curve_excision_Type1_DDM,.*");
  Pset_default(P_"B0_U2","XCTS_curve_excision_Type1_DDM,.*");
  
  /* external functions */
  eq_field_update  = 0;
  eq_source_update = 0;
  eq_stop_criteria = 0;
  eq_analyze_solution = 0;
  
  /* please add all of the possible equations below as shown.
  // there is no problem to have same equation with multiple names.
  // also note the convention used for name of the equations. */
  
  /* free data base of equations if exists */
  free_db_eqs(eq_global_field_eq);
  free_db_eqs(eq_global_bc_eq);
  free_db_eqs(eq_global_jacobian_field_eq);
  free_db_eqs(eq_global_jacobian_bc_eq);
  
  /* add equations data base */
  eq_global_field_eq          = init_eq();
  eq_global_bc_eq             = init_eq();
  eq_global_jacobian_field_eq = init_eq();
  eq_global_jacobian_bc_eq    = init_eq();
  
  /* XCTS_curve_T1_ddm_phi equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_T1_ddm_eq_phi,
         "eq_XCTS_curve_Type1_ddm_phi");
  add_eq(&eq_global_bc_eq ,eq_XCTS_curve_T1_ddm_bc_phi,
         "bc_XCTS_curve_Type1_ddm_phi");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_T1_ddm_jacobian_eq_phi,
         "jacobian_eq_XCTS_curve_Type1_ddm_phi");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_T1_ddm_jacobian_bc_phi,
         "jacobian_bc_XCTS_curve_Type1_ddm_phi");
  
  /* XCTS_curve_T1_ddm_phi1 equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_T1_ddm_eq_phi,
         "eq_XCTS_curve_Type1_ddm_phi1");
  add_eq(&eq_global_bc_eq,eq_XCTS_curve_T1_ddm_bc_phi,
         "bc_XCTS_curve_Type1_ddm_phi1");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_T1_ddm_jacobian_eq_phi,
         "jacobian_eq_XCTS_curve_Type1_ddm_phi1");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_T1_ddm_jacobian_bc_phi,
         "jacobian_bc_XCTS_curve_Type1_ddm_phi1");
  
  /* XCTS_curve_T1_ddm_phi2 equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_T1_ddm_eq_phi,
         "eq_XCTS_curve_Type1_ddm_phi2");
  add_eq(&eq_global_bc_eq,eq_XCTS_curve_T1_ddm_bc_phi,
         "bc_XCTS_curve_Type1_ddm_phi2");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_T1_ddm_jacobian_eq_phi,
         "jacobian_eq_XCTS_curve_Type1_ddm_phi2");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_T1_ddm_jacobian_bc_phi,
         "jacobian_bc_XCTS_curve_Type1_ddm_phi2");
  
  /* XCTS_curve_T2_ddm_phi equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_T2_ddm_eq_phi,
         "eq_XCTS_curve_Type2_ddm_phi");
  add_eq(&eq_global_bc_eq ,eq_XCTS_curve_T1_ddm_bc_phi,
         "bc_XCTS_curve_Type2_ddm_phi");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_T1_ddm_jacobian_eq_phi,
         "jacobian_eq_XCTS_curve_Type2_ddm_phi");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_T1_ddm_jacobian_bc_phi,
         "jacobian_bc_XCTS_curve_Type2_ddm_phi");
         
  /* XCTS_curve_T2_ddm_phi1 equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_T2_ddm_eq_phi,
         "eq_XCTS_curve_Type2_ddm_phi1");
  add_eq(&eq_global_bc_eq ,eq_XCTS_curve_T1_ddm_bc_phi,
         "bc_XCTS_curve_Type2_ddm_phi1");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_T1_ddm_jacobian_eq_phi,
         "jacobian_eq_XCTS_curve_Type2_ddm_phi1");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_T1_ddm_jacobian_bc_phi,
         "jacobian_bc_XCTS_curve_Type2_ddm_phi1");
  
  /* XCTS_curve_T2_ddm_phi2 equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_T2_ddm_eq_phi,
         "eq_XCTS_curve_Type2_ddm_phi2");
  add_eq(&eq_global_bc_eq ,eq_XCTS_curve_T1_ddm_bc_phi,
         "bc_XCTS_curve_Type2_ddm_phi2");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_T1_ddm_jacobian_eq_phi,
         "jacobian_eq_XCTS_curve_Type2_ddm_phi2");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_T1_ddm_jacobian_bc_phi,
         "jacobian_bc_XCTS_curve_Type2_ddm_phi2");
  
  /* XCTS_curve_excision_Type1_ddm_psi equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_exc_T1_ddm_eq_psi,
         "eq_XCTS_curve_excision_Type1_ddm_psi");
  add_eq(&eq_global_bc_eq,eq_XCTS_curve_exc_T1_ddm_bc_psi,
         "bc_XCTS_curve_excision_Type1_ddm_psi");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_eq_psi,
         "jacobian_eq_XCTS_curve_excision_Type1_ddm_psi");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_bc_psi,
         "jacobian_bc_XCTS_curve_excision_Type1_ddm_psi");

  /* XCTS_curve_excision_Type1_ddm_alphaPsi equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_exc_T1_ddm_eq_alphaPsi,
         "eq_XCTS_curve_excision_Type1_ddm_alphaPsi");
  add_eq(&eq_global_bc_eq,eq_XCTS_curve_exc_T1_ddm_bc_alphaPsi,
        "bc_XCTS_curve_excision_Type1_ddm_alphaPsi");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_eq_alphaPsi,
        "jacobian_eq_XCTS_curve_excision_Type1_ddm_alphaPsi");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_bc_alphaPsi,
        "jacobian_bc_XCTS_curve_excision_Type1_ddm_alphaPsi");
  
  /* XCTS_curve_excision_Type2_ddm_alphaPsi equations(only AH bc changed) */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_exc_T1_ddm_eq_alphaPsi,
         "eq_XCTS_curve_excision_Type2_ddm_alphaPsi");
  add_eq(&eq_global_bc_eq,eq_XCTS_curve_exc_T2_ddm_bc_alphaPsi,
        "bc_XCTS_curve_excision_Type2_ddm_alphaPsi");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_eq_alphaPsi,
        "jacobian_eq_XCTS_curve_excision_Type2_ddm_alphaPsi");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_exc_T2_ddm_jacobian_bc_alphaPsi,
        "jacobian_bc_XCTS_curve_excision_Type2_ddm_alphaPsi");
  
  /* XCTS_curve_excision_Type1_ddm_B0_U0 equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_exc_T1_ddm_eq_B0_U0,
         "eq_XCTS_curve_excision_Type1_ddm_B0_U0");
  add_eq(&eq_global_bc_eq ,eq_XCTS_curve_exc_T1_ddm_bc_B0_U0,
         "bc_XCTS_curve_excision_Type1_ddm_B0_U0");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_eq_B0_U0,
         "jacobian_eq_XCTS_curve_excision_Type1_ddm_B0_U0");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_bc_B0_U0,
         "jacobian_bc_XCTS_curve_excision_Type1_ddm_B0_U0");

  /* XCTS_curve_excision_Type1_ddm_B0_U1 equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_exc_T1_ddm_eq_B0_U1,
         "eq_XCTS_curve_excision_Type1_ddm_B0_U1");
  add_eq(&eq_global_bc_eq,eq_XCTS_curve_exc_T1_ddm_bc_B0_U1,
         "bc_XCTS_curve_excision_Type1_ddm_B0_U1");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_eq_B0_U1,
         "jacobian_eq_XCTS_curve_excision_Type1_ddm_B0_U1");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_bc_B0_U1,
         "jacobian_bc_XCTS_curve_excision_Type1_ddm_B0_U1");

  /* XCTS_curve_excision_Type1_ddm_B0_U2 equations */
  add_eq(&eq_global_field_eq,eq_XCTS_curve_exc_T1_ddm_eq_B0_U2,
         "eq_XCTS_curve_excision_Type1_ddm_B0_U2");
  add_eq(&eq_global_bc_eq,eq_XCTS_curve_exc_T1_ddm_bc_B0_U2,
         "bc_XCTS_curve_excision_Type1_ddm_B0_U2");
  add_eq(&eq_global_jacobian_field_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_eq_B0_U2,
         "jacobian_eq_XCTS_curve_excision_Type1_ddm_B0_U2");
  add_eq(&eq_global_jacobian_bc_eq,eq_XCTS_curve_exc_T1_ddm_jacobian_bc_B0_U2,
         "jacobian_bc_XCTS_curve_excision_Type1_ddm_B0_U2");
  

  FUNC_TOC
  return EXIT_SUCCESS;
}

/* adding fields. */
static int add_equation_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid = phys->grid;
  char **fields_name = read_separated_items_in_string(Pgets("solve_Order"),',');
  Uint f;
  
  f = 0;
  while (fields_name[f])
  {
    /* add backup fields */
    char field_backup[STR_LEN];
    sprintf(field_backup,P_Backup_"%s",fields_name[f]);
    FOR_ALL_p(grid->np)
    {
      Patch_T *patch = grid->patch[p];
      add_field(field_backup,0,patch,NO);
    }
    
    /* add residual fields */
    char field_res[STR_LEN];
    sprintf(field_res,"%s"EQ_Residual_Suffix,fields_name[f]);
    FOR_ALL_p(grid->np)
    {
      Patch_T *patch = grid->patch[p];
      add_field(field_res,0,patch,NO);
    }
    
    f++;
  }
  free_2d(fields_name);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* solve fields. */
static int solve_equation(Physics_T *const phys)
{
  FUNC_TIC
  
  if (Pcmps(P_"type","elliptic"))
    eq_solve_elliptic_equation(phys);
  else
    Error0(NO_OPTION);
    
  FUNC_TOC
  return EXIT_SUCCESS;
}

