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

/* set default parameters. */
static int set_equation_params(Physics_T *const phys)
{
  FUNC_TIC
  
  /* equation name and region to be solved
  // format: par_name = equation_name,region.
  // for instance, Eq_psi=XCTS_curve_excision_KerrSchild_DDM, .*
  // which shows Eq_psi is "XCTS_curve_excision_KerrSchild_DDM" and
  // it is supposed to be solved everywhere. */
  
  /* params:
  // ==============
  // phi,phi1,ph2: irrotaional piece of NS fluid
  // psi         : conformal factor
  // alphaPsi:   : lapse * psi
  // beta_?      : shifts
  //
  // options:
  // XCTS        : extended conformal thin sandwich method
  // curve       : metric is not flat
  // excision    : BH is excised so we have B.C on AH
  // DDM         : DDM_Schur_Complement method for solve.
  // Wolf        : Wolfgang Tichy's formalism
  // KerrSchild  : Using KerrSchild values for BC
  // .*          : everywhere
  // NS?          : only for patches covering NS? */
  
  Pset_default(P_"phi","XCTS_curve_Wolf_DDM,NS");
  Pset_default(P_"phi1","XCTS_curve_Wolf_DDM,NS1");
  Pset_default(P_"phi2","XCTS_curve_Wolf_DDM,NS2");
  Pset_default(P_"psi","XCTS_curve_excision_KerrSchild_DDM,.*");
  Pset_default(P_"alphaPsi","XCTS_curve_excision_KerrSchild_DDM,.*");
  Pset_default(P_"beta_x","XCTS_curve_excision_KerrSchild_DDM,.*");
  Pset_default(P_"beta_y","XCTS_curve_excision_KerrSchild_DDM,.*");
  Pset_default(P_"beta_z","XCTS_curve_excision_KerrSchild_DDM,.*");
  
  /* please add all of the possible equations below as shown.
  // there is no problem to have same equation with multiple names.
  // also note the convention used for name of the equations. */
  
  /* add equations data base */
  field_eq          = init_eq();
  bc_eq             = init_eq();
  jacobian_field_eq = init_eq();
  jacobian_bc_eq    = init_eq();
  
  /* XCTS_curve_Wolf_ddm_phi equations */
  add_eq(&field_eq,eq_XCTS_curve_Wolf_ddm_eq_phi,
         "eq_XCTS_curve_Wolf_ddm_phi");
  add_eq(&bc_eq ,eq_XCTS_curve_Wolf_ddm_bc_phi,
         "bc_XCTS_curve_Wolf_ddm_phi");
  add_eq(&jacobian_field_eq,eq_XCTS_curve_Wolf_ddm_jacobian_eq_phi,
         "jacobian_eq_XCTS_curve_Wolf_ddm_phi");
  add_eq(&jacobian_bc_eq,eq_XCTS_curve_Wolf_ddm_jacobian_bc_phi,
         "jacobian_bc_XCTS_curve_Wolf_ddm_phi");
  
  /* XCTS_curve_Wolf_ddm_phi1 equations */
  add_eq(&field_eq,eq_XCTS_curve_Wolf_ddm_eq_phi,
         "eq_XCTS_curve_Wolf_ddm_phi1");
  add_eq(&bc_eq,eq_XCTS_curve_Wolf_ddm_bc_phi,
         "bc_XCTS_curve_Wolf_ddm_phi1");
  add_eq(&jacobian_field_eq,eq_XCTS_curve_Wolf_ddm_jacobian_eq_phi,
         "jacobian_eq_XCTS_curve_Wolf_ddm_phi1");
  add_eq(&jacobian_bc_eq,eq_XCTS_curve_Wolf_ddm_jacobian_bc_phi,
         "jacobian_bc_XCTS_curve_Wolf_ddm_phi1");
  
  /* XCTS_curve_Wolf_ddm_phi2 equations */
  add_eq(&field_eq,eq_XCTS_curve_Wolf_ddm_eq_phi,
         "eq_XCTS_curve_Wolf_ddm_phi2");
  add_eq(&bc_eq,eq_XCTS_curve_Wolf_ddm_bc_phi,
         "bc_XCTS_curve_Wolf_ddm_phi2");
  add_eq(&jacobian_field_eq,eq_XCTS_curve_Wolf_ddm_jacobian_eq_phi,
         "jacobian_eq_XCTS_curve_Wolf_ddm_phi2");
  add_eq(&jacobian_bc_eq,eq_XCTS_curve_Wolf_ddm_jacobian_bc_phi,
         "jacobian_bc_XCTS_curve_Wolf_ddm_phi2");
  
  /* XCTS_curve_excision_KerrSchild_ddm_psi equations */
  add_eq(&field_eq,eq_XCTS_curve_exc_KS_ddm_eq_psi,
         "eq_XCTS_curve_excision_KerrSchild_ddm_psi");
  add_eq(&bc_eq,eq_XCTS_curve_exc_KS_ddm_bc_psi,
         "bc_XCTS_curve_excision_KerrSchild_ddm_psi");
  add_eq(&jacobian_field_eq,eq_XCTS_curve_exc_KS_ddm_jacobian_eq_psi,
         "jacobian_eq_XCTS_curve_excision_KerrSchild_ddm_psi");
  add_eq(&jacobian_bc_eq,eq_XCTS_curve_exc_KS_ddm_jacobian_bc_psi,
         "jacobian_bc_XCTS_curve_excision_KerrSchild_ddm_psi");

  /* XCTS_curve_excision_KerrSchild_ddm_alphaPsi equations */
  add_eq(&field_eq,eq_XCTS_curve_exc_KS_ddm_eq_alphaPsi,
         "eq_XCTS_curve_excision_KerrSchild_ddm_alphaPsi");
  add_eq(&bc_eq,eq_XCTS_curve_exc_KS_ddm_bc_alphaPsi,
        "bc_XCTS_curve_excision_KerrSchild_ddm_alphaPsi");
  add_eq(&jacobian_field_eq,eq_XCTS_curve_exc_KS_ddm_jacobian_eq_alphaPsi,
        "jacobian_eq_XCTS_curve_excision_KerrSchild_ddm_alphaPsi");
  add_eq(&jacobian_bc_eq,eq_XCTS_curve_exc_KS_ddm_jacobian_bc_alphaPsi,
        "jacobian_bc_XCTS_curve_excision_KerrSchild_ddm_alphaPsi");
  
  /* XCTS_curve_excision_KerrSchild_ddm_beta_x equations */
  add_eq(&field_eq,eq_XCTS_curve_exc_KS_ddm_eq_beta_U0,
         "eq_XCTS_curve_excision_KerrSchild_ddm_beta_x");
  add_eq(&bc_eq ,eq_XCTS_curve_exc_KS_ddm_bc_beta_U0,
         "bc_XCTS_curve_excision_KerrSchild_ddm_beta_x");
  add_eq(&jacobian_field_eq,eq_XCTS_curve_exc_KS_ddm_jacobian_eq_beta_U0,
         "jacobian_eq_XCTS_curve_excision_KerrSchild_ddm_beta_x");
  add_eq(&jacobian_bc_eq,eq_XCTS_curve_exc_KS_ddm_jacobian_bc_beta_U0,
         "jacobian_bc_XCTS_curve_excision_KerrSchild_ddm_beta_x");

  /* XCTS_curve_excision_KerrSchild_ddm_beta_y equations */
  add_eq(&field_eq,eq_XCTS_curve_exc_KS_ddm_eq_beta_U1,
         "eq_XCTS_curve_excision_KerrSchild_ddm_beta_y");
  add_eq(&bc_eq,eq_XCTS_curve_exc_KS_ddm_bc_beta_U1,
         "bc_XCTS_curve_excision_KerrSchild_ddm_beta_y");
  add_eq(&jacobian_field_eq,eq_XCTS_curve_exc_KS_ddm_jacobian_eq_beta_U1,
         "jacobian_eq_XCTS_curve_excision_KerrSchild_ddm_beta_y");
  add_eq(&jacobian_bc_eq,eq_XCTS_curve_exc_KS_ddm_jacobian_bc_beta_U1,
         "jacobian_bc_XCTS_curve_excision_KerrSchild_ddm_beta_y");

  /* XCTS_curve_excision_KerrSchild_ddm_beta_z equations */
  add_eq(&field_eq,eq_XCTS_curve_exc_KS_ddm_eq_beta_U2,
         "eq_XCTS_curve_excision_KerrSchild_ddm_beta_z");
  add_eq(&bc_eq,eq_XCTS_curve_exc_KS_ddm_bc_beta_U2,
         "bc_XCTS_curve_excision_KerrSchild_ddm_beta_z");
  add_eq(&jacobian_field_eq,eq_XCTS_curve_exc_KS_ddm_jacobian_eq_beta_U2,
         "jacobian_eq_XCTS_curve_excision_KerrSchild_ddm_beta_z");
  add_eq(&jacobian_bc_eq,eq_XCTS_curve_exc_KS_ddm_jacobian_bc_beta_U2,
         "jacobian_bc_XCTS_curve_excision_KerrSchild_ddm_beta_z");
  

  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* adding fields. */
static int add_equation_fields(Physics_T *const phys)
{
  FUNC_TIC
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* solve fields. */
static int solve_equation(Physics_T *const phys)
{
  FUNC_TIC
  
  /* populating solution managment */
  initialize_solving_man(grid,field_eq,bc_eq,jacobian_field_eq,jacobian_bc_eq);
  
  /* solving equation(s) */
  Solve_Equations_T *SolveEqs = init_solve_equations(grid);
  SolveEqs->solving_order = Pgets("Solving_Order");
  SolveEqs->FieldUpdate   = bbn_SolveEqs_FieldUpdate;
  SolveEqs->SourceUpdate  = bbn_SolveEqs_SourceUpdate;
  SolveEqs->StopCriteria  = bbn_stop_criteria;
  SolveEqs->umfpack_refine= PgetdEZ("Solving_UMFPACK_refinement_step");
  SolveEqs->umfpack_size  = PgetiEZ("Solving_UMFPACK_size");
  Grid_T *phi_grid = bbn_phi_grid(grid);/* phi needed to be solved only in NS */
  add_special_grid_solve_equations(phi_grid,"phi",SolveEqs);
  
  /* saving the field being solved for relaxation scheme purposes */
  save_fields(grid);

  const Uint max_iter = (Uint)Pgeti("Solving_Max_Number_of_Iteration");
  const int max_newton_step = Pgeti("Solving_Max_Number_of_Newton_Step");
  Uint iter = 0;
  
  while (iter < max_iter)
  {
    /* some prints */
    pr_line_custom('=');
    printf("{ Solving XCTS equations at a fixed resolution ...\n");
    printf("  |---> Iteration      = %d / %d\n",iter+1,max_iter);
    printf("  |---> Equation(s)    = %s\n",SolveEqs->solving_order);
    printf("  |---> Newton step(s) = %d\n",max_newton_step);
    fflush(stdout);
    
    /* solve equations */
    solve_eqs(SolveEqs);
    
    /* study the solution */
    if (Pcmps("Elliptic_Convergence_Test","yes"))
    {
      calculate_equation_residual(SolveEqs);
      bbn_study_initial_data(grid);
      bbn_write_checkpoint(grid);
    }
    
    /* some prints */
    
    printf("} Solving XCTS equations at a fixed resolution ==> Done.\n ");
    pr_clock();
    
    ++iter;
  }
  pr_line_custom('=');

  /* updating the fields using relaxed scheme */
  update_fields_relaxed_scheme(grid);

  /* calculate the field residual for diagnostic purposes */
  calculate_equation_residual(SolveEqs);

  /* free SolveEqs and phi grid */
  free_solve_equations(SolveEqs);
  bbn_free_phi_grid(phi_grid);

  /* free data base of equations */
  free_db_eqs(field_eq);
  free_db_eqs(bc_eq);
  free_db_eqs(jacobian_field_eq);
  free_db_eqs(jacobian_bc_eq);
  
  field_eq = 0;
  bc_eq    = 0;
  jacobian_field_eq = 0;
  jacobian_bc_eq    = 0;
  
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}


