/*
// Alireza Rashti
// October 2019
*/

#include "sbh_solve_initial_data_eqs.h"

/* solving initial data equations for the given grid */
void sbh_solve_initial_data_eqs(Grid_T *const grid)
{
  pr_line_custom('='); 
  printf("{ Solving initial data equations for Single BH ...\n");
  
  sEquation_T **field_eq/* field equation */,
              **bc_eq/* B.C. for the field */,
              **jacobian_field_eq/* jacobian for field equation */,
              **jacobian_bc_eq/* jacobian for B.C. */;

  /* filling db equations of XCTS */
  sbh_XCTS_fill_db_eqs(&field_eq,&bc_eq,&jacobian_field_eq,&jacobian_bc_eq);
  
  /* populating solution managment */
  initialize_solving_man(grid,field_eq,bc_eq,jacobian_field_eq,jacobian_bc_eq);
  
  /* solving equation(s) */
  Solve_Equations_T *SolveEqs = init_solve_equations(grid);
  SolveEqs->solving_order = GetParameterS_E("Solving_Order");
  SolveEqs->FieldUpdate  = sbh_SolveEqs_FieldUpdate;
  SolveEqs->StopCriteria = sbh_stop_criteria;
  
  const int max_iter = GetParameterI_E("Solving_Max_Number_of_Iteration");
  int iter = 0;
  while (iter < max_iter)
  {
    /* some prints */
    pr_line_custom('=');
    printf("{ Iteration %d For Solving XCTS Equations at a Fixed Resolution ...\n",iter);
    printf("        |---> %s Equations ...\n",SolveEqs->solving_order);
    
    solve_eqs(SolveEqs);
    
    /* some prints */
    printf("} Iteration %d For Solving XCTS Equations at a Fixed Resolution ==> Done.\n",iter);
    pr_clock();
    
    ++iter;
  }
  pr_line_custom('=');
  
  /* free SolveEqs */
  free_solve_equations(SolveEqs);
  
  /* free data base of equations */
  free_db_eqs(field_eq);
  free_db_eqs(bc_eq);
  free_db_eqs(jacobian_field_eq);
  free_db_eqs(jacobian_bc_eq);
  
  printf("} Solving initial data equations for Single BH ==> Done.\n");
  pr_clock();
  pr_line_custom('='); 
}

/* stop criteria for solver, namely, if some conditions satisfied, 
// it stops. one can give specific criteria according to the given field's name.
// ->return value: 0 means stop, 1 means continue to solve */
int sbh_stop_criteria(Grid_T *const grid,const char *const name)
{
  int stop = 1;
  int stop_max = 1;
  int stop_res = 0;
  int stop_backtrack = 1;
  const double res_TOLERANCE = 1E-10;/* this is the tolerance that solver allowed to increase residual */
  const double res_bckt = GetParameterD_E("Solving_Allowed_Relative_Residual_Backtrack_Tolerance");
  const double res_d    = GetParameterD_E("Solving_Residual");/* desired residual */
  const int max_step    = GetParameterI_E("Solving_Max_Number_of_Newton_Step");
  const double res_fac  = GetParameterD_E("Solving_Residual_Factor");
  const unsigned npatch = grid->np;
  unsigned p;
  
  /* if no step should be taken */
  if (max_step  == 0)
  {
    printf("%s equation:\n"
           "---> Newton solver reached maximum step number so existing ...\n",name);
    fflush(stdout);
    return 0;
  }
    
  
  /* NOTE: due to the break command, the order of ifs are important */
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    double res      = patch->solving_man->Frms;/* current residual */
    double res_last;
    int solver_step = patch->solving_man->settings->solver_step;/* iteration number */
    
    /* if this is the very first step, don't check the following */
    if (solver_step  == 0)
      continue;
    
    /* if residual increased stop */
    res_last = patch->solving_man->settings->HFrms[solver_step-1];
    if (res > res_last*(1. + res_bckt) && GRT(res,res_TOLERANCE))
    {
      stop_backtrack = 0;
      break;
    }
    
    /* note: all patches have same solver_step */
    if (solver_step >= max_step)
    {
      stop_max = 0;
      break;
    }
  }
  
  if (!stop_backtrack)
  {
    /* get the value of last solution */
    printf("%s equation:\n"
           "---> Newton solver increased the residual so backtrack and exist ...\n",name);
    fflush(stdout);
    sbh_backtrack(grid,name);
    return stop_backtrack;
  }
  
  if (!stop_max)
  {
    printf("%s equation:\n"
           "---> Newton solver reached maximum step number so existing ...\n",name);
    fflush(stdout);
    return stop_max;
  }
  
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    double res      = patch->solving_man->Frms;/* current residual */
    double res_i    = patch->solving_man->settings->Frms_i;/* initial residual */
    
    /* since one of them is enough to continue */
    if (res > res_d && res > res_fac*res_i)
    {
      stop_res = 1;
      break;
    }
    
  }
  
  if (!stop_res)  
  {
    printf("%s equation:\n"
           "---> Newton solver satisfies demanding residual so existing ...\n",name);
    fflush(stdout);
    return stop_res;
  }
  
  return stop;
}

/* backtrack to restore to the last solution */
static void sbh_backtrack(Grid_T *const grid,const char *const name)
{
  const unsigned npatch = grid->np;
  unsigned p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    Field_T *f      = patch->pool[Ind(name)];
    double *v = f->v;
    const double *last_sol = patch->solving_man->settings->last_sol;
    unsigned ijk;
    
    free_coeffs(f);
    for(ijk = 0; ijk < patch->nn; ++ijk)
      v[ijk] = last_sol[ijk];
    
    sbh_SolveEqs_FieldUpdate(patch,name);
  }
  
}

/* updating field after they were solved */
void sbh_SolveEqs_FieldUpdate(Patch_T *const patch,const char *const name)
{
  if (!strcmp(name,"psi"))
  {
    sbh_update_derivative_psi(patch);
  }
  else if (!strcmp(name,"eta"))
  {
    sbh_update_derivative_eta(patch);
  }
  else if (!strcmp(name,"B0_U0"))
  {
    sbh_update_Beta_U0(patch);
    sbh_update_derivative_Beta_U0(patch);
    sbh_update_psi10A_UiUj(patch);
  }
  else if (!strcmp(name,"B0_U1"))
  {
    sbh_update_Beta_U1(patch);
    sbh_update_derivative_Beta_U1(patch);
    sbh_update_psi10A_UiUj(patch);
  }
  else if (!strcmp(name,"B0_U2"))
  {
    sbh_update_Beta_U2(patch);
    sbh_update_derivative_Beta_U2(patch);
    sbh_update_psi10A_UiUj(patch);
  }
  
}

static void sbh_XCTS_fill_db_eqs(sEquation_T ***const field_eq, 
                          sEquation_T ***const bc_eq, 
                          sEquation_T ***const jacobian_field_eq,
                          sEquation_T ***const jacobian_bc_eq)
{
  /* adding field and boundary condition equations of XCTS formalism */
  *field_eq      = init_eq();
  *bc_eq         = init_eq();
  *jacobian_field_eq   = init_eq();
  *jacobian_bc_eq      = init_eq();

  /* psi equations */
  add_eq(field_eq,sbh_eq_psi,"eq_psi");
  add_eq(bc_eq   ,sbh_bc_psi,"bc_psi");
  add_eq(jacobian_field_eq,sbh_jacobian_eq_psi,"jacobian_eq_psi");
  add_eq(jacobian_bc_eq   ,sbh_jacobian_bc_psi,"jacobian_bc_psi");
  
  /* eta equations */
  add_eq(field_eq,sbh_eq_eta,"eq_eta");
  add_eq(bc_eq   ,sbh_bc_eta,"bc_eta");
  add_eq(jacobian_field_eq,sbh_jacobian_eq_eta,"jacobian_eq_eta");
  add_eq(jacobian_bc_eq   ,sbh_jacobian_bc_eta,"jacobian_bc_eta");
  
  /* shift equations, remember we solve for B0_U? rather than Beta_U? */
  /* Beta_U0 equations */
  add_eq(field_eq,sbh_eq_Beta_U0,"eq_B0_U0");
  add_eq(bc_eq   ,sbh_bc_Beta_U0,"bc_B0_U0");
  add_eq(jacobian_field_eq,sbh_jacobian_eq_Beta_U0,"jacobian_eq_B0_U0");
  add_eq(jacobian_bc_eq   ,sbh_jacobian_bc_Beta_U0,"jacobian_bc_B0_U0");
  
  /* Beta_U1 equations */
  add_eq(field_eq,sbh_eq_Beta_U1,"eq_B0_U1");
  add_eq(bc_eq   ,sbh_bc_Beta_U1,"bc_B0_U1");
  add_eq(jacobian_field_eq,sbh_jacobian_eq_Beta_U1,"jacobian_eq_B0_U1");
  add_eq(jacobian_bc_eq   ,sbh_jacobian_bc_Beta_U1,"jacobian_bc_B0_U1");
  
  /* Beta_U2 equations */
  add_eq(field_eq,sbh_eq_Beta_U2,"eq_B0_U2");
  add_eq(bc_eq   ,sbh_bc_Beta_U2,"bc_B0_U2");
  add_eq(jacobian_field_eq,sbh_jacobian_eq_Beta_U2,"jacobian_eq_B0_U2");
  add_eq(jacobian_bc_eq   ,sbh_jacobian_bc_Beta_U2,"jacobian_bc_B0_U2");
  
}