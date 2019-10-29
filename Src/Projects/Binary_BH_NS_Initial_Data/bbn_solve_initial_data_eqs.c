/*
// Alireza Rashti
// August 2019
*/

#include "bbn_solve_initial_data_eqs.h"

/* solving initial data equations for the given grid */
void bbn_solve_initial_data_eqs(Grid_T *const grid)
{
  pr_line_custom('='); 
  printf("{ Solving initial data equations for Binary BH and NS ...\n");
  
  sEquation_T **field_eq/* field equation */,
              **bc_eq/* B.C. for the field */,
              **jacobian_field_eq/* jacobian for field equation */,
              **jacobian_bc_eq/* jacobian for B.C. */;

  /* filling db equations of XCTS */
  bbn_XCTS_fill_db_eqs(&field_eq,&bc_eq,&jacobian_field_eq,&jacobian_bc_eq);
  
  /* populating solution managment */
  initialize_solving_man(grid,field_eq,bc_eq,jacobian_field_eq,jacobian_bc_eq);
  
  /* solving equation(s) */
  Solve_Equations_T *SolveEqs = init_solve_equations(grid);
  SolveEqs->solving_order = GetParameterS_E("Solving_Order");
  SolveEqs->FieldUpdate   = bbn_SolveEqs_FieldUpdate;
  SolveEqs->SourceUpdate  = bbn_SolveEqs_SourceUpdate;
  SolveEqs->StopCriteria  = bbn_stop_criteria;
  
  Grid_T *phi_grid = bbn_phi_grid(grid);/* phi needed to be solved only in NS */
  add_special_grid_solve_equations(phi_grid,"phi",SolveEqs);
  
  /* saving the field being solved for relaxation scheme purposes */
  save_fields(grid);
  
  const unsigned max_iter = (unsigned)GetParameterI_E("Solving_Max_Number_of_Iteration");
  unsigned iter = 0;
  while (iter < max_iter)
  {
    /* some prints */
    pr_line_custom('=');
    printf("{ Iteration %d For Solving XCTS Equations at a Fixed Resolution ...\n",iter);
    printf("        |---> %s Equations ...\n",SolveEqs->solving_order);
    
    /* test if jacobian of equations written correctly */
    if (0) test_Jacobian_of_equations(SolveEqs);
    
    /* solve equations */
    solve_eqs(SolveEqs);
    
    /* some prints */
    printf("} Iteration %d For Solving XCTS Equations at a Fixed Resolution ==> Done.\n",iter);
    pr_clock();
    
    ++iter;
  }
  pr_line_custom('=');
  
  /* free SolveEqs and phi grid */
  free_solve_equations(SolveEqs);
  bbn_free_phi_grid(phi_grid);
  
  /* updating the fields using relaxed scheme */
  update_fields_relaxed_scheme(grid);
  
  /* free data base of equations */
  free_db_eqs(field_eq);
  free_db_eqs(bc_eq);
  free_db_eqs(jacobian_field_eq);
  free_db_eqs(jacobian_bc_eq);
  
  printf("} Solving initial data equations for Binary BH and NS ==> Done.\n");
  pr_clock();
  pr_line_custom('='); 
}

/* using initial fields and solved fields update every thing 
// using relaxed scheme, note, we also update source fields */ 
static void update_fields_relaxed_scheme(Grid_T *const grid)
{
  const unsigned npatch = grid->np;
  const char *const solving_order = GetParameterS_E("Solving_Order");
  const double W1  = GetParameterD_E("Solving_Field_Update_Weight");
  const double W2  = 1-W1;
  char **field_name;
  unsigned p,nf,f;
  
  /* no need to update it W1 = 1 */
  if (EQL(W1,1))
    return;
  
  field_name = get_solving_field_name(solving_order,&nf);
  
  /* update all of the fields were solved according to the given weight */
  for (f = 0; f < nf; ++f)
  {
    const char *field_new = field_name[f];
    char field_old[100];
    sprintf(field_old,"OLD_%s",field_new);
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < npatch; ++p)
    {
      Patch_T *patch = grid->patch[p];
      unsigned nn    = patch->nn;
      unsigned ijk;
      
      /* if the field is not defined in this patch */
      if (_Ind(field_new) < 0)
        continue;
        
      Field_T *f_old  = patch->pool[Ind(field_old)];
      Field_T *f_new  = patch->pool[Ind(field_new)];
      free_coeffs(f_new);
      
      for (ijk = 0; ijk < nn; ++ijk)
        f_new->v[ijk] = W1*f_new->v[ijk]+W2*f_old->v[ijk];
      
      bbn_SolveEqs_FieldUpdate(patch,field_new);
    }
  }/* end of for (f = 0; f < nf; ++f) */
  Tij_IF_CTS_psi6Sources(grid);
  
  /* free names */
  free_2d_mem(field_name,nf);
  
}

/* saving the field with the given name for iterative purposes. */
static void save_fields(Grid_T *const grid)
{
  const unsigned npatch = grid->np;
  const char *const solving_order = GetParameterS_E("Solving_Order");  
  char **field_name;
  unsigned p,nf,f;
  
  field_name = get_solving_field_name(solving_order,&nf);
  
  /* save all of the fields being solved */
  for (f = 0; f < nf; ++f)
  {
    const char *fname0 = field_name[f];
    char fname_old[100];
    sprintf(fname_old,"OLD_%s",fname0);
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < npatch; ++p)
    {
      Patch_T *patch = grid->patch[p];
      unsigned nn    = patch->nn;
      unsigned ijk;
      
      /* if no field defined in this patch */
      if (_Ind(fname0) < 0)
        continue;
        
      Field_T *f_old  = patch->pool[Ind(fname_old)];
      Field_T *f0     = patch->pool[Ind(fname0)];
      free_coeffs(f_old);
      
      for (ijk = 0; ijk < nn; ++ijk)
        f_old->v[ijk] = f0->v[ijk];
   
    }
  }/* end of for (f = 0; f < nf; ++f)*/
}

/* stop criteria for solver, namely, if some conditions satisfied, 
// it stops. one can give specific criteria according to the given field's name.
// ->return value: 0 means stop, 1 means continue to solve */
int bbn_stop_criteria(Grid_T *const grid,const char *const name)
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
    bbn_backtrack(grid,name);
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
static void bbn_backtrack(Grid_T *const grid,const char *const name)
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
    
    bbn_SolveEqs_FieldUpdate(patch,name);
  }
  
}

/* updating sources after field is solved */
void bbn_SolveEqs_SourceUpdate(Grid_T *const grid,const char *const name)
{
  Tij_IF_CTS_psi6Sources(grid);
  
  /*if (!strcmp(name,"phi"))
  {
    unsigned p;
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      
      if (!IsItNSPatch(patch))
        continue;
        
      Tij_IF_CTS_enthalpy(patch);
      bbn_update_derivative_enthalpy(patch);
      bbn_update_rho0(patch);
      bbn_update_derivative_rho0(patch);
    }
  }*/
  
  UNUSED(name);
}

/* updating field after they were solved */
void bbn_SolveEqs_FieldUpdate(Patch_T *const patch,const char *const name)
{
  if (!strcmp(name,"phi"))
  {
    bbn_update_derivative_phi(patch);
    //Tij_IF_CTS_enthalpy(patch);
    //bbn_update_derivative_enthalpy(patch);
  }
  else if (!strcmp(name,"psi"))
  {
    bbn_update_derivative_psi(patch);
  }
  else if (!strcmp(name,"eta"))
  {
    bbn_update_derivative_eta(patch);
  }
  else if (!strcmp(name,"B0_U0"))
  {
    bbn_update_Beta_U0(patch);
    bbn_update_derivative_Beta_U0(patch);
    bbn_update_psi10A_UiUj(patch);
  }
  else if (!strcmp(name,"B0_U1"))
  {
    bbn_update_Beta_U1(patch);
    bbn_update_derivative_Beta_U1(patch);
    bbn_update_psi10A_UiUj(patch);
  }
  else if (!strcmp(name,"B0_U2"))
  {
    bbn_update_Beta_U2(patch);
    bbn_update_derivative_Beta_U2(patch);
    bbn_update_psi10A_UiUj(patch);
  }
  
}

/* at Euler's equation for phi field we need to confine the grid 
// to only NS patches. this function collect all of NS patches and
// make a new grid out of them and return it.
// NOTE, WE ARE NOT ALLOWED TO ADD FIELD TO THIS GRID, generally
// since, we want to keep track of the value of field, we do not
// deep copy in structures.
// ->return value: set of all NS patches as a separate grid */
static Grid_T *bbn_phi_grid(Grid_T *const grid)
{
  Grid_T *phi_grid = 0;
  
  if (strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
  {
    phi_grid       = alloc_grid();
    phi_grid->kind = grid->kind;
    phi_grid->gn   = grid->gn;
    /* NS at left composed of 6 cubed spherical + 1 Cartesian: */
    phi_grid->np = 7;
    bbn_phi_grid_CS(phi_grid,grid);
  }
  else
    abortEr(NO_JOB);
  
  return phi_grid;
}

/* free only those thing we allocate for this particular grid */
static void bbn_free_phi_grid(Grid_T *grid)
{
  if (!grid)
    return;
  
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    free_patch_interface(patch);
    free_patch_SolMan_jacobian(patch);
    free_patch_SolMan_method_Schur(patch);
    free(patch->solving_man);
  }
  free_2d_mem(grid->patch,grid->np);
  free(grid);
}

/* confining the whole grid to only NS grid for phi in cubed spherical coords. */
static void bbn_phi_grid_CS(Grid_T *const phi_grid,Grid_T *const grid)
{
  unsigned p,i;
  
  /* NS at left composed of 6 cubed spherical + 1 Cartesian,
  // and 1 more to be Null = 8 */
  phi_grid->patch = calloc(8,sizeof(*phi_grid->patch));
  pointerEr(phi_grid->patch);
  
  i = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (IsItNSPatch(patch))/* since we have only one NS we can use IsItNSPatch */
    {
      phi_grid->patch[i] = calloc(1,sizeof(*phi_grid->patch[i]));
      pointerEr(phi_grid->patch[i]);
      /* note, for the following all of the pointers inside the structures 
      // will be equal, since this is not a deep copy. */
      phi_grid->patch[i][0]    = patch[0];
      phi_grid->patch[i]->pn   = i;
      phi_grid->patch[i]->grid = phi_grid;
      phi_grid->nn            += patch->nn;
      /* the following needs to be constructed from scratch */
      phi_grid->patch[i]->interface = 0;
      phi_grid->patch[i]->solving_man = calloc(1,sizeof(*phi_grid->patch[i]->solving_man));
      pointerEr(phi_grid->patch[i]->solving_man);
      phi_grid->patch[i]->solving_man[0] = patch->solving_man[0];
      phi_grid->patch[i]->solving_man->patch = phi_grid->patch[i];
      phi_grid->patch[i]->solving_man->jacobian = 0;
      phi_grid->patch[i]->solving_man->nj       = 0;
      phi_grid->patch[i]->solving_man->method->Schur_Complement = 0;
      phi_grid->patch[i]->solving_man->method->SchurC = 0;
      
      ++i;
    }
  }
  
  assert(i==7);
  
  /* no let's fill up phi_grid->patch[?]->interface */
  realize_geometry(phi_grid);      
}

static void bbn_XCTS_fill_db_eqs(sEquation_T ***const field_eq, 
                          sEquation_T ***const bc_eq, 
                          sEquation_T ***const jacobian_field_eq,
                          sEquation_T ***const jacobian_bc_eq)
{
  /* adding field and boundary condition equations of XCTS formalism */
  *field_eq      = init_eq();
  *bc_eq         = init_eq();
  *jacobian_field_eq   = init_eq();
  *jacobian_bc_eq      = init_eq();

 /* phi equations */
  add_eq(field_eq,bbn_eq_phi,"eq_phi");
  add_eq(bc_eq   ,bbn_bc_phi,"bc_phi");
  add_eq(jacobian_field_eq,bbn_jacobian_eq_phi,"jacobian_eq_phi");
  add_eq(jacobian_bc_eq   ,bbn_jacobian_bc_phi,"jacobian_bc_phi");
  
  /* psi equations */
  add_eq(field_eq,bbn_eq_psi,"eq_psi");
  add_eq(bc_eq   ,bbn_bc_psi,"bc_psi");
  add_eq(jacobian_field_eq,bbn_jacobian_eq_psi,"jacobian_eq_psi");
  add_eq(jacobian_bc_eq   ,bbn_jacobian_bc_psi,"jacobian_bc_psi");
  
  /* eta equations */
  add_eq(field_eq,bbn_eq_eta,"eq_eta");
  add_eq(bc_eq   ,bbn_bc_eta,"bc_eta");
  add_eq(jacobian_field_eq,bbn_jacobian_eq_eta,"jacobian_eq_eta");
  add_eq(jacobian_bc_eq   ,bbn_jacobian_bc_eta,"jacobian_bc_eta");
  
  /* shift equations, remember we solve for B0_U? rather than Beta_U? */
  /* Beta_U0 equations */
  add_eq(field_eq,bbn_eq_Beta_U0,"eq_B0_U0");
  add_eq(bc_eq   ,bbn_bc_Beta_U0,"bc_B0_U0");
  add_eq(jacobian_field_eq,bbn_jacobian_eq_Beta_U0,"jacobian_eq_B0_U0");
  add_eq(jacobian_bc_eq   ,bbn_jacobian_bc_Beta_U0,"jacobian_bc_B0_U0");
  
  /* Beta_U1 equations */
  add_eq(field_eq,bbn_eq_Beta_U1,"eq_B0_U1");
  add_eq(bc_eq   ,bbn_bc_Beta_U1,"bc_B0_U1");
  add_eq(jacobian_field_eq,bbn_jacobian_eq_Beta_U1,"jacobian_eq_B0_U1");
  add_eq(jacobian_bc_eq   ,bbn_jacobian_bc_Beta_U1,"jacobian_bc_B0_U1");
  
  /* Beta_U2 equations */
  add_eq(field_eq,bbn_eq_Beta_U2,"eq_B0_U2");
  add_eq(bc_eq   ,bbn_bc_Beta_U2,"bc_B0_U2");
  add_eq(jacobian_field_eq,bbn_jacobian_eq_Beta_U2,"jacobian_eq_B0_U2");
  add_eq(jacobian_bc_eq   ,bbn_jacobian_bc_Beta_U2,"jacobian_bc_B0_U2");
  
}
