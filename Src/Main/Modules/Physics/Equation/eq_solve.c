/*
// Alireza Rashti
// December 2020
*/

/* issue solve equations and related */

#include "eq_solve.h"

/* solve elliptic equations fields.
//
// algorithm: 
// ==========
// setup Solve_Equations_T
// backup fields
// solve fields
// update fields. */
void eq_solve_elliptic_equation(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid = phys->grid;
  
  /* populating solution managment */
  initialize_solving_man
    (grid,eq_global_field_eq,eq_global_bc_eq,eq_global_jacobian_field_eq,eq_global_jacobian_bc_eq,P_);
  
  /* solving equation(s) */
  Solve_Equations_T *SolveEqs = init_solve_equations(grid);
  if (eq_field_update)  SolveEqs->FieldUpdate  = eq_field_update;
  if (eq_source_update) SolveEqs->SourceUpdate = eq_source_update;
  if (eq_stop_criteria) SolveEqs->StopCriteria = eq_stop_criteria;
  SolveEqs->residual_suffix = EQ_Residual_Suffix;
  
  /* set region of each equation */
  Grid_T **lgrid = set_equation_grid(phys,SolveEqs);

  /* populate backup fields */
  backup_fields(phys);

  const int max_iter = Pgeti("solve_Max_Iteration");
  int iter = 0;
  while (iter < max_iter)
  {
    printf(Pretty0"Iteration      = %d / %d\n",iter+1,max_iter);
    printf(Pretty0"Equation(s)    = %s\n",SolveEqs->solving_order);
    printf(Pretty0"Newton step(s) = %d\n",Pgeti("solve_Max_Newton_Step"));
    fflush(stdout);
    
    /* solve equations */
    solve_eqs(SolveEqs);
    
    /* updating the fields using relaxed scheme */
    if (Pcmps(P_"update_method","relaxed_scheme"))
      update_fields_relaxed_scheme(phys);
    else
      Error0(NO_OPTION);

    /* update backup fields.
    // for optimization purposes if it is only a single run, 
    // don't need to backup. in fact, a single run is the most used 
    // case in the initial data construction. */
    if (max_iter > 1)
      backup_fields(phys);
      
    /* study the solution */
    if (Pcmps(P_"elliptic_test","yes"))
    {
      calculate_equation_residual(SolveEqs);
      if (eq_analyze_solution) 
        eq_analyze_solution(phys,iter);
    }
    
    ++iter;
  }
  
  /* calculate the field residual for diagnostic purposes */
  calculate_equation_residual(SolveEqs);

  /* free */
  free_solve_equations(SolveEqs);
  free_equation_grid(lgrid);

  /* free data base of equations */
  free_db_eqs(eq_global_field_eq);
  free_db_eqs(eq_global_bc_eq);
  free_db_eqs(eq_global_jacobian_field_eq);
  free_db_eqs(eq_global_jacobian_bc_eq);
  eq_global_field_eq = 0;
  eq_global_bc_eq    = 0;
  eq_global_jacobian_field_eq = 0;
  eq_global_jacobian_bc_eq    = 0;
  
  FUNC_TOC
}

/* ->: return list of added special grid for solve, ends with NULL.
// set region(grid) of each equation for SolveEqs 
// according to the parameter value. ex: Eq_phi = XCTS_curve,NS 
// => equation region is all the patches which cover NS */
static Grid_T **set_equation_grid(Physics_T *const phys,
                                   Solve_Equations_T *const SolveEqs)
{
  FUNC_TIC
  
  /* what fields to be solved */
  Grid_T **lgrid = 0;
  char **fields_name = read_separated_items_in_string(Pgets("solve_Order"),',');
  char par_val[STR_LEN];
  char par_name[STR_LEN];
  char *aux = 0;
  Uint Nl = 0;
  Uint p,Np,fn;
  
  fn = 0;
  while (fields_name[fn])
  {
    sprintf(par_name,P_"%s",fields_name[fn]);
    sprintf(par_val,"%s",Pgets(par_name));
    aux = strchr(par_val,',');
    if (!aux)
    {
      Errors("No region specified in parameter '%s'!",par_name);
    }
    aux++;
    
    /* collect the pertinent patches, unless it is every patches (.*) */
    if (!strstr(aux,".*"))
    {
      lgrid = realloc(lgrid,(Nl+2)*sizeof(*lgrid));IsNull(lgrid);
      Grid_T *eq_grid   = calloc(1,sizeof(*eq_grid));IsNull(eq_grid);
      Patch_T **patches = collect_patches(phys->grid,aux,&Np);
      lgrid[Nl]   = eq_grid;
      lgrid[Nl+1] = 0;
      Nl++;
      
      eq_grid->kind = phys->grid->kind;
      eq_grid->gn   = phys->grid->gn;
      eq_grid->np   = Np;
      
      /* populate patches */
      eq_grid->patch = calloc(Np+1,sizeof(*eq_grid->patch));
      IsNull(eq_grid->patch);
      for(p = 0; p < Np; ++p)
      {
        Patch_T *patch = patches[p];
        eq_grid->patch[p] = calloc(1,sizeof(*eq_grid->patch[p]));
        IsNull(eq_grid->patch[p]);
        /* note, for the following all of the pointers inside the structures 
        // will be equal, since this is not a deep copy. */
        eq_grid->patch[p][0]    = patch[0];
        eq_grid->patch[p]->pn   = p;
        eq_grid->patch[p]->grid = eq_grid;
        /* the following needs to be constructed from scratch */
        eq_grid->patch[p]->interface = 0;
        eq_grid->patch[p]->solving_man = calloc(1,sizeof(*eq_grid->patch[p]->solving_man));
        IsNull(eq_grid->patch[p]->solving_man);
        eq_grid->patch[p]->solving_man[0] = patch->solving_man[0];
        eq_grid->patch[p]->solving_man->patch = eq_grid->patch[p];
        eq_grid->patch[p]->solving_man->jacobian = 0;
        eq_grid->patch[p]->solving_man->nj       = 0;
        eq_grid->patch[p]->solving_man->method->Schur_Complement = 0;
        eq_grid->patch[p]->solving_man->method->SchurC = 0;
      }
      /* now let's fill up eq_grid->patch[?]->interface */
      realize_interfaces(eq_grid);
      
      /* add to special grid for solve */
      add_special_grid_solve_equations(eq_grid,fields_name[fn],SolveEqs);
      
      Free(patches);
    }
    fn++;
  }
  
  /* free */
  free_2d(fields_name);
  
  FUNC_TOC
  return lgrid;
}

/* free only those thing we allocate for equation grid */
static void free_equation_grid(Grid_T **lgrid)
{
  if (!lgrid)
    return;
  
  Uint ng = 0;
  while (lgrid[ng])
  {
    Grid_T *grid = lgrid[ng];
    Uint p;
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      
      free_patch_interface(patch);
      //free_patch_SolMan_jacobian(patch);
      free_patch_SolMan_method_Schur(patch);
      free(patch->solving_man);
    }
    free_2d_mem(grid->patch,grid->np);
    Free(grid);
    ng++;
  }
  
  Free(lgrid);
}

/* saving the field with the given name for iterative purposes. */
static void backup_fields(Physics_T *const phys)
{
  Grid_T *const grid = phys->grid;
  const Uint npatch  = grid->np;
  char **field_name = read_separated_items_in_string(Pgets("solve_Order"),',');
  Uint fn;
  
  fn = 0;
  while (field_name[fn])
  {
    char field_backup[STR_LEN];
    sprintf(field_backup,P_Backup_"%s",field_name[fn]);
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (Uint p = 0; p < npatch; ++p)
    {
      Patch_T *patch = grid->patch[p];
  
      /* if field not defined in this patch */
      if (_Ind(field_name[fn]) < 0) continue;
        
      Field_T *f      = patch->fields[Ind(field_name[fn])];
      Field_T *f_bckp = patch->fields[Ind(field_backup)];
      
      /* if no field value defined in this patch */
      if (!f->v) continue;
        
      empty_field(f_bckp);
      f_bckp->v = alloc_double(patch->nn);
      
      FOR_ALL_ijk
        f_bckp->v[ijk] = f->v[ijk];
   
    }
    fn++;
  }
  
  free_2d(field_name);
}

/* using initial fields and solved fields then 
// update every thing using relaxed scheme: new = W1*solved+W2*old. 
// NOTE: it updates source and fields too. */ 
static void update_fields_relaxed_scheme(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid = phys->grid;
  const Uint npatch  = grid->np;
  char **field_name = read_separated_items_in_string(Pgets("solve_Order"),',');
  Uint f;
  
  /* update all of the fields were solved according to the given weight */
  f = 0;
  while (field_name[f])
  {
    const char *field_new = field_name[f];
    char field_old[STR_LEN];
    char par[STR_LEN];
    
    sprintf(field_old,P_Backup_"%s",field_new);
    sprintf(par,P_"update_weight_%s",field_new);
    double W1  = Pgetd(par);
    double W2  = 1.-W1;
    
    /* don't bother if w = 1. (optimize) */
    if (EQL(W1,1.))
    {
      f++;
      continue;
    }
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (Uint p = 0; p < npatch; ++p)
    {
      Patch_T *patch = grid->patch[p];
      
      /* if the field is not defined in this patch */
      if (_Ind(field_new) < 0)
        continue;
        
      Field_T *f_old  = patch->fields[Ind(field_old)];
      Field_T *f_new  = patch->fields[Ind(field_new)];
      
      /* if the field value is not defined in this patch */
      if (!f_new->v)
        continue;
        
      free_coeffs(f_new);
      
      FOR_ALL_ijk
        f_new->v[ijk] = W1*f_new->v[ijk]+W2*f_old->v[ijk];
      
      if (eq_field_update)
        eq_field_update(patch,field_new);
    }
    f++;
  }
  if (eq_source_update)
    eq_source_update(grid,0);
  
  /* free names */
  free_2d(field_name);
  
  FUNC_TOC
}

