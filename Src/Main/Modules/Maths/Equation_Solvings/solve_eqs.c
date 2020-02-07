/*
// Alireza Rashti
// August 2018
*/

#include "solve_eqs.h"

/* solving equations
// ->return value: EXIT_SUCCESS
*/
int solve_eqs(Solve_Equations_T *const SolveEqs)
{
  fSolve_T *fsolve = 0;
  
  /* choosing solving method */
  if (strcmp_i(Pgets_E("Solving_Method"),"DDM_Schur_Complement"))
    fsolve = ddm_schur_complement;
  else
    abortEr_s("No such method \"%s\" defined for this function.\n",
      PgetsEZ("Solving_Method"));
  
  /* call the specific solving method */
  fsolve(SolveEqs);
  
  return EXIT_SUCCESS;
}

// ->return value: 0 means stop, 1 means continue; */
int default_stop_criteria_solve_equations(Grid_T *const grid,const char *const name)
{
  int stop = 1;
  int stop_max = 1;
  int stop_res = 0;
  const double res_d    = Pgetd_E("Solving_Residual");/* desired residual */
  const int max_step    = Pgeti_E("Solving_Max_Number_of_Newton_Step");
  const unsigned npatch = grid->np;
  unsigned p;

  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    double res      = patch->solving_man->Frms;/* current residual */
    int solver_step = patch->solving_man->settings->solver_step;/* iteration number */

    /* NOTE: due to the break command, the order of ifs are important */
    
    /* note: all patches have same solver_step */
    if (solver_step >= max_step)
    {
      stop_max = 0;
      break;
    }
    
    /* since one of them is enough to continue */
    if (res > res_d)
    {
      stop_res = 1;
      break;
    }
  }

  if (!stop_max)
  {
    printf("%s equation:\n"
           "---> Newton solver reached maximum step number so existing ...\n",name);
    fflush(stdout);
    return stop_max;
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

/* ->return value : get the double value of relaxation factor 
// in relaxation scheme the default value is 1, 
// which means no relaxation at all. */
double get_relaxation_factor_solve_equations(Solve_Equations_T *const solve)
{
  const char *f_name = solve->field_name;
  double factor = Pgetd("Solving_Newton_Update_Weight");/* relaxation factor */
  char par[400] = {'\0'};
  
  if (factor == DBL_MAX)/* if no such parameter defined */
    factor = 1;
  
  if (f_name)
  {
    sprintf(par,"Solving_Newton_Update_Weight_%s",f_name);
    double factor2 = Pgetd(par);
    if (factor2 != DBL_MAX)
      factor = factor2;
  }
  
  return factor;
}

/* initialize Solve Equations struct */
Solve_Equations_T *init_solve_equations(Grid_T *const grid)
{
  Solve_Equations_T *solve = calloc(1,sizeof(*solve));
  pointerEr(solve);
  
  solve->grid         = grid;
  solve->StopCriteria = default_stop_criteria_solve_equations; 
  
  return solve;
}

/* free Solve Equations struct */
void free_solve_equations(Solve_Equations_T *solve)
{
  if (!solve)
    return;
  
  if (solve->Sgrid)
  {
    unsigned i = 0;
    while (solve->Sgrid[i])
    {
      free(solve->Sgrid[i]->name);
      free(solve->Sgrid[i]);
      ++i;
    }
    free(solve->Sgrid);
  }
  
  if(solve)
    free(solve);
}

/* if no function for GetGrid is defined, it gives 
// the default value which is the original grid */
Grid_T *get_grid_solve_equations(Solve_Equations_T *const solve)
{
  Grid_T *grid = solve->grid;/* default grid (the whole grid) */
  const char *const name = solve->field_name;
  
  if (solve->Sgrid)
  {
    unsigned i = 0;
    while (solve->Sgrid[i])
    {
      if (!strcmp(solve->Sgrid[i]->name,name))
        return solve->Sgrid[i]->sgrid;
      ++i;
    }
  }
  return grid;
}

/* give grid and name, add them to the Sgrid */
void add_special_grid_solve_equations(Grid_T *const grid,const char *const name, Solve_Equations_T *const solve)
{
  unsigned N = 0;
  if (solve->Sgrid)/* if there is any Sgrid then count */
    N = countf(solve->Sgrid);
  
  solve->Sgrid = realloc(solve->Sgrid,(N+2)*sizeof(*solve->Sgrid));
  pointerEr(solve->Sgrid);
  solve->Sgrid[N+1] = 0;/* terminated with null */
  solve->Sgrid[N] = calloc(1,sizeof(*solve->Sgrid[N]));
  pointerEr(solve->Sgrid[N]);
  solve->Sgrid[N]->sgrid = grid;
  solve->Sgrid[N]->name  = dup_s(name);
}

/* since each field might be solved in a special grid 
// so each patch of the grid shares the same pool but different
// interface structure; however, pointer to pool might get changes
// due to free or add field while solving and calculating some equations.
// these changes must be considered while moving to the next spcial gird. 
// this function sync all of the pool according to the latest grid. */
void sync_patch_pools(const Grid_T*const latest_grid,Solve_Equations_T *const solve)
{
  /* if there is no special grid so don't bother */
  if (!solve->Sgrid)
    return;
    
  unsigned p1;
  
  FOR_ALL_PATCHES(p1,latest_grid)
  {
    const Patch_T *patch1 = latest_grid->patch[p1];
    const char *name1 = strstr(patch1->name,"_");
    Grid_T *outdated_grid;
    Patch_T *patch2;
    const char *name2;
    unsigned i,p2;
    name1++;
    
    outdated_grid = solve->grid;/* default grid */
    FOR_ALL_PATCHES(p2,outdated_grid)
    {
      patch2 = outdated_grid->patch[p2];
      name2  = strstr(patch2->name,"_");
      name2++;
  
      if (!strcmp(name2,name1))
      {
        patch2->pool = patch1->pool;
        patch2->nfld = patch1->nfld;
        
        /* if patch1 has any jacobians */
        if (patch1->solving_man)
        {
          /* if patch2 need memory */
          if (!patch2->solving_man)
          {
            patch2->solving_man = calloc(1,sizeof(*patch2->solving_man));
            pointerEr(patch2->solving_man);
          }
          patch2->solving_man->jacobian = patch1->solving_man->jacobian;
          patch2->solving_man->nj       = patch1->solving_man->nj;
        }
        
        break;
      }
    }/* FOR_ALL_PATCHES(p2,outdated_grid) */
    
    /* for other special grid */
    i = 0;
    while (solve->Sgrid[i])
    {
      outdated_grid = solve->Sgrid[i]->sgrid;
      FOR_ALL_PATCHES(p2,outdated_grid)
      {
        patch2 = outdated_grid->patch[p2];
        name2  = strstr(patch2->name,"_");
        name2++;
    
        if (!strcmp(name2,name1))
        {
          patch2->pool = patch1->pool;
          patch2->nfld = patch1->nfld;
          /* if patch1 has any jacobians */
          if (patch1->solving_man)
          {
            /* if patch2 need memory */
            if (!patch2->solving_man)
            {
              patch2->solving_man = calloc(1,sizeof(*patch2->solving_man));
              pointerEr(patch2->solving_man);
            }
            patch2->solving_man->jacobian = patch1->solving_man->jacobian;
            patch2->solving_man->nj       = patch1->solving_man->nj;
          }
          break;
        }
      }/* FOR_ALL_PATCHES(p2,outdated_grid) */
      i++;
    }/* end of while (solve->Sgrid[i]) */
    
  }/* FOR_ALL_PATCHES(p1,latest_grid) */
  
}
