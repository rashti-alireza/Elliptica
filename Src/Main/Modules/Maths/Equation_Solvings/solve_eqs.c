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
  if (strcmp_i(GetParameterS_E("Solving_Method"),"DDM_Schur_Complement"))
    fsolve = ddm_schur_complement;
  else
    abortEr_s("No such method \"%s\" defined for this function.\n",
      GetParameterS("Solving_Method"));
  
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
  const double res_d    = GetParameterD_E("Solving_Residual");/* desired residual */
  const int max_step    = GetParameterI_E("Solving_Max_Number_of_Newton_Step");
  const unsigned npatch = grid->np;
  unsigned p;

  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    double res      = patch->solving_man->Frms;/* current residual */
    int solver_step = patch->solving_man->settings->solver_step;/* iteration number */

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
  double factor = GetParameterD("Solving_Relaxation_Factor");/* relaxation factor */
  char par[400] = {'\0'};
  
  if (factor == DBL_MAX)/* if no such parameter defined */
    factor = 1;
  
  if (f_name)
  {
    sprintf(par,"Solving_Relaxation_Factor_%s",f_name);
    double factor2 = GetParameterD(par);
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
