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

/* initialize Solve Equations struct */
Solve_Equations_T *init_solve_equations(Grid_T *const grid)
{
  Solve_Equations_T *solve = calloc(1,sizeof(*solve));
  pointerEr(solve);
  
  solve->grid    = grid;
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
