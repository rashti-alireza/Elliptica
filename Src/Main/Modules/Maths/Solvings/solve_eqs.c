/*
// Alireza Rashti
// August 2018
*/

#include "solve_eqs.h"

/* solving equations
// ->return value: EXIT_SUCCESS
*/
int solve_eqs(Grid_T *const grid)
{
  fSolve_T *solve = 0;
  
  /* choosing solving method */
  if (strcmp_i(GetParameterS_E("Solving_Method"),"Parallel_Patch"))
    solve = parallel_patch_method;
  else
    abortEr_s("No such method \"%s\" defined for this function.\n",
      GetParameterS("Solving_Method"));
  
  /* call the specific solving method */
  solve(grid);
  
  return EXIT_SUCCESS;
}

