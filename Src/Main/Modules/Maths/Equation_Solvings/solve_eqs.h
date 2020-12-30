#include "core_lib.h"
#include "error_handling_lib.h"
#include "maths_equation_solvings_lib.h"
#include "utilities_lib.h"

typedef int fSolve_T (Solve_Equations_T *const SolveEqs);

int solve_eqs(Solve_Equations_T *const SolveEqs);
int ddm_schur_complement(Solve_Equations_T *const SolveEqs);
Solve_Equations_T *init_solve_equations(Grid_T *const grid);
void free_solve_equations(Solve_Equations_T *solve);
Grid_T *get_grid_solve_equations(Solve_Equations_T *const solve);
void add_special_grid_solve_equations(Grid_T *const grid,const char *const name, Solve_Equations_T *const solve);
double get_relaxation_factor_solve_equations(Solve_Equations_T *const solve);
void sync_patch_pools(const Grid_T*const latest_grid,Solve_Equations_T *const solve);
static int default_stop_criteria_solve_equations(Grid_T *const grid,const char *const name);

