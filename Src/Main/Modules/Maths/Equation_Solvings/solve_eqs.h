#include "core_lib.h"
#include "error_handling_lib.h"

typedef int fSolve_T (Grid_T *const grid);

int solve_eqs(Grid_T *const grid);
int ddm_schur_complement(Grid_T *const grid);

