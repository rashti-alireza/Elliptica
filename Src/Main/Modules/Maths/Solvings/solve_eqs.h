#include "core_lib.h"
//#include "utilities_lib.h"
#include "error_handling_lib.h"
//#include "maths_calculus_lib.h"
//#include "coordinates_lib.h"
//#include "maths_approximation_lib.h"
//#include "memory_managing_lib.h"
//#include "maths_solvers_lib.h"

typedef int fSolve_T (Grid_T *const grid);

int solve_eqs(Grid_T *const grid);
int parallel_patch_method (Grid_T *const grid);
