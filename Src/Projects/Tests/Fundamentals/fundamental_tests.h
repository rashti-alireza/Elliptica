#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "maths_calculus_lib.h"
#include "coordinates_lib.h"
#include "maths_approximation_lib.h"
#include "memory_managing_lib.h"
#include "maths_solvers_lib.h"
#include "maths_solvings_lib.h"
#include "maths_general_lib.h"

#define OpenMP(x) _Pragma(#x)

int Fundamental_Tests(void);
Grid_T *fundamental_tests_make_grid(void);
