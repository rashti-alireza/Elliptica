#include "core_lib.h"
#include "coordinates_lib.h"
#include "macros_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_analytic_lib.h"
#include "maths_calculus_lib.h"
#include "maths_solvers_lib.h"
#include "maths_linear_algebra_lib.h"

#define DO 1
#define NOT_DO 0

static int test_solver_umfpack_di(void);
static int test_solver_series_umfpack_di(void);
void solver_tests(void);
