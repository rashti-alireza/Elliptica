#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "maths_calculus_lib.h"
#include "coordinates_lib.h"
#include "maths_approximation_lib.h"
#include "memory_managing_lib.h"
#include "maths_solvers_lib.h"
#include "maths_general_lib.h"
#include "maths_linear_algebra_lib.h"

int ddm_schur_complement(Grid_T *const grid);
static void preparing_ingredients(Grid_T *const grid);
static void make_map_and_inv(Patch_T *const patch);

