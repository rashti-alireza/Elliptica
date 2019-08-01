#include "core_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"
#include "memory_managing_lib.h"
#include "utilities_lib.h"
#include "maths_calculus_lib.h"

void populate_free_data(Grid_T *const grid);
static void _gammas(Grid_T *const grid);
static void free_conformal_metric_derivatives(Patch_T *const patch);
static void preparing_conformal_metric_derivatives(Patch_T *const patch);
static void _Gamma(Grid_T *const grid);
static void _dGamma(Grid_T *const grid);



