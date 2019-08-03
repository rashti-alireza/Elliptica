#include "core_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"
#include "memory_managing_lib.h"
#include "utilities_lib.h"
#include "maths_calculus_lib.h"

void bbn_populate_free_data(Grid_T *const grid);
static void _gammas(Grid_T *const grid);
static void free_conformal_metric_derivatives(Patch_T *const patch);
static void preparing_conformal_metric_derivatives(Patch_T *const patch);
static void _Gamma(Grid_T *const grid);
static void _dGamma(Grid_T *const grid);
static void _Ricci(Grid_T *const grid);
static void populate_KS_trKij(Patch_T *const patch);
static void populating_KSGamma(Patch_T *const patch);
static void populate_KSgammas_KSalpha_KSBeta(Patch_T *const patch);
static void partial_derivative_KSBeta(Patch_T *const patch);
static void free_KSfields(Patch_T *const patch);
static void tr_KSKij(Grid_T *const grid);





