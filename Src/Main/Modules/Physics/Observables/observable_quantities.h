#include "core_lib.h"
#include "memory_managing_lib.h"
#include "physics_observables_lib.h"
#include "utilities_lib.h"
#include "coordinates_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"

Observable_T *init_observable(void *grid);
void plan_observable(Observable_T *const obs);
static void populate_observable_BBN_CS(Observable_T *const obs);
void free_observable(Observable_T *obs);
void populate_ADM_momentums_integrand_PdS_GdV(const Observable_T *const obs);
