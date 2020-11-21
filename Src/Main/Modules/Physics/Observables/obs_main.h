#include "core_lib.h"
#include "physics_observables_lib.h"
#include "utilities_lib.h"

Observable_T *init_observable(void *grid,void (*plan_items)(struct OBSERVABLE_T *obs),void (*free_items)(struct OBSERVABLE_T *obs));
void plan_observable(Observable_T *const obs);
void free_observable(Observable_T *obs);


