#include "obs_header.h"

int observe(Physics_T *const phys,const char *const sq,
            const char *const method,double *const ret);
void obs_calculate(Observe_T *const obs);
static void free_obs(Observe_T *obs);
static int set_observe_params(Physics_T *const phys);
static int add_observe_fields(Physics_T *const phys);


