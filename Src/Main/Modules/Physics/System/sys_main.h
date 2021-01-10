#include "sys_header.h"
#include "physics_star_lib.h"


int sys_main(Physics_T *const phys);
static int tune_system_ADM_momenta(Physics_T *const phys);
static int initialize_fields(Physics_T *const phys);
static int set_system_params(Physics_T *const phys);
static int add_system_fields(Physics_T *const phys);


