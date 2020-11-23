#include "core_lib.h"
#include "utilities_lib.h"
#include "physics_lib.h"
#include "physics_StressEnergyTensor_lib.h"
#include "physics_star_lib.h"

#define STR_LEN (999)

int physics(Physics_T *const phys,const cmd_T cmd,
            const char *const file, const int line);

Physics_T *
init_physics
 (
 Grid_T *const grid/* computation grid */,
 const Com_Obj_T type/* object type NS,BH,etc */
 );

void free_physics(Physics_T *obj);
static void set_phys_default_region(Physics_T *const phys);
const char *phys_return_correct_stype(Physics_T *const phys,
                               const char *const stype);
                              








