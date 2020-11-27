#include "core_lib.h"
#include "utilities_lib.h"
#include "physics_lib.h"
#include "physics_StressEnergyTensor_lib.h"
#include "physics_star_lib.h"

#define STR_LEN (999)

int physics_main(Physics_T *const phys,const cmd_T cmd,
            const char *const file, const int line);

Physics_T *
init_physics
 (
 Physics_T *const parent_phys/* if null, it means this is a parent physics */,
 const Com_Obj_T type/* object type NS,BH,etc */
 );

void free_physics(Physics_T *obj);
void phys_set_region(Physics_T *const phys);
const char *phys_autoindex_stype(Physics_T *const phys,
                               const char *const stype);

Grid_T *mygrid(Physics_T *const phys,const char *const region);










