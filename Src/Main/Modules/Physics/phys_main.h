#include "core_lib.h"
#include "utilities_lib.h"
#include "physics_lib.h"
#include "physics_stress_energy_lib.h"
#include "physics_star_lib.h"
#include "physics_blackhole_lib.h"
#include "physics_freedata_lib.h"
#include "physics_system_lib.h"
#include "physics_adm_lib.h"
#include "physics_observe_lib.h"
#include "physics_equation_lib.h"

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
const char *phys_autoindex_stype(Physics_T *const phys,
                               const char *const stype);

Grid_T *mygrid(Physics_T *const phys,const char *const region);










