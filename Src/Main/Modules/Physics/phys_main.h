#include "core_lib.h"
#include "utilities_lib.h"
#include "physics_lib.h"
#include "physics_StressEnergyTensor_lib.h"
#include "physics_star_lib.h"

#define STR_LEN (999)

int physics(Obj_Man_T *const obj,const cmd_T cmd,
            const char *const file, const int line);

Obj_Man_T *
init_obj_man
 (
 Grid_T *const grid/* computation grid */,
 const Com_Obj_T type/* object type NS,BH,etc */
 );

void free_obj_man(Obj_Man_T *obj);



