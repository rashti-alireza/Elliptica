#include "core_lib.h"
#include "utilities_lib.h"
#include "physics_compact_object_lib.h"

Compact_Obj_T *
init_compact_obj
 (
 Grid_T *const grid/* computation grid */,
 const Com_Obj_T type/* object type NS,BH,etc */
 );

void free_compact_obj(Compact_Obj_T *obj);


