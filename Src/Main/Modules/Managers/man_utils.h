#include "core_lib.h"
#include "utilities_lib.h"
#include "managers_lib.h"

Obj_Man_T *
init_obj_man
 (
 Grid_T *const grid/* computation grid */,
 const Com_Obj_T type/* object type NS,BH,etc */
 );

void free_obj_man(Obj_Man_T *obj);


