#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"


void init_func_PtoV(sFunc_PtoV_T ***const func);
void add_func_PtoV(sFunc_PtoV_T ***const func,void (*f)(Patch_T *const patch),const char *const task,const Coord_T coord);
void run_func_PtoV(sFunc_PtoV_T **const func,const char *const task,Patch_T *const patch);
Coord_T find_coord(const Patch_T *const patch);

