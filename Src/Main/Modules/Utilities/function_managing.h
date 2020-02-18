#include "core_lib.h"
#include "error_handling_lib.h"


void init_func_PtoV(sFunc_PtoV_T ***const func);
void add_func_PtoV(sFunc_PtoV_T ***const func,void (*f)(Patch_T *const patch),const char *const task,const Coord_T coord);
void run_func_PtoV(sFunc_PtoV_T **const func,const char *const task,Patch_T *const patch);
void init_func_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const func);
void add_func_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const func,double *(*f)(Patch_T *const patch),const char *const name);
sFunc_Patch2Pdouble_T *get_func_Patch2Pdouble(const char *const name,sFunc_Patch2Pdouble_T **const func);
void *alloc_sFunc_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const mem);
void free_func_PtoV(sFunc_PtoV_T **func);
void *alloc_sFunc_PtoV(sFunc_PtoV_T ***const mem);


