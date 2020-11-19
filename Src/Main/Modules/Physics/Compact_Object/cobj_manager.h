#include "core_lib.h"
#include "utilities_lib.h"
#include "physics_compact_object_lib.h"


/* NOTE this MUST be as the same order as cmd_T */
static const char *cmd_enum_str[] =
{
 "CMD_UNDEFINED",
 "STRESS_ENERGY",
 "EULER_CONST",
 "AH_RADIUS",
 "P_ADM",
 "AH_OMEGA",
 "FORCE_BALANCE",
 "FIX_CENTER",
 "FIND_SURFACE",
 "AH_NORMAL_VECTOR",
 "EXTRAPOLATE_OUTSIDE",
 0
};


int NS_update(Compact_Obj_T *const obj,const cmd_T cmd);
int BH_update(Compact_Obj_T *const obj,const cmd_T cmd);
int system_update(Compact_Obj_T *const obj,const cmd_T cmd);


