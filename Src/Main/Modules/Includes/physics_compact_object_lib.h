#ifndef physics_COMPACT_OBJ_LIB_H
#define physics_COMPACT_OBJ_LIB_H

/* commands */
typedef enum CMD_T
{
 CMD_UNDEFINED = -1,
 STRESS_ENERGY = 0,
 EULER_CONST   = 1,
 AH_RADIUS     = 2,
 P_ADM         = 3,
 AH_OMEGA      = 4,
 FORCE_BALANCE = 5,
 FIX_CENTER    = 6,
 FIND_SURFACE  = 7,
 AH_NORMAL_VECTOR    = 8,
 EXTRAPOLATE_OUTSIDE = 9,
 CMD_END
}cmd_T;

/* struct for compact object */
typedef struct COMPACT_OBJ_T
{
 Grid_T *grid;
 Grid_Char_T *grid_char;/* grid character when used for surface finder */
 Flag_T pos;/* position of the object (LEFT|RIGHT|CENTER|NONE) */
 
 
}Compact_Obj_T;

/* NOTE this MUST be as the same order as cmd_T */
static const char *cmd_enum_str[] =
{
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

Compact_Obj_T *init_compact_obj(Grid_T *const grid,const char *const obj);
void free_compact_obj(Grid_T *const grid *obj);

#endif


