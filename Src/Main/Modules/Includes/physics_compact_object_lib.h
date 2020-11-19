#ifndef physics_COMPACT_OBJ_LIB_H
#define physics_COMPACT_OBJ_LIB_H

/* commands, DON'T change the numeration and ADD after one to the last */
typedef enum CMD_T
{
 CMD_UNDEFINED = 0,
 STRESS_ENERGY = 1,
 EULER_CONST   = 2,
 AH_RADIUS     = 3,
 P_ADM         = 4,
 AH_OMEGA      = 5,
 FORCE_BALANCE = 6,
 FIX_CENTER    = 7,
 FIND_SURFACE  = 8,
 AH_NORMAL_VECTOR    = 9,
 EXTRAPOLATE_OUTSIDE = 10,
 CMD_END
}cmd_T;

/* compact object type */
typedef enum COMP_OBJ_T
{
 OBJ_UNDEFINED = 0,
 NS,
 NS1,
 NS2,
 BH,
 BH1,
 BH2,
 BHNS
 BHBH,
 NSNS,
}Com_Obj_T;

/* struct for compact object */
typedef struct COMPACT_OBJ_T
{
 Grid_T *grid;
 Grid_Char_T *grid_char;/* grid character when used for surface finder */
 Com_Obj_T type;/* BH1,NS2, NS, etc */
 const char *stype;/* string of type (above) */
 Com_Obj_T sys;/* system: BHNS, NSNS, BHBH, etc */
 Flag_T pos;/* position of the object with respect to the whole grid
            // i.e. (LEFT|RIGHT|CENTER|NONE) */
 const char *spos;/* string type for position above, 
                  // used for patch collections */
 
 
}Compact_Obj_T;

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

Compact_Obj_T *init_compact_obj(Grid_T *const grid,const char *const sobj);
void free_compact_obj(Grid_T *const grid *obj);

#endif


