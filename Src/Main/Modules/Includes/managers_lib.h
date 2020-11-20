#ifndef managers_LIB_H
#define managers_LIB_H

/* size of an object parameter */
#define OPAR_SIZE (99)

/* if the given patch does NOT cover the region issue continue.
// this is generally used in a loop over all patches */
#define ONLY_IF_COVER(patch,obj) \
 if (!IsItCovering(patch,obj->region,obj->pos)) {continue;}

/* get double parameter for this given object. ex: 
// obj->stype = "NS1",par = "enthalpy_update_weight" => 
// s = "NS1_enthalpy_update_weight"; thus, we can have various parameter
// calls with the same concept for different compact objects. */
#define OPgetd(opar,obj,par) \
 (sprintf(opar,"%s_%s",obj->stype,par) ? Pgetd(opar) : DBL_MAX)

/* same as OPgetd but for integer type */ 
#define OPgeti(opar,obj,par) \
 (sprintf(opar,"%s_%s",obj->stype,par) ? Pgeti(opar) : INT_MAX)

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
 BHNS,
 BHBH,
 NSNS
}Com_Obj_T;

/* struct for object manager */
typedef struct OBJ_MAN_T
{
 Grid_T *grid;
 Grid_Char_T *grid_char;/* grid character when used for surface finder */
 cmd_T cmd;/* current command */
 const char *region;/* grid region you want to issue the command/ */
 Com_Obj_T type;/* BH1,NS2, NS, etc */
 const char *stype;/* string of type (above) used for parameter prefix */
 Com_Obj_T sys;/* system: BHNS, NSNS, BHBH, etc */
 Flag_T pos;/* position of the object with respect to the whole grid
            // i.e. (LEFT|RIGHT|CENTER|NONE) */
 const char *spos;/* string type for position above, 
                  // used for patch collections */
 
 
}Obj_Man_T;


Obj_Man_T *init_obj_man(Grid_T *const grid,const Com_Obj_T type);
void free_obj_man(Obj_Man_T *obj);

#endif


