#ifndef physics_LIB_H
#define physics_LIB_H

/* string length for parameter in Physics */
#define PAR_LEN (99)

/* if the given patch DOES cover the region. 
// this is generally used in a loop over all patches */ 
#define if_cover(patch,phys)     \
 if(IsItCovering(patch,phys->region))
 
/* if the given patch DOES NOT cover the region. 
// this is generally used in a loop over all patches */ 
#define if_not_cover(patch,phys) \
 if(!IsItCovering(patch,phys->region))

/* get double parameter for the given physics. ex: 
// phys->stype = "NS1",par = "enthalpy_update_weight" => 
// s = "NS1_enthalpy_update_weight"; thus, we can have various parameter
// calls with the same concept for different compact objects.
// NOTE: phys MUST be defined. 
// NOTE: not thread safe. */
#define Getd(param_name) \
 (sprintf(phys->par,"%s_%s",phys->stype,param_name) ? Pgetd(phys->par) : DBL_MAX)

/* same as Getd but for integer type.
// NOTE: not thread safe. */
#define Geti(param_name) \
 (sprintf(phys->par,"%s_%s",phys->stype,param_name) ? Pgeti(phys->par) : INT_MAX)
 
/* same as Getd but for string type
// NOTE: not thread safe. */
#define Gets(param_name) \
 (sprintf(phys->par,"%s_%s",phys->stype,param_name) ? Pgets(phys->par) : NULL)


/* get double parameter for the given system. ex: 
// phys->ssys = "BHNS",par = "x_CM" => 
// s = "BHNS_x_CM"; thus, we can have various parameter
// calls with the same concept for different systems.
// NOTE: phys MUST be defined.
// NOTE: not thread safe. */
#define sysGetd(param_name) \
 (sprintf(phys->par,"%s_%s",phys->stype,param_name) ? Pgetd(phys->par) : DBL_MAX)

/* same as Getd but for integer type
// NOTE: not thread safe. */
#define sysGeti(param_name) \
 (sprintf(phys->par,"%s_%s",phys->stype,param_name) ? Pgeti(phys->par) : INT_MAX)
 
/* same as Getd but for string type
// NOTE: not thread safe. */
#define sysGets(param_name) \
 (sprintf(phys->par,"%s_%s",phys->stype,param_name) ? Pgets(phys->par) : NULL)


/* set double parameter for the given physics. ex: 
// phys->stype = "NS1",par = "enthalpy_update_weight" => 
// s = "NS1_enthalpy_update_weight"; thus, we can have various parameter
// calls with the same concept for different compact objects.
// NOTE: phys MUST be defined.
// NOTE: not thread safe. */
#define Setd(param_name,val) \
{sprintf(phys->par,"%s_%s",phys->stype,param_name); Psetd(phys->par,(val));}

/* same as Setd but for integer.
// NOTE: not thread safe. */
#define Seti(param_name,val) \
{sprintf(phys->par,"%s_%s",phys->stype,param_name); Pseti(phys->par,(val));}

/* same as Setd but for string.
// NOTE: not thread safe. */
#define Sets(param_name,val) \
{sprintf(phys->par,"%s_%s",phys->stype,param_name); Psets(phys->par,(val));}


/* commands, DON'T change the numeration. */
typedef enum CMD_T
{
 CMD_UNDEFINED = 0,
 STRESS_ENERGY,
 EULER_CONST,
 AH_RADIUS,
 P_ADM,
 AH_OMEGA,
 FORCE_BALANCE,
 FIX_CENTER,
 FIND_SURFACE,
 AH_NORMAL_VECTOR,
 EXTRAPOLATE_OUTSIDE,
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

/* struct for physics manager */
typedef struct PHYSICS_T
{
 Grid_T *grid;
 Grid_Char_T *grid_char;/* grid character when used for surface finder */
 cmd_T cmd;/* current command */
 const char *region;/* grid region you want to issue the command/ */
 
 Com_Obj_T type;/* BH1,NS2, NS, etc */
 const char *stype;/* string of type (above) used for parameter prefix */
 
 Com_Obj_T sys;/* system: BHNS, NSNS, BHBH, etc */
 const char *ssys;/* string type of system: "BHNS", "NSNS", etc */
 
 Flag_T pos;/* position of the object with respect to the whole grid
            // i.e. (LEFT|RIGHT|CENTER|NONE) */
 const char *spos;/* string type for position above, 
                  // used for patch collections */
 
 char par[PAR_LEN];/* related parameter for to be inquired. */
}Physics_T;

#undef PAR_LEN

Physics_T *init_physics(Grid_T *const grid,const Com_Obj_T type);
int physics(Physics_T *const phys,const cmd_T cmd,
            const char *const file, const int line);
void free_physics(Physics_T *phys);

#endif


