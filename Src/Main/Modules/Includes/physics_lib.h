#ifndef physics_LIB_H
#define physics_LIB_H

/* string length for string in Physics */
#define PAR_LEN   (99)
#define STEMP_LEN (99)

/* if the given patch DOES cover the region. 
// this is generally used in a loop over all patches */ 
#define IF_cover(patch,phys)     \
 if(IsItCovering(patch,phys->region))
 
/* if the given patch DOES NOT cover the region. 
// this is generally used in a loop over all patches */ 
#define IF_not_cover(patch,phys) \
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
 (sprintf(phys->par,"%s_%s",phys->ssys,param_name) ? Pgetd(phys->par) : DBL_MAX)

/* same as sysGetd but for integer type
// NOTE: not thread safe. */
#define sysGeti(param_name) \
 (sprintf(phys->par,"%s_%s",phys->ssys,param_name) ? Pgeti(phys->par) : INT_MAX)
 
/* same as sysGetd but for string type
// NOTE: not thread safe. */
#define sysGets(param_name) \
 (sprintf(phys->par,"%s_%s",phys->ssys,param_name) ? Pgets(phys->par) : NULL)


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

/* set double parameter for physical system. ex: 
// phys->ssys = "SBH",par = "Px_ADM" => 
// s = "SBH_Px_ADM"; thus, we can have various parameter
// calls with the same concept for different situations.
// NOTE: phys MUST be defined.
// NOTE: not thread safe. */
#define sysSetd(param_name,val) \
{sprintf(phys->par,"%s_%s",phys->ssys,param_name); Psetd(phys->par,(val));}

/* same as sysSetd but for integer.
// NOTE: not thread safe. */
#define sysSeti(param_name,val) \
{sprintf(phys->par,"%s_%s",phys->ssys,param_name); Pseti(phys->par,(val));}

/* same as sysSetd but for string.
// NOTE: not thread safe. */
#define sysSets(param_name,val) \
{sprintf(phys->par,"%s_%s",phys->ssys,param_name); Psets(phys->par,(val));}


/* handy string comparison against various values */
#define IF_sval(par,val) if(strcmp_i(Gets(par),val))

/* handy macro for "phy_return_correct_stype" function.
// NOTE: phys must be defined. */
#define Ftype(s) phys_autoindex_stype(phys,s)

/* physics function */
#define physics(phys,cmd) physics_main(phys,cmd,__FILE__,__LINE__)

/* commands, DON'T change the numeration. */
typedef enum CMD_T
{
 CMD_UNDEFINED = 0,
 
 /* star related */
 STAR_TUNE_EULER_CONST,
 STAR_TUNE_FORCE_BALANCE,
 STAR_TUNE_CENTER,
 STAR_EXTRAPOLATE_MATTERS,
 STAR_START,
 STAR_ADD_PARAMS,
 STAR_ADD_FIELDS,
 STAR_FIND_SURFACE,
 
 /* stress energy related */
 STRESS_ENERGY_UPDATE,
 STRESS_ENERGY_ADD_PARAMS,
 STRESS_ENERGY_ADD_FIELDS,
 
 /* BH related */
 BH_TUNE_RADIUS,
 BH_TUNE_SPIN,
 BH_FILL,
 BH_FIND_SURFACE,
 BH_START,
 BH_ADD_PARAMS,
 BH_ADD_FIELDS,
 
 /* free data related */
 FREE_DATA_ADD_PARAMS,
 FREE_DATA_ADD_FIELDS,
 FREE_DATA_POPULATE,
 
 /* system related */
 SYS_TUNE_P_ADM,
 
 CMD_END
}cmd_T;

/* various compact object types situation */
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
 NSNS,
 SNS,
 SBH
}Com_Obj_T;

/* forward declaration */
struct GRID_T;

/* struct for physics manager */
typedef struct PHYSICS_T
{
 struct GRID_T *grid;
 Grid_Char_T *grid_char;/* grid character when used for surface finder */
 unsigned igc;/* index of grid_char for this physics. */
 
 cmd_T cmd;/* current command */
 unsigned IsThisParent:1;/* if this is a parent physics 1, otherwise 0. */
 const char *region;/* grid region you want to issue the command */
 const char *Uregion;/* grid region specifed by User before to issue 
                     // the command, this is a backup, in case we need 
                     // to change region from original value. */
 
 Com_Obj_T ctype;/* handy for avoid many ifs; for instance when everything 
                 // is the same for NS, NS1 and NS2 one can check only 
                 // if (ctype == NS) ... . */
 Com_Obj_T type;/* BH1,NS2, NS, etc */
 const char *stype;/* string of type (above) used for parameter prefix */
 
 Com_Obj_T sys;/* system: BHNS, NSNS, BHBH, etc */
 const char *ssys;/* string type of system: "BHNS", "NSNS", etc */
 
 Flag_T pos;/* position of the object with respect to the whole grid
            // i.e. (LEFT|RIGHT|CENTER|NONE) */
 const char *spos;/* string type for position above, 
                  // used for patch collections */
 
 char par[PAR_LEN];/* related parameter for to be inquired. */
 /* some temp variables */
 char stemp[STEMP_LEN];/* a temperory string for various use. */
 struct GRID_T **gridtemp;/* temporaty grid for mygrid function */
 unsigned Ngridtemp;/* number of gridtemp */
 
}Physics_T;
#undef PAR_LEN
#undef STEMP_LEN

Physics_T *init_physics(Physics_T *const parent_phys,const Com_Obj_T type);
int physics_main(Physics_T *const phys,const cmd_T cmd,
            const char *const file, const int line);
void free_physics(Physics_T *phys);
void phys_set_region(Physics_T *const phys);
const char *phys_autoindex_stype(Physics_T *const phys,
                               const char *const stype);

struct GRID_T *mygrid(Physics_T *const phys,const char *const region);

#endif


