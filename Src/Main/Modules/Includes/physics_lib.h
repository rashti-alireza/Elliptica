#ifndef physics_LIB_H
#define physics_LIB_H
#include "elliptica_system_lib.h"

/* string length for string in Physics */
#define PAR_LEN   (99)
#define STEMP_LEN (99)

 
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

/* property print double */
#define PR_PROPERTY_IN_FILE_d(x,file,pr_screen) \
  { \
    const double par___temp = Getd(x); \
    fprintf(file,"%-30s = %+0.15f\n",phys->par,par___temp);\
    if (pr_screen) printf("%-30s = %+0.15f\n",phys->par,par___temp); \
  }

/* property print string */
#define PR_PROPERTY_IN_FILE_s(x,file,pr_screen) \
   {\
     const char *par___temp = Gets(x); \
     fprintf(file,"%-30s = %s\n",phys->par,par___temp);\
     if (pr_screen) printf("%-30s = %s\n",phys->par,par___temp); \
   }


/* rest assure the criterion is met for physics type
// usage ex: AssureType(phys->ctype == NS). */
#define AssureType(x) \
{if(!(x))Errors("Wrong physics type passed to function '%s'.\n",__func__);}

/* commands, DON'T change the numeration. */
typedef enum CMD_T
{
 CMD_UNDEFINED = 0,
 
 /* star related */
 STAR_TUNE_EULER_CONST,/* Euler const. in fluid eq. */
 STAR_TUNE_FORCE_BALANCE,/* dh/d? at star center */
 STAR_TUNE_CENTER,/* avoid star drifting */
 STAR_EXTRAPOLATE_MATTERS,/* extrapolate outside star */
 STAR_START,/* initiation before first grid */
 STAR_SET_PARAMS,/* set params */
 STAR_ADD_FIELDS,/* add fields */
 STAR_FIND_SURFACE,/* find new star surface */
 
 /* stress energy related */
 STRESS_ENERGY_UPDATE,/* update whole matter related fields */
 STRESS_ENERGY_SET_PARAMS,/* set params */
 STRESS_ENERGY_ADD_FIELDS,/* add fields */
 
 /* BH related */
 BH_TUNE_RADIUS,/* adjust BH radius */
 BH_TUNE_SPIN,/* adjust BH spin */
 BH_FILL,/* BH-filler */
 BH_FIND_SURFACE,/* find surface of BH to create grid */
 BH_START,/* initiation before first grid */
 BH_SET_PARAMS,/* set params */
 BH_ADD_FIELDS,/* add fields */
 BH_UPDATE_sConf,/* update sConf^i and dsConf^i_j conformal normal on AH */
 BH_UPDATE_INNER_BC,/* update inner BC values on AH */
 
 /* free data related */
 FREE_DATA_SET_PARAMS,/* set params */
 FREE_DATA_ADD_FIELDS,/* add fields */
 FREE_DATA_POPULATE,/* populate free data */
 
 /* system related */
 SYS_SET_PARAMS,/* set params */
 SYS_ADD_FIELDS,/* add fields */
 SYS_TUNE_P_ADM,/* adjust P_adm of the system */
 SYS_INITIALIZE_FIELDS,/* start off with these fields */
 
 /* adm related */
 ADM_SET_PARAMS,/* set params */
 ADM_ADD_FIELDS,/* add fields */
 ADM_UPDATE_Kij,/* adm_K_{ij} */
 ADM_UPDATE_KIJ,/* adm_K^{ij} */
 ADM_UPDATE_gij,/* adm_g_{ij} = psi^4 *gConf_{ij} */
 ADM_UPDATE_AConfIJ,/* conformal traceless part of K^{ij} */
 ADM_UPDATE_beta,/* shift */
 ADM_UPDATE_B1I,/* beta = B0+B1 in which B1 can be inspiral piece */
 ADM_COMPUTE_CONSTRAINTS,/* ham and mom constrains */
 ADM_DOCTEST,/* some internal tests */
 
 /* observe related */
 OBSERVE_SET_PARAMS,/* set params */
 OBSERVE_ADD_FIELDS,/* add fields */
 
 /* equation related */
 EQ_SET_PARAMS,/* set params */
 EQ_ADD_FIELDS,/* add fields */
 EQ_SOLVE,/* solve eq(s) */
 
 
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
 Uint igc;/* index of grid_char for this physics. */
 
 cmd_T cmd;/* current command */
 Uint IsThisParent:1;/* if this is a parent physics 1, otherwise 0. */
 
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
 Uint Ngridtemp;/* number of gridtemp */
 
}Physics_T;
#undef PAR_LEN
#undef STEMP_LEN

Physics_T *init_physics(Physics_T *const parent_phys,const Com_Obj_T type);
int physics_main(Physics_T *const phys,const cmd_T cmd,
            const char *const file, const int line);
void free_physics(Physics_T *phys);
const char *phys_autoindex_stype(Physics_T *const phys,
                               const char *const stype);

struct GRID_T *mygrid(Physics_T *const phys,const char *const region);

#endif


