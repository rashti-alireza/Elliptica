/*
// Alireza Rashti
// November 2020
*/

/* call the project pertinent to system update and adjustments among
// others. this is supposed to be a physics manager of various physics
// projects and concepts. this helps to stay on track for further 
// developments and easily using of physics between the projects. */

#include "phys_main.h"

/* call the requested function */
int physics_main(Physics_T *const phys,const cmd_T cmd,
            const char *const file, const int line)
{
  int ret = EXIT_SUCCESS;
  char msg[STR_LEN] = {'\0'};
  
  phys->cmd = cmd;
  
  switch (cmd)
  {
    case STAR_TUNE_EULER_CONST:
    case STAR_TUNE_FORCE_BALANCE:
    case STAR_TUNE_CENTER:
    case STAR_FIND_SURFACE:
    case STAR_START:
    case STAR_ADD_FIELDS:
    case STAR_SET_PARAMS:
    case STAR_EXTRAPOLATE_MATTERS:
      ret = star_main(phys);
    break;
    
    case STRESS_ENERGY_UPDATE:
    case STRESS_ENERGY_SET_PARAMS:
    case STRESS_ENERGY_ADD_FIELDS:
      ret = Tij_main(phys);
    break;
    
    case BH_FIND_SURFACE:
    case BH_TUNE_RADIUS:
    case BH_TUNE_SPIN:
    case BH_FILL:
    case BH_START:
    case BH_SET_PARAMS:
    case BH_ADD_FIELDS:
    case BH_UPDATE_sConf:
    case BH_UPDATE_INNER_BC:
      ret = bh_main(phys);
    break;
    
    case FREE_DATA_SET_PARAMS:
    case FREE_DATA_ADD_FIELDS:
    case FREE_DATA_POPULATE:
      ret = fd_main(phys);
    break;
    
    case SYS_SET_PARAMS:
    case SYS_ADD_FIELDS:
    case SYS_TUNE_P_ADM:
    case SYS_INITIALIZE_FIELDS:
      ret = sys_main(phys);
    break;
    
    case ADM_SET_PARAMS:
    case ADM_ADD_FIELDS:
    case ADM_UPDATE_Kij:
    case ADM_UPDATE_KIJ:
    case ADM_UPDATE_gij:
    case ADM_UPDATE_B1I:
    case ADM_UPDATE_beta:
    case ADM_UPDATE_AConfIJ:
    case ADM_COMPUTE_CONSTRAINTS:
    case ADM_DOCTEST:
      ret = adm_main(phys);
    break;
    
    case OBSERVE_SET_PARAMS:
    case OBSERVE_ADD_FIELDS:
      ret = observe_main(phys);
    break;
    
    case EQ_SET_PARAMS:
    case EQ_ADD_FIELDS:
    case EQ_SOLVE:
      ret = eq_main(phys);
    break;
    
    default:
      sprintf(msg,"No such command found!\n"
              "Incident triggered at:\n"
              "file = %s\nline = %d",file,line);
      Error0(msg);
  }
  
  /* set to 0 to trap bugs */
  phys->cmd = CMD_UNDEFINED;
  
  return ret;
}

/* initialize a phys from a parent_phys(if any).
// note: for parant physics, one must set phys->grid manually 
// after initialization. */
Physics_T *
init_physics
 (
 Physics_T *const parent_phys/* if null, it means this is a parent physics */,
 const Com_Obj_T type/* object type NS,BH,etc */
 )
{
  Physics_T *phys = calloc(1,sizeof(*phys)); IsNull(phys);
  const char *spos  = 0;
  
  /* if given parent_phys in not null => we have a child physics */
  if (parent_phys)
  {
    phys->grid         = parent_phys->grid;
    phys->IsThisParent = 0;
  }
  else
  {
    phys->IsThisParent = 1;
  }
  
  phys->type = type;
  
  if (Pcmps("project","BH_NS_binary_initial_data"))
  {
    phys->sys  = BHNS;
    phys->ssys = "BHNS";
  }
  else if (Pcmps("project","NS_NS_binary_initial_data"))
  {
    phys->sys  = NSNS;
    phys->ssys = "NSNS";
  }
  else if (Pcmps("project","BH_BH_binary_initial_data"))
  {
    phys->sys  = BHBH;
    phys->ssys = "BHBH";
  }
  else if (Pcmps("project","Single_BH_initial_data"))
  {
    phys->sys  = SBH;
    phys->ssys = "SBH";/* important to have different name for system */
  }
  else if (Pcmps("project","Sinlge_NS_initial_data"))
  {
    phys->sys  = SNS;
    phys->ssys = "SNS";/* important to have different name for system */
  }
  else if (Pcmps("project","Modules_Test"))
  {
    phys->sys  = OBJ_UNDEFINED;
    phys->ssys = "OBJ_UNDEFINED";/* important to have different name for system */
  }
  else if (Pcmps("project","TOV_Star"))
  {
    phys->sys  = OBJ_UNDEFINED;
    phys->ssys = "OBJ_UNDEFINED";/* important to have different name for system */
  }
  else
  {
    Error0(NO_OPTION);
  }
  
  switch(type)
  {
    case NS:
      phys->ctype = NS;
      phys->stype = "NS";
      spos = Pgets("grid_set_NS");
      if (strstr_i(spos,"left"))
      {
        phys->pos  = LEFT;
        phys->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        phys->pos  = RIGHT;
        phys->spos = "right";
      }
      else if (strstr_i(spos,"center"))
      {
        phys->pos  = CENTER;
        phys->spos = "center";
      }
      else
      {
        Error0(NO_OPTION);
      }
    break;
    
    case NS1:
      phys->ctype = NS;
      phys->stype = "NS1";
      spos = Pgets("grid_set_NS1");
      if (strstr_i(spos,"left"))
      {
        phys->pos  = LEFT;
        phys->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        phys->pos  = RIGHT;
        phys->spos = "right";
      }
      else
      {
        Error0(NO_OPTION);
      }
    break;
    
    case NS2:
      phys->ctype = NS;
      phys->stype = "NS2";
      spos = Pgets("grid_set_NS2");
      if (strstr_i(spos,"left"))
      {
        phys->pos  = LEFT;
        phys->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        phys->pos  = RIGHT;
        phys->spos = "right";
      }
      else
      {
        Error0(NO_OPTION);
      }
    break;
    
    case BH:
      phys->ctype = BH;
      phys->stype = "BH";
      spos = Pgets("grid_set_BH");
      if (strstr_i(spos,"left"))
      {
        phys->pos  = LEFT;
        phys->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        phys->pos  = RIGHT;
        phys->spos = "right";
      }
      else if (strstr_i(spos,"center"))
      {
        phys->pos  = CENTER;
        phys->spos = "center";
      }
      else
      {
        Error0(NO_OPTION);
      }
    break;
    
    case BH1:
      phys->ctype = BH;
      phys->stype = "BH1";
      spos = Pgets("grid_set_BH1");
      if (strstr_i(spos,"left"))
      {
        phys->pos  = LEFT;
        phys->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        phys->pos  = RIGHT;
        phys->spos = "right";
      }
      else
      {
        Error0(NO_OPTION);
      }
    break;
    
    case BH2:
      phys->ctype = BH;
      phys->stype = "BH2";
      spos = Pgets("grid_set_BH2");
      if (strstr_i(spos,"left"))
      {
        phys->pos  = LEFT;
        phys->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        phys->pos  = RIGHT;
        phys->spos = "right";
      }
      else
      {
        Error0(NO_OPTION);
      }
    break;

    case SBH:
      phys->ctype = SBH;
      phys->stype = "SBH";
      spos = Pgets("grid_set_BH");
      if (strstr_i(spos,"center"))
      {
        phys->pos  = CENTER;
        phys->spos = "center";
      }
      else
      {
        Error0(NO_OPTION);
      }
    break;
    
    case BHNS:
      phys->ctype = BHNS;
      phys->stype = "BHNS";
      phys->pos   = NONE;
    break;
    
    case NSNS:
      phys->ctype = NSNS;
      phys->stype = "NSNS";
      phys->pos   = NONE;
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return phys;
}

/* free  */
void free_physics(Physics_T *phys)
{
  Uint i;
  
  if (phys)
  {
    /* free temp grid */
    for (i = 0; i < phys->Ngridtemp; ++i)
    {
      Free(phys->gridtemp[i]->patch);
      Free(phys->gridtemp[i]);
    }
    Free(phys->gridtemp);
  
    /* free grid and its parameters */
    if (phys->IsThisParent)
    {
      free_grid_params(phys->grid);
      free_grid(phys->grid);
    }
    
    free(phys);
  }
}


/* ->: stype with correct indices (kind of auto spell).
// since we have different indices for objects like NS1, BH2 etc
// and some function are very similar but the stype indices are different
// this function gets a prototype and with respect to the given
// physics adjust the correct stype. eg:
// if phys->type = NS2 => 
// phys_autoindex_stype(phys,"NS,NS_around") = "NS2,NS2_around". 
// NOTE: it's not very limited, please see the if's in the function.
// NOTE: not thread safe */
const char *phys_autoindex_stype(Physics_T *const phys,
                               const char *const stype)
{
  /* if matches all */
  if (!strcmp(stype,".*"))
    return ".*";
    
  /* if there is no NS or BH for instance it is outermost */
  if (!strstr(stype,"NS") && !strstr(stype,"BH"))
    return stype;
  
  /* if two different object asked */
  if (strstr(stype,"NS") && strstr(stype,"BH"))
  {
    if (phys->grid->kind == Grid_SplitCubedSpherical_BHNS)
      return stype;
    else
      Errors("%s is not supported!\n",stype);
  }
  
  if (strstr(stype,"BH1") && strstr(stype,"BH2"))
  {
    if (phys->grid->kind == Grid_SplitCubedSpherical_BHBH)
      return stype;
    else
      Errors("%s is not supported!\n",stype);
  }
  
  if (strstr(stype,"NS1") && strstr(stype,"NS2"))
  {
    if (phys->grid->kind == Grid_SplitCubedSpherical_NSNS)
      return stype;
    else
      Errors("%s is not supported!\n",stype);
  }
  
  /* if indices 1 or 2 is already taken into account 
  // or the object doesn't have any index, don't change anything. */
  if (phys->type == NS  || phys->type == BH ||
      strchr(stype,'1') || strchr(stype,'2'))
    return stype;
  
  /* having assured everything is fine, we now do a simple autoindex.
  // it replaces NS (BH) with NSi (BHi) in which i is the correct index. */
  AssureType(phys->ctype == NS || phys->ctype == BH);
  regex_replace(stype,"(NS|BH)",phys->stype,phys->stemp);
  
  return phys->stemp;
}


/* a handy function to gather pertinent patches to the given region
// into a grid. 
// NOTE: no need to free, it is freed when free_physics is called. 
// NOTE: at each call it adds one grid to phys->gridtemp. */
Grid_T *mygrid(Physics_T *const phys,const char *const region)
{
  const Uint ng = phys->Ngridtemp;
  Uint Np;
  Patch_T **patches = collect_patches(phys->grid,Ftype(region),&Np);
  
  phys->gridtemp = realloc(phys->gridtemp,(ng+2)*sizeof(*phys->gridtemp));
  IsNull(phys->gridtemp);
  
  phys->Ngridtemp      = ng+1;
  phys->gridtemp[ng+1] = 0;
  phys->gridtemp[ng]   = calloc(1,sizeof(*phys->gridtemp[ng]));
  IsNull(phys->gridtemp[ng]);
  
  phys->gridtemp[ng]->kind  = phys->grid->kind;
  phys->gridtemp[ng]->patch = patches;
  phys->gridtemp[ng]->np    = Np;
  phys->gridtemp[ng]->gn    = phys->grid->gn;
  
  return phys->gridtemp[ng];
}

