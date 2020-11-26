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
  int ret = -1;
  char msg[STR_LEN] = {'\0'};
  
  phys->cmd   = cmd;
  
  /* if initialy user already specifed a region, make a backup */
  if (phys->region)
    phys->Uregion = phys->region;
    
  /* note: first must set phys->cmd */
  phys_set_region(phys);
  
  switch (cmd)
  {
    case TUNE_EULER_CONST:
      ret = star_main(phys);
    break;
    case TUNE_FORCE_BALANCE:
      ret = star_main(phys);
    break;
    case TUNE_NS_CENTER:
      ret = star_main(phys);
    break;
    case FIND_STAR_SURFACE:
      ret = star_main(phys);
    break;
    case UPDATE_STRESS_ENERGY:
      ret = Tij_main(phys);
    break;
    case EXTRAPOLATE_MATTERS:
      ret = star_main(phys);
    break;
   /* case AH_RADIUS:
      ret = update_apparent_horizon_radius(phys);
    break;
    case AH_OMEGA:
      ret = update_apparent_horizon_omega(phys);
    break;
    case AH_NORMAL_VECTOR:
      ret = update_apparent_horizon_normal(phys);
    break;*/
    /*case P_ADM:
      ret = adjust_ADM_momentum(phys);
    break;*/
    default:
      sprintf(msg,"No such command found!\n"
              "Incident triggered at\n"
              "file = %s\nline = %d",file,line);
      Error0(msg);
  }
  
  /* set to 0 to trap bugs */
  phys->cmd    = CMD_UNDEFINED;
  phys->region = 0;
  
  return ret;
}

/* call the requested function */
//int Mount(Physics_T *const phys,const cmd_T cmd)
//{
  //int ret = -1;
  
  //phys->cmd = cmd;
  
  //UNUSED(phys);
  //switch (cmd)
  //{
    //case UPDATE_STRESS_ENERGY:
      //ret = Tij_mount(phys);
    //break;
    /*case TUNE_EULER_CONST:
      ret = update_Euler_constant(phys);
    break;
    case FORCE_BALANCE:
      ret = adjust_force_balance_eq(phys);
    break;
    case FIX_CENTER:
      ret = fix_star_center(phys);
    break;
    case FIND_SURFACE:
      ret = find_star_surface(phys);
    break;
    case EXTRAPOLATE_OUTSIDE:
      ret = extrapolate_matter_outside_star(phys);
    break;*/
    /*case AH_RADIUS:
      ret = update_apparent_horizon_radius(phys);
    break;
    case AH_OMEGA:
      ret = update_apparent_horizon_omega(phys);
    break;
    case AH_NORMAL_VECTOR:
      ret = update_apparent_horizon_normal(phys);
    break;*/
    /*case P_ADM:
      ret = adjust_ADM_momentum(phys);
    break;*/
    //default:
      //Error0(NO_OPTION);
  //}
  
  /* set to no command to catch bug */
//  phys->cmd = CMD_UNDEFINED;
  
  //return ret;
//}

Physics_T *
init_physics
 (
 Grid_T *const grid/* computation grid */,
 const Com_Obj_T type/* object type NS,BH,etc */
 )
{
  Physics_T *phys = calloc(1,sizeof(*phys)); IsNull(phys);
  const char *spos  = 0;
  
  phys->grid = grid;
  phys->type = type;
  
  if (Pcmps("project","BH_NS_initial_data"))
  {
    phys->sys  = BHNS;
    phys->ssys = "BHNS";
  }
  else if (Pcmps("project","NS_NS_initial_data"))
  {
    phys->sys  = NSNS;
    phys->ssys = "NSNS";
  }
  else if (Pcmps("project","BH_BH_initial_data"))
  {
    phys->sys  = BHBH;
    phys->ssys = "BHBH";
  }
  else if (Pcmps("project","BH_initial_data"))
  {
    phys->sys  = SBH;
    phys->ssys = "SBH";/* important to have different name for system */
  }
  else if (Pcmps("project","NS_initial_data"))
  {
    phys->sys  = SNS;
    phys->ssys = "SNS";/* important to have different name for system */
  }
  else
    Error0(NO_OPTION);
  
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
        Error0(NO_OPTION);
        
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
        Error0(NO_OPTION);
      
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
        Error0(NO_OPTION);
      
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
        Error0(NO_OPTION);
        
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
        Error0(NO_OPTION);
      
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
        Error0(NO_OPTION);
      
    break;
    case BHNS:
      phys->stype = "BHNS";
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
  unsigned i;
  
  if (phys)
  {
    for (i = 0; i < phys->Ngridtemp; ++i)
    {
      _free(phys->gridtemp[i]->patch);
      _free(phys->gridtemp[i]);
    }
  _free(phys->gridtemp);
  }
  _free(phys);
}

/* set default phys->region  */
void phys_set_region(Physics_T *const phys)
{
  /* set natural region, later one can make it more sophisticated */
  if (phys->Uregion)
    phys->region = phys->Uregion;
  else
    phys->region = phys->stype;
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
  /* some checks: */
  /* if there is no NS or BH */
  if (!strstr(stype,"NS") && !strstr(stype,"BH"))
      Error1("Argument (%s) is not supported!",stype);
  /* if two different object asked, right now this is not supported! */
  if (strstr(stype,"NS") && strstr(stype,"BH"))
      Error1("Two different objects (%s) are not supported!",stype);
  if (strstr(stype,"BH1") && strstr(stype,"BH2"))
      Error1("Two different objects (%s) are not supported!",stype);
  if (strstr(stype,"NS1") && strstr(stype,"NS2"))
      Error1("Two different objects (%s) are not supported!",stype);
  
  /* if indices 1 or 2 is already taken into account 
  // or the object doesn't have any index, don't change anything. */
  if (phys->type == NS  || phys->type == BH ||
      strchr(stype,'1') || strchr(stype,'2'))
    return stype;
  
  /* having made sure everything is find now do a simple autoindex
  // it replace NS (BH) with NSi (BHi) in which i is the correct index. */
  regex_replace(stype,"(NS|BH)",phys->stype,phys->stemp);
  
  return phys->stemp;
}


/* a handy function to gather pertinent patches to the given region
// into a grid. 
// NOTE: no need to free, it is freed when free_physics is called. 
// NOTE: at each call it adds one grid to phys->gridtemp. */
Grid_T *mygrid(Physics_T *const phys,const char *const region)
{
  const unsigned ng = phys->Ngridtemp;
  unsigned Np;
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
  
  
  return phys->gridtemp[ng];
}

