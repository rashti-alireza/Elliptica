/*
// Alireza Rashti
// November 2020
*/

/* call the project pertinent to system update and adjustments */

#include "phys_main.h"

/* call the requested function */
int physics(Obj_Man_T *const obj,const cmd_T cmd,
            const char *const file, const int line)
{
  int ret = -1;
  char msg[STR_LEN] = {'\0'};
  
  obj->cmd = cmd;
  
  if(!obj->region)
  {
    Error0("NO region specified!\n");
  }
  
  switch (cmd)
  {
    case STRESS_ENERGY:
      ret = Tij_tune(obj);
    break;
    case EULER_CONST:
      ret = star_tune(obj);
    break;
    /*case FORCE_BALANCE:
      ret = star_tunes(obj);
    break;
    case FIX_CENTER:
      ret = star_tunes(obj);
    break;
    case FIND_SURFACE:
      ret = star_tunes(obj);
    break;
    case EXTRAPOLATE_OUTSIDE:
      ret = star_tunes(obj);
    break;
    case AH_RADIUS:
      ret = update_apparent_horizon_radius(obj);
    break;
    case AH_OMEGA:
      ret = update_apparent_horizon_omega(obj);
    break;
    case AH_NORMAL_VECTOR:
      ret = update_apparent_horizon_normal(obj);
    break;*/
    /*case P_ADM:
      ret = adjust_ADM_momentum(obj);
    break;*/
    default:
      sprintf(msg,"No such command found!\n"
              "Incident triggered at\n"
              "file = %s\nline = %d",file,line);
      Error0(msg);
  }
  
  /* set to 0 to catch bug and miss you */
  obj->cmd    = CMD_UNDEFINED;
  obj->region = 0;
  
  return ret;
}

/* call the requested function */
//int Mount(Obj_Man_T *const obj,const cmd_T cmd)
//{
  //int ret = -1;
  
  //obj->cmd = cmd;
  
  //UNUSED(obj);
  //switch (cmd)
  //{
    //case STRESS_ENERGY:
      //ret = Tij_mount(obj);
    //break;
    /*case EULER_CONST:
      ret = update_Euler_constant(obj);
    break;
    case FORCE_BALANCE:
      ret = adjust_force_balance_eq(obj);
    break;
    case FIX_CENTER:
      ret = fix_star_center(obj);
    break;
    case FIND_SURFACE:
      ret = find_star_surface(obj);
    break;
    case EXTRAPOLATE_OUTSIDE:
      ret = extrapolate_matter_outside_star(obj);
    break;*/
    /*case AH_RADIUS:
      ret = update_apparent_horizon_radius(obj);
    break;
    case AH_OMEGA:
      ret = update_apparent_horizon_omega(obj);
    break;
    case AH_NORMAL_VECTOR:
      ret = update_apparent_horizon_normal(obj);
    break;*/
    /*case P_ADM:
      ret = adjust_ADM_momentum(obj);
    break;*/
    //default:
      //Error0(NO_OPTION);
  //}
  
  /* set to no command to catch bug */
//  obj->cmd = CMD_UNDEFINED;
  
  //return ret;
//}

Obj_Man_T *
init_obj_man
 (
 Grid_T *const grid/* computation grid */,
 const Com_Obj_T type/* object type NS,BH,etc */
 )
{
  Obj_Man_T *obj = calloc(1,sizeof(*obj)); IsNull(obj);
  const char *spos  = 0;
  
  Error0("/* what should i do for BHNS region? */");
  
  assert(obj->region);
  
  obj->grid = grid;
  obj->type = type;
  
  if (Pcmps("project","BH_NS_initial_data"))
  {
    obj->sys  = BHNS;
    obj->ssys = "BHNS";
    
  }
  else if (Pcmps("project","NS_NS_initial_data"))
  {
    obj->sys  = NSNS;
    obj->ssys = "NSNS";
    
  }
  else if (Pcmps("project","BH_BH_initial_data"))
  {
    obj->sys  = BHBH;
    obj->ssys = "BHBH";
    
  }
  else
    Error0(NO_OPTION);
  
  switch(type)
  {
    case NS:
      obj->stype = "NS";
      spos = Pgets("grid_set_NS");
      if (strstr_i(spos,"left"))
      {
        obj->pos  = LEFT;
        obj->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        obj->pos  = RIGHT;
        obj->spos = "right";
      }
      else if (strstr_i(spos,"center"))
      {
        obj->pos  = CENTER;
        obj->spos = "center";
      }
      else
        Error0(NO_OPTION);
        
    break;
    case NS1:
      obj->stype = "NS1";
      spos = Pgets("grid_set_NS1");
      if (strstr_i(spos,"left"))
      {
        obj->pos  = LEFT;
        obj->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        obj->pos  = RIGHT;
        obj->spos = "right";
      }
      else
        Error0(NO_OPTION);
      
    break;
    case NS2:
      obj->stype = "NS2";
      spos = Pgets("grid_set_NS2");
      if (strstr_i(spos,"left"))
      {
        obj->pos  = LEFT;
        obj->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        obj->pos  = RIGHT;
        obj->spos = "right";
      }
      else
        Error0(NO_OPTION);
      
    break;
    case BH:
      obj->stype = "BH";
      spos = Pgets("grid_set_BH");
      if (strstr_i(spos,"left"))
      {
        obj->pos  = LEFT;
        obj->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        obj->pos  = RIGHT;
        obj->spos = "right";
      }
      else if (strstr_i(spos,"center"))
      {
        obj->pos  = CENTER;
        obj->spos = "center";
      }
      else
        Error0(NO_OPTION);
        
    break;
    case BH1:
      obj->stype = "BH1";
      spos = Pgets("grid_set_BH1");
      if (strstr_i(spos,"left"))
      {
        obj->pos  = LEFT;
        obj->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        obj->pos  = RIGHT;
        obj->spos = "right";
      }
      else
        Error0(NO_OPTION);
      
    break;
    case BH2:
      obj->stype = "BH2";
      spos = Pgets("grid_set_BH2");
      if (strstr_i(spos,"left"))
      {
        obj->pos  = LEFT;
        obj->spos = "left";
      }
      else if (strstr_i(spos,"right"))
      {
        obj->pos  = RIGHT;
        obj->spos = "right";
      }
      else
        Error0(NO_OPTION);
      
    break;
    case BHNS:
      obj->stype = "BHNS";
      obj->pos   = NONE;
      
    break;
    default:
      Error0(NO_OPTION);
  }
  
  return obj;
}

/* free  */
void free_obj_man(Obj_Man_T *obj)
{
  _free(obj);
}

