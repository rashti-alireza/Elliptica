/*
// Alireza Rashti
// November 2020
*/

/* call the project pertinent to system update and adjustments */

#include "man_call.h"

/* call the requested function */
int man_tune_call(Obj_Man_T *const obj,const cmd_T cmd,
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
    //case EULER_CONST:
      //ret = star_tune(obj);
    //break;
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

