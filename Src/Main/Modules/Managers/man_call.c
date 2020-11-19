/*
// Alireza Rashti
// November 2020
*/

/* call the project pertinent to system update and adjustments */

#include "man_call.h"

/* call the requested function */
int Update(Obj_Man_T *const obj,const cmd_T cmd)
{
  int ret = -1;
  
  obj->cmd = cmd;
  
  UNUSED(obj);
  switch (cmd)
  {
    case STRESS_ENERGY:
      ret = Tij_update(obj);
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
    default:
      Error0(NO_OPTION);
  }
  
  /* set to no command to catch bug */
  obj->cmd = CMD_UNDEFINED;
  
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

