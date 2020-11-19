/*
// Alireza Rashti
// November 2020
*/

/* call the algorithms pertinent to system update and adjustments */

#include "manager.h"

/* call the requested function */
int NS_update(Compact_Obj_T *const obj,const cmd_T cmd)
{
  int ret = -1;
  
  switch (cmd)
  {
    case STRESS_ENERGY:
      ret = update_stress_energy_tensor(obj);
    break;
    case EULER_CONST:
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
    break;
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}


/* call the requested function */
int BH_update(Compact_Obj_T *const obj,const cmd_T cmd)
{
  int ret = -1;
  
  switch (cmd)
  {
    case AH_RADIUS:
      ret = update_apparent_horizon_radius(obj);
    break;
    case AH_OMEGA:
      ret = update_apparent_horizon_omega(obj);
    break;
    case AH_NORMAL_VECTOR:
      ret = update_apparent_horizon_normal(obj);
    break;
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* call the requested function */
int system_update(Compact_Obj_T *const obj,const cmd_T cmd)
{
  int ret = -1;
  
  switch (cmd)
  {
    case P_ADM:
      ret = adjust_ADM_momentum(obj);
    break;
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

