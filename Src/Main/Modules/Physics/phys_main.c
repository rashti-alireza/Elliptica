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
int physics(Physics_T *const phys,const cmd_T cmd,
            const char *const file, const int line)
{
  int ret = -1;
  char msg[STR_LEN] = {'\0'};
  
  phys->cmd   = cmd;
  
  /* note: first must set phys->cmd */
  set_phys_default_region(phys);
  
  switch (cmd)
  {
    case STRESS_ENERGY:
      ret = Tij_tune(phys);
    break;
    case EULER_CONST:
      ret = star_tune(phys);
    break;
    /*case FORCE_BALANCE:
      ret = star_tunes(phys);
    break;
    case FIX_CENTER:
      ret = star_tunes(phys);
    break;
    case FIND_SURFACE:
      ret = star_tunes(phys);
    break;
    case EXTRAPOLATE:
      ret = star_tunes(phys);
    break;
    case AH_RADIUS:
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
  
  /* set to 0 to catch bug and miss you */
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
    //case STRESS_ENERGY:
      //ret = Tij_mount(phys);
    //break;
    /*case EULER_CONST:
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
  
  Error0("/* what should i do for BHNS region? */");
  
  assert(phys->region);
  
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
  else
    Error0(NO_OPTION);
  
  switch(type)
  {
    case NS:
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
  _free(phys);
}

/* set phys->region, if it's already set by the user skip  */
static void set_phys_default_region(Physics_T *const phys)
{
  if (phys->region)
    return;
  
  /* set natural region, later one can make it more sophisticated */
  phys->region = phys->stype;
    
}

/* ->: stype with correct indices.
// since we have different indices for objects like NS1, BH2 etc
// and some function are very similar but the stype indices are different
// this function gets a prototype and with respect to the given
// physics adjust the correct stype. eg:
// if phys->type = NS2 => 
// phy_return_correct_stype(phys,"NS_around") = "NS2_around". */
const char *phys_return_correct_stype(Physics_T *const phys,
                               const char *const stype)
{
  if (strcmp_i(stype,"NS"))
  {
    if (phys->type == NS)
      return "NS";
    else if (phys->type == NS1)
      return "NS1";
    else if (phys->type == NS2)
      return "NS2";
    else
      Error0(NO_OPTION);
    
  }
  else if (strcmp_i(stype,"NS_around"))
  {
    if (phys->type == NS)
      return "NS_around";
    else if (phys->type == NS1)
      return "NS1_around";
    else if (phys->type == NS2)
      return "NS2_around";
    else
      Error0(NO_OPTION);
    
  }
  else if (strcmp_i(stype,"NS_around_IB"))
  {
    if (phys->type == NS)
      return "NS_around_IB";
    else if (phys->type == NS1)
      return "NS1_around_IB";
    else if (phys->type == NS2)
      return "NS2_around_IB";
    else
      Error0(NO_OPTION);
  }
  else if (strcmp_i(stype,"NS_around_OB"))
  {
    if (phys->type == NS)
      return "NS_around_OB";
    else if (phys->type == NS1)
      return "NS1_around_OB";
    else if (phys->type == NS2)
      return "NS2_around_OB";
    else
      Error0(NO_OPTION);
  }
  else if (strcmp_i(stype,"BH"))
  {
    if (phys->type == BH)
      return "BH";
    else if (phys->type == BH1)
      return "BH1";
    else if (phys->type == BH2)
      return "BH2";
    else
      Error0(NO_OPTION);
    
  }
  else if (strcmp_i(stype,"BH_around"))
  {
    if (phys->type == BH)
      return "BH_around";
    else if (phys->type == BH1)
      return "BH1_around";
    else if (phys->type == BH2)
      return "BH2_around";
    else
      Error0(NO_OPTION);
    
  }
  else if (strcmp_i(stype,"BH_around_IB"))
  {
    if (phys->type == BH)
      return "BH_around_IB";
    else if (phys->type == BH1)
      return "BH1_around_IB";
    else if (phys->type == BH2)
      return "BH2_around_IB";
    else
      Error0(NO_OPTION);
  }
  else if (strcmp_i(stype,"BH_around_OB"))
  {
    if (phys->type == BH)
      return "BH_around_OB";
    else if (phys->type == BH1)
      return "BH1_around_OB";
    else if (phys->type == BH2)
      return "BH2_around_OB";
    else
      Error0(NO_OPTION);
  }
  else if (strcmp_i(stype,"NS1"))
  {
    if (phys->type == NS1)
      return "NS1";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"NS2"))
  {
    if (phys->type == NS2)
      return "NS2";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"NS1_around"))
  {
    if (phys->type == NS1)
      return "NS1_around";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"NS2_around"))
  {
    if (phys->type == NS2)
      return "NS2_around";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"NS1_around_IB"))
  {
    if (phys->type == NS1)
      return "NS1_around_IB";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"NS2_around_IB"))
  {
    if (phys->type == NS2)
      return "NS2_around_IB";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"NS1_around_OB"))
  {
    if (phys->type == NS1)
      return "NS1_around_OB";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"NS2_around_OB"))
  {
    if (phys->type == NS2)
      return "NS2_around_OB";
    else
      Error0("Wrong stype!");
  }
  
  else if (strcmp_i(stype,"BH1"))
  {
    if (phys->type == BH1)
      return "BH1";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"BH2"))
  {
    if (phys->type == BH2)
      return "BH2";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"BH1_around"))
  {
    if (phys->type == BH1)
      return "BH1_around";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"BH2_around"))
  {
    if (phys->type == BH2)
      return "BH2_around";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"BH1_around_IB"))
  {
    if (phys->type == BH1)
      return "BH1_around_IB";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"BH2_around_IB"))
  {
    if (phys->type == BH2)
      return "BH2_around_IB";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"BH1_around_OB"))
  {
    if (phys->type == BH1)
      return "BH1_around_OB";
    else
      Error0("Wrong stype!");
  }
  else if (strcmp_i(stype,"BH2_around_OB"))
  {
    if (phys->type == BH2)
      return "BH2_around_OB";
    else
      Error0("Wrong stype!");
  }
  else
    Error1("No replacement for %s\n",stype);
  
  return 0;
}

