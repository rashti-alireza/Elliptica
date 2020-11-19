/*
// Alireza Rashti
// November 2020
*/

/* utilities for various purposes */

#include "man_utils.h"

Obj_Man_T *
init_obj_man
 (
 Grid_T *const grid/* computation grid */,
 const Com_Obj_T type/* object type NS,BH,etc */
 )
{
  Obj_Man_T *obj = calloc(1,sizeof(*obj)); IsNull(obj);
  const char *spos  = 0;
  
  /* what should i do for BHNS? */
  assert(obj->region);
  
  obj->grid = grid;
  obj->type = type;
  
  if (Pcmps("project","BH_NS_initial_data"))
  {
    obj->sys = BHNS;
  }
  else if (Pcmps("project","NS_NS_initial_data"))
  {
    obj->sys = NSNS;
  }
  else if (Pcmps("project","BH_BH_initial_data"))
  {
    obj->sys = BHBH;
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
      obj->pos = NONE;
      
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

