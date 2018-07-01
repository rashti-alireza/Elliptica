/*
// Alireza Rashti
// June 2018
*/

#include "function_managing.h"

/* initializing func */
void init_func_PtoV(sFunc_PtoV_T ***const func)
{
  (*func) = 0;
}

/* add a general patch to void function to func struct*/
void add_func_PtoV(sFunc_PtoV_T ***const func,
                    void (*f)(Patch_T *const patch),
                      const char *const task,
                        Coord_T coord)
{
  sFunc_PtoV_T *new_func;
  
  new_func = alloc_sFunc_PtoV(func);
  new_func->f = f;
  new_func->task = dup_s(task);
  new_func->coord = coord;
}

/* execute func based on coord and task */
void run_func_PtoV(sFunc_PtoV_T **const func,const char *const task,Patch_T *const patch)
{
  int i;
  Coord_T coord = find_coord(patch);
  Flag_T flg = NONE;
  
  FOR_ALL(i,func)
    if (strcmp_i(func[i]->task,task) && func[i]->coord == coord)
    {
      func[i]->f(patch);
      flg = FOUND;
    }
  
  if (flg != FOUND)
    abortEr_s("There is not %s task.\n",task);
}

/* find coord enum based in patch */
Coord_T find_coord(const Patch_T *const patch)
{
  Coord_T coord;
  
  if(strcmp_i(patch->coordsys,"Cartesian"))
    coord = Cartesian;
  else
    abortEr_s("There is no such %s coordinates.\n",patch->coordsys);  
    
  return coord;
}
