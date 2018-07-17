/*
// Alireza Rashti
// June 2018
*/

#include "function_managing.h"

/* initializing func patch to void */
void init_func_PtoV(sFunc_PtoV_T ***const func)
{
  (*func) = 0;
}

/* add a general patch to void function to func struct*/
void add_func_PtoV(sFunc_PtoV_T ***const func,
                    void (*f)(Patch_T *const patch),
                      const char *const task,
                        const Coord_T coord)
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
  unsigned i;
  Coord_T coord = patch->coordsys;
  Flag_T flg = NONE;
  
  if (!func) abortEr("The func is null!\n");
  
  FOR_ALL(i,func)
    if (strcmp_i(func[i]->task,task) && func[i]->coord == coord)
    {
      func[i]->f(patch);
      flg = FOUND;
    }
  
  if (flg != FOUND)
    abortEr_s("There is not %s task.\n",task);
}

/* initiatin a sFunc_Grid2Pdouble_T struct */
void init_func_Grid2Pdouble(sFunc_Grid2Pdouble_T ***const func)
{
  *func = 0;
}

/* add a general grid to pointer to double function to func struct*/
void add_func_Grid2Pdouble(sFunc_Grid2Pdouble_T ***const func,
                    double *(*f)(Grid_T *const grid),
                      const char *const name)
{
  sFunc_Grid2Pdouble_T *new_func;
  
  new_func = alloc_sFunc_Grid2Pdouble(func);
  new_func->func = f;
  new_func->name = dup_s(name);
  new_func->flg = 0;
}

/* given name and a data base of grid to pointer to double function 
// it finds the function with that name and return the a pointer to
// its structure.
// ->return value: found sFunc_Grid2Pdouble_T *, null otherwise.
*/
sFunc_Grid2Pdouble_T *get_func_Grid2Pdouble(const char *const name,
                                sFunc_Grid2Pdouble_T **const func)
{
  unsigned i;
  
  if (!func) return 0;
  
  FOR_ALL(i,func)
    if (strcmp_i(func[i]->name,name))
    {
      return func[i];
    }
  
  return 0;
}
