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
  Uint i;
  Coord_T coord = patch->coordsys;
  Flag_T flg = NONE;
  
  if (!func) Error0("The func is null!\n");
  
  FOR_ALL(i,func)
    if (strcmp_i(func[i]->task,task) && func[i]->coord == coord)
    {
      func[i]->f(patch);
      flg = FOUND;
    }
  
  if (flg != FOUND)
    Errors("There is not %s task.\n",task);
}

/* initiatin a sFunc_Grid2Pdouble_T struct */
void init_func_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const func)
{
  *func = 0;
}

/* add a general grid to pointer to double function to func struct*/
void add_func_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const func,
                    double *(*f)(Patch_T *const patch),
                      const char *const name)
{
  sFunc_Patch2Pdouble_T *new_func;
  
  new_func = alloc_sFunc_Patch2Pdouble(func);
  new_func->func = f;
  new_func->name = dup_s(name);
  new_func->flg = 0;
}

/* given name and a data base of grid to pointer to double function 
// it finds the function with that name and return the a pointer to
// its structure.
// ->return value: found sFunc_Patch2Pdouble_T *, null otherwise.
*/
sFunc_Patch2Pdouble_T *get_func_Patch2Pdouble(const char *const name,
                                sFunc_Patch2Pdouble_T **const func)
{
  Uint i;
  
  if (!func) return 0;
  
  FOR_ALL(i,func)
    if (strcmp_i(func[i]->name,name))
    {
      return func[i];
    }
  
  return 0;
}

/* allocating 2 block of memory for sFunc_Patch2Pdouble_T 
// and putting the last block to NULL and returning
// the new available pointer.
// ->return value: a pointer to a ready sFunc_Patch2Pdouble
*/
void *alloc_sFunc_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const mem)
{
  Uint i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  IsNull((*mem));
  
  (*mem)[i] = calloc(1,sizeof(*(*mem)[i]));
  IsNull((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}

/* free function structure form patch to void */
void free_func_PtoV(sFunc_PtoV_T **func)
{
  Uint i;
  
  for (i = 0; func[i] != 0; ++i)
  {
    free(func[i]->task);
    free(func[i]);
  }
  
  free(func);
}

/* allocating 2 block of memory for sFunc_PtoV_T 
// and putting the last block to NULL and returning
// the new available pointer.
// ->return value: a pointer to a ready sFunc_PtoV
*/
void *alloc_sFunc_PtoV(sFunc_PtoV_T ***const mem)
{
  Uint i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  IsNull((*mem));
  
  (*mem)[i] = malloc(sizeof(*(*mem)[i]));
  IsNull((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}




