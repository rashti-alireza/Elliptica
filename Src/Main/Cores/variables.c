/*
// Alireza Rashti
// July 2018
*/

#include "variables.h"

/* add variable with the specified name and attribute to
// the pool in the given patch.
// if alloc_flg == YES, it also allocates memroy for v on the patch.
// note: if attribute is null, the variable won't have any attribute.
// ->return value: a pointer to the new made variable
*/
Variable_T *add_variable(const char *const name,const char *const attribute,Patch_T *const patch,const Flag_T alloc_flg)
{
  assert(name);
  assert(patch);
  
  Variable_T *var = 0;
  
  if (LookUpVar(name,patch) >= 0)
    abortEr_s("There is already a variable with the same name %s.\n",name);
  else
  {  
    
    var = calloc(1,sizeof(*var));
    pointerEr(var);
    var->patch = patch;
    add_attribute(var,attribute);
    
    if (alloc_flg == YES)
    {
      const unsigned nn = total_nodes_patch(patch);
      var->v = calloc(nn,sizeof(*var->v));
      pointerEr(var->v);
    }
    
    patch->pool = 
      realloc(patch->pool,(patch->nv+1)*sizeof(*patch->pool));
    pointerEr(patch->pool);
    patch->pool[patch->nv] = var;
    ++patch->nv;
  }
  
  return var;
}

/* remove the variable with the given name from the pool.
// and then shrink the pool.
// NOTE: THE INDEX OF VARIABLES ARE DYNAMIC,
// SO IT IS UNSAFE TO SAVE AN INDEX AND USED IT LATER. ONLY THE VALUES AND
// POINTER TO VALUES ARE STATIC AND SAFE TO SAVE.
*/
void remove_variable(const char *const name,Patch_T *const patch)
{
  /* if patch has no variable */
  if (!patch->nv) return;
  
  const int ind_remove = LookUpVar(name,patch);
  Variable_T *v_remove = patch->pool[ind_remove];
  Variable_T *v_last   = patch->pool[patch->nv-1];
  patch->pool[ind_remove] = v_last;
 
  free_variable(v_remove);
  patch->pool = realloc(patch->pool,(patch->nv)*sizeof(*patch->pool));
  pointerEr(patch->pool);
}

/* adding an attribute or info to var */
void add_attribute(Variable_T *const var,const char *const attribute)
{
  unsigned l1 = 0,l2 = 1;
  
  assert(var);
  /* if no attribute, return */
  if (!attribute)
    return;
  /* if the attribute already exists, return */
  if (strstr(var->info,attribute))
    return;
  
  if (var->info)
    l1 = (unsigned)strlen(var->info);
  l2 += (unsigned)strlen(attribute);
  
  var->info = realloc(var->info,l1+l2);
  pointerEr(var->info);
  var->info[l1] = '\0';
  strcat(var->info,attribute);
  
}

/* given name and patch find the index of a variable in the pool.
// ->return value: index of variable in the pool. INT_MIN if not found.
*/
int LookUpVar(const char *const name,Patch_T *const patch)
{
  int ind = INT_MIN;
  int i;
  
  for (i = 0; i < (int)patch->nv; ++i )
    if (!strcmp(patch->pool[i]->name,name))
      ind = i;
  
  return ind;
}
