/*
// Alireza Rashti
// May 2018
*/

#include "memory_free.h"

/* freeing 2 dimensions block of memory
// knowing the last column is pointing to null
*/
void free_2d(void *mem0)
{
  int i;
  void **mem = mem0;
    
  for (i = 0; mem[i] != 0; i++)
    free(mem[i]);
    
  free(mem);
    
}

/* freeing 2 dimensions block of memory
// knowing the number of columns is c
*/
void free_matrix(void *mem0, const unsigned long c)
{
  unsigned long i;
  void **mem = mem0;
  
  for (i = 0; i < c; i++)
    free(mem[i]);
    
  free(mem);
    
}

/* free needle */
void free_needle(Needle_T *needle)
{
  if (needle == 0) return;
  
  else
  {
    if (needle->Nin  != 0) free (needle->in);
    if (needle->Nex  != 0) free (needle->ex);
    if (needle->Ng   != 0) free (needle->guess);
    if (needle->Nans != 0) free (needle->ans);
    
  }
  
  free(needle);
}

/* feeing memory of Point_T inside grid */
void free_points(Grid_T *const grid)
{
  unsigned pa;
  
  FOR_ALL(pa,grid->patch)
  {
    Interface_T **face = grid->patch[pa]->interface;
    unsigned f;
    
    FOR_ALL(f,face)
    {
      free_matrix(face[f]->point,face[f]->np);
      face[f]->np = 0;
      face[f]->point = 0;
    }

  }
}

/* free function structure form patch to void */
void free_func_PtoV(sFunc_PtoV_T **func)
{
  unsigned i;
  
  for (i = 0; func[i] != 0; ++i)
  {
    free(func[i]->task);
    free(func[i]);
  }
  
  free(func);
}

/* free function structure form grid to pointer to double */
void free_func_Patch2Pdouble(sFunc_Patch2Pdouble_T **func)
{
  unsigned i;
  
  for (i = 0; func[i] != 0; ++i)
  {
    free(func[i]->name);
    free(func[i]);
  }
  
  free(func);
}

/* freeing v2 of field and put it to 0 */
void free_v2(Field_T *f)
{
  if (f->v2)
    free(f->v2);
  f->v2 = 0;
}

/* freeing a variable */
void free_field(Field_T *fld)
{
  if (!fld)
    return;
  
  if (fld->name)
    free(fld->name);
  if (fld->v)
    free(fld->v);
  if (fld->v2)
    free(fld->v2);
  if (fld->info)
    free(fld->info);
    
  free(fld);
}
