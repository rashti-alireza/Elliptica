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
// knowing the number of rows is c
*/
void free_2d_mem(void *mem0, const unsigned long c)
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
      free_2d_mem(face[f]->point,face[f]->np);
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

/* freeing v of field and put it to 0 */
void free_v(Field_T *f)
{
  if (f->v)
    free(f->v);
  f->v = 0;
}


/* freeing info of field and put it to 0 */
void free_info(Field_T *f)
{
  if (f->info)
    free(f->info);
  f->info = 0;
}

/* freeing a the content of a given field */
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
  if (fld->attr)
    free(fld->attr);
  if (fld->info)
    free(fld->info);
    
  free(fld);
}

/* freeing v2 and info of a field */
void free_coeffs(Field_T *fld)
{
  free_info(fld);
  free_v2(fld);
}

/* freeing Equation_T data base */
void free_db_eqs(sEquation_T **db)
{
  unsigned ndb;
  
  if (!db) return;
  
  for (ndb = 0; db[ndb] != 0; ++ndb)
  {
    free(db[ndb]);
  }
  
  free(db);
}

/* freeing interpolation structure.
// NOTE: it frees 
// o. interp_s->point struct 
// o. interp_s
*/
void free_interpolation(Interpolation_T *interp_s)
{
  if (!interp_s)
    return;
    
  if (interp_s->point)
    free(interp_s->point);
    
  free(interp_s);
}

/* freeing matrix memroy */
void free_matrix(Matrix_T *m)
{
  if (!m)
    return;
    
  if (m->reg_f)
  {
    if (m->reg->A)
      free_2d_mem(m->reg->A,(long unsigned)m->row);
  }
  else if (m->tri_f)
  {
    if (m->tri->row)
      free(m->tri->row);
    if (m->tri->col)
      free(m->tri->col);
    if (m->tri->a)
      free(m->tri->a); 
  }
  else if (m->ccs_f)
  {
    if (m->ccs->Ap)
      free(m->ccs->Ap);
    if (m->ccs->Ai)
      free(m->ccs->Ai);
    if (m->ccs->Ax)
      free(m->ccs->Ax);
  }
  else if (m->crs_f)
  {
    if (m->crs->Ap)
      free(m->crs->Ap);
    if (m->crs->Aj)
      free(m->crs->Aj);
    if (m->crs->Ax)
      free(m->crs->Ax);
  }
  else if (m->tri_l_f)
  {
    if (m->tri_long->row)
      free(m->tri_long->row);
    if (m->tri_long->col)
      free(m->tri_long->col);
    if (m->tri_long->a)
      free(m->tri_long->a); 
  }
  else if (m->ccs_l_f)
  {
    if (m->ccs_long->Ap)
      free(m->ccs_long->Ap);
    if (m->ccs_long->Ai)
      free(m->ccs_long->Ai);
    if (m->ccs_long->Ax)
      free(m->ccs_long->Ax);
  }
  else if (m->crs_l_f)
  {
    if (m->crs_long->Ap)
      free(m->crs_long->Ap);
    if (m->crs_long->Aj)
      free(m->crs_long->Aj);
    if (m->crs_long->Ax)
      free(m->crs_long->Ax);
  }
  else
    abortEr("No matrix format is defined for this given matrix.\n");
    
  free(m);
}
