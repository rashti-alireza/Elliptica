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
  if (!mem0)
    return;
    
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
  if (!mem0)
    return;
    
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

/* freeing attr of field and put it to 0 */
void free_attr(Field_T *f)
{
  if (f->attr)
    free(f->attr);
  f->attr = 0;
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

/* freeing field */
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

/* free v, v2 and info to update values of a field */
void empty_field(Field_T *fld)
{
  if (!fld)
    return;
  
  if (fld->v)
    free(fld->v);
  if (fld->v2)
    free(fld->v2);
  if (fld->info)
    free(fld->info);
  
  fld->v    = 0;
  fld->v2   = 0;
  fld->info = 0; 
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

/* freeing interpolation structure. */
void free_interpolation(Interpolation_T *interp_s)
{
  if (!interp_s)
    return;
  if (strstr_i(interp_s->method,"Natural_Cubic_Spline_1D"))
  {
    if (interp_s->N_cubic_spline_1d->b)
      free(interp_s->N_cubic_spline_1d->b);
    if (interp_s->N_cubic_spline_1d->c)
      free(interp_s->N_cubic_spline_1d->c);
    if (interp_s->N_cubic_spline_1d->d)
      free(interp_s->N_cubic_spline_1d->d);
    if (interp_s->N_cubic_spline_1d->Alloc_Mem)
    {
      free(interp_s->N_cubic_spline_1d->x);
      free(interp_s->N_cubic_spline_1d->f);
    }
  }
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

/* free thoroughly patch->interface */
void free_patch_interface(Patch_T *const patch)
{
  Interface_T **face = patch->interface;
  SubFace_T **subface = 0;
  unsigned f,i;
  
  if (face)
  FOR_ALL(f,face)
  {
    /* free point */
    free_2d_mem(face[f]->point,face[f]->np);
    face[f]->np = 0;
    
    /* free subface */
    subface = face[f]->subface;
    for (i = 0; i < face[f]->ns; ++i)
    {
      _free(subface[i]->flags_str);
      _free(subface[i]->id);
      _free(subface[i]->adjid);
      _free(subface[i]);
    }
    _free(subface);
    
    /* free face */
    _free(face[f]);
  }
  _free(face);
  
  patch->interface = 0;
}

/* free thoroughly patch->interface */
void free_patch_SolMan_jacobian(Patch_T *const patch)
{
  Solving_Man_T *const SolMan = patch->solving_man;
  unsigned i;
  
  if (!SolMan)
    return;
  
  for (i = 0; i < SolMan->nj; ++i)
  {
    free_matrix(SolMan->jacobian[i]->J);
    _free(SolMan->jacobian[i]);
  }
  _free(SolMan->jacobian);
  
  SolMan->jacobian = 0;
  
}

/* free thoroughly patch->interface */
void free_patch_SolMan_method_Schur(Patch_T *const patch)
{
  DDM_Schur_Complement_T *s = patch->solving_man->method->SchurC;
  if (!s)
    return;
    
  Sewing_T **se             = s->sewing;
  Pair_T **p                = 0;
  unsigned i,j;
  
  /* note, those matrices and double populated during solver
  // won't needed to be freed */
  _free(s->map);
  _free(s->inv);
  _free(s->Imap);
  _free(s->Iinv);
  _free(s->NS_p);
  _free(s->NI_p);
  
  for (i = 0; i < s->nsewing; ++i)
  {
    p = 0;
    if (se[i])
    {
      if (se[i]->patchN != patch->pn)
      {
        _free(se[i]->map);
        _free(se[i]->inv);
        _free(se[i]->Imap);
        _free(se[i]->Iinv);
      }
      p = se[i]->pair;
      for (j = 0; j < se[i]->npair; ++j)
      {
        _free(p[j]->ip);
        _free(p[j]->nv);
        
        /* because only for the following patches we allocate memory */
        if (se[i]->patchN != patch->pn)
        {
          _free(p[j]->subface->flags_str);
          _free(p[j]->subface->id);
          _free(p[j]->subface->adjid);
          _free(p[j]->subface);  
          p[j]->subface = 0;
        }
        free(p[j]);
      }
      _free(p);
    }
    _free(se[i]);
  }
  _free(se);
  _free(s);
  
  patch->solving_man->method->SchurC = 0;
}

/* free only if p != NULL */
void _free(void *p)
{
  if (p)
    free(p);
}

/* free the given grid completely */
void free_grid(Grid_T *grid)
{
  unsigned p,ijk,nn,f;
  
  if (!grid)
    return;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    nn             = patch->nn;
    
    _free(patch->name);
    
    if (patch->coordsys != Cartesian)
      for (ijk = 0; ijk < nn; ++ijk)
        _free(patch->node[ijk]->X);
    
    if (patch->node)    
      free_2d_mem(patch->node,nn);
    
    for (f = 0; f < patch->nfld; ++f)
    {
      Field_T *field = patch->pool[f];
      free_field(field);
    }
    _free(patch->pool);
    _free(patch->JacobianT);
    free_patch_interface(patch);
    if (patch->solving_man)
    {
      free_patch_SolMan_jacobian(patch);
      free_patch_SolMan_method_Schur(patch);
      _free(patch->solving_man->field_eq);
      _free(patch->solving_man->bc_eq);
      _free(patch->solving_man->jacobian_field_eq);
      _free(patch->solving_man->jacobian_bc_eq);
      free_2d_mem(patch->solving_man->field_name,patch->solving_man->nf);
      free(patch->solving_man);
    }
  }
  free_2d_mem(grid->patch,grid->np);
  _free(grid->kind);
  free(grid);
}

/* free the patch completely */
void free_patch(Patch_T *patch)
{
  unsigned ijk,nn,f;
  
  if (!patch)
    return;
  
  nn = patch->nn;
  
  _free(patch->name);
  
  if (patch->coordsys != Cartesian)
    for (ijk = 0; ijk < nn; ++ijk)
      _free(patch->node[ijk]->X);
  
  if (patch->node)    
    free_2d_mem(patch->node,nn);
  
  for (f = 0; f < patch->nfld; ++f)
  {
    Field_T *field = patch->pool[f];
    free_field(field);
  }
  _free(patch->pool);
  _free(patch->JacobianT);
  free_patch_interface(patch);
  if (patch->solving_man)
  {
    free_patch_SolMan_jacobian(patch);
    free_patch_SolMan_method_Schur(patch);
    _free(patch->solving_man->field_eq);
    _free(patch->solving_man->bc_eq);
    _free(patch->solving_man->jacobian_field_eq);
    _free(patch->solving_man->jacobian_bc_eq);
    free_2d_mem(patch->solving_man->field_name,patch->solving_man->nf);
    free(patch->solving_man);
  }
  
  free(patch);
}

/* given the parameter name, free the parameter data base from it 
// and shrink the data base and put the last parameter in place of
// the deleted parameter. */
void free_parameter(const char *const par_name)
{
  Parameter_T *last_par = 0;
  unsigned np,i;
  
  /* count total number of parameters */
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
  
  if (np == 0)
    return;
    
  for (i = 0; i < np; ++i)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
    {
      last_par = parameters_global[np-1];
      
      _free(parameters_global[i]->lv);
      _free(parameters_global[i]->rv);
      _free(parameters_global[i]->rv_ip);
      _free(parameters_global[i]->rv_array);
      free(parameters_global[i]);
      
      parameters_global[i] = last_par;
      
      parameters_global = 
        realloc(parameters_global,np*sizeof(*parameters_global));
      pointerEr(parameters_global);
      
      parameters_global[np-1] = 0;
      break;
    }
  }
  
}

/* free the whole parameter data base */
void free_parameter_db(void)
{
  unsigned np;
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
  {
  
    _free(parameters_global[np]->lv);
    _free(parameters_global[np]->rv);
    _free(parameters_global[np]->rv_ip);
    _free(parameters_global[np]->rv_array);
    free(parameters_global[np]);
    np++;
  }
  
  _free(parameters_global);
}

