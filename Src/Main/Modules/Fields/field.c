/*
// Alireza Rashti
// July 2018
*/

#include "field.h"
#define MAX_STR 400
#define PR_FORMAT "(p%u,d%u,c%d,b%d)"

/* add field with the specified name and attribute to
// the pool in the given patch.
// if alloc_flg == YES, it also allocates memroy for v on the patch.
// note: if attribute is null, the field attribute is (3dim).
// note: one can simulate 2-d fields by replicating the value of
// the field in 2-d for each value of the coordiantes which this field doesn't depend on it.
// for example,
// if the field is 2-d like F = F(y,z) one can populate F = F(i,y,z) for all i's which
// is 3-d but the x direction is not effective.
// ->return value: a pointer to the new made field
*/
Field_T *add_field(const char *const name,const char *attribute,Patch_T *const patch,const Flag_T alloc_flg)
{
  assert(name);
  assert(patch);
  
  Field_T *fld = 0;
  
  if (alloc_flg != NO && alloc_flg != YES)
    abortEr("Wrong Flag was used. Flag for allocation is either YES or NO.\n");
  
  if (LookUpField(name,patch) >= 0)
    abortEr_s("There is already a field with the same name \"%s\".\n",name);
  else
  {  
    
    fld = calloc(1,sizeof(*fld));
    pointerEr(fld);
    fld->patch = patch;
    
    fld->name = dup_s(name);
    
    /* if attribute is null */
    if (!attribute)
      attribute = "(3dim)";
    add_attribute(fld,attribute);
    
    if (alloc_flg == YES)
    {
      const unsigned nn = total_nodes_patch(patch);
      fld->v = calloc(nn,sizeof(*fld->v));
      pointerEr(fld->v);
    }
    
    patch->pool = 
      realloc(patch->pool,(patch->nfld+1)*sizeof(*patch->pool));
    pointerEr(patch->pool);
    patch->pool[patch->nfld] = fld;
    ++patch->nfld;
  }
  
  return fld;
}

/* remove the given field from the pool and then shrink the pool.
// NOTE: THE INDEX OF VARIABLES ARE DYNAMIC,
// SO IT IS UNSAFE TO SAVE AN INDEX AND USED IT LATER. ONLY THE VALUES AND
// POINTER TO VALUES ARE STATIC AND SAFE TO SAVE.
*/
void remove_field(Field_T *f)
{
  if (!f) return;
  
  Patch_T *const patch = f->patch;
  char *const name = f->name;
  
  const int remove_ind = LookUpField(name,patch);
  /* if no such field exists */
  if (remove_ind < 0) 
    return;
  
  Field_T *remove_fld = patch->pool[remove_ind];
  Field_T *last_fld   = patch->pool[patch->nfld-1];
  patch->pool[remove_ind] = last_fld;
 
  free_field(remove_fld);
  patch->pool = realloc(patch->pool,(patch->nfld)*sizeof(*patch->pool));
  pointerEr(patch->pool);
  
  --patch->nfld;
}

/* adding an attribute or info to fld. each attribute must be writen in
// format = (attribute).
// attribute could be: (3dim), (2dim) etc.
*/
void add_attribute(Field_T *const fld,const char *const attribute)
{
  unsigned l1 = 0,l2 = 1;
  
  assert(fld);
  /* if no attribute, return */
  if (!attribute)
    return;
  if(!strchr(attribute,'(') || !strchr(attribute,')'))
    abortEr_s("Each attribute must be written in parentheses.\n"
    "This attribute %s doesn't have.\n",attribute);
    
  /* if the attribute already exists, return */
  if (fld->attr)
    if (strstr(fld->attr,attribute))
      return;
  
  if (fld->attr)
    l1 = (unsigned)strlen(fld->attr);
  l2 += (unsigned)strlen(attribute);
  
  fld->attr = realloc(fld->attr,l1+l2);
  pointerEr(fld->attr);
  fld->attr[l1] = '\0';
  strcat(fld->attr,attribute);
  
}

/* given name and patch find the index of a field in the pool.
// ->return value: index of field in the pool. INT_MIN if doesn't exist. */
int LookUpField(const char *const name,const Patch_T *const patch)
{
  int ind = INT_MIN;
  int i;
  
  if (!patch->pool)
    return ind;
    
  for (i = 0; i < (int)patch->nfld; ++i )
  {
    if (!patch->pool[i])
      continue;
    else if(!patch->pool[i]->name)
      continue;
    else if (!strcmp(patch->pool[i]->name,name))
    {
      ind = i;
      break;
    }
  }
  
  return ind;
}

/* given name and patch find the index of a field in the pool.
// ->return value: index of field in the pool, OR gives error if it could not find it. */
int LookUpField_E(const char *const name,const Patch_T *const patch)
{
  int ind = LookUpField(name,patch);
  
  if (ind == INT_MIN)
  {
    fprintf(stderr,"Field %s could not be found in patch %s.\n",name,patch->name);
    abortEr("Field could not be found.\n");
  }
  
  return ind;
}

/* values -> 1-d coeffs for the specified direction.
// making coeffs of a field based on basis used for expansion.
// arguments of function:
// ======================
//
// o. f is demanded field
// o. dir refers to directions which this coefficients may be found.
//	0 for a-direction, 1 for b-direction, 2 for c-direction.
//
// note1: each patch uses it own basis type and each part of coeffs
// corresponds to that basis type.
// note2: it DOES ALLOCATED MEMORY for coeffs on the "the patch" if not available.
// and puts it inside f->v2
// ->return value: coeffs
*/
double *make_coeffs_1d(Field_T *const f,const unsigned dir)
{
  assert(f);
  assert(dir < 3);/* note that we cannot have more that 3 direction, 
                  // since the grid is 3 dimension.
                  */
  assert(strstr(f->attr,"(3dim)"));
  
  const unsigned N = f->patch->nn;
  
  if (IsAvailable_1d(f,dir))
    return f->v2;
  else
  {
    if (!f->v2)
      f->v2 = alloc_double(N);
    
    find_1d_coeffs_in_patch(f,dir);
  }
  
  return f->v2;
}

/* values -> 2-d coeffs in specified direction and patch.
// for more info refere to make_coeffs_1d
// ->return value: coeffs.
*/
double *make_coeffs_2d(Field_T *const f,const unsigned dir1,const unsigned dir2)
{
  assert(dir1 < 3);
  assert(dir2 < 3);
  assert(strstr(f->attr,"(3dim)"));
  
  if (dir1 == dir2)
    abortEr_s("What does it mean to have twice transformation "
        "in the same direction for the field \"%s\"!\n",f->name);
    
  const unsigned N = f->patch->nn;
  int dir;
  Patch_T p_tmp1 = make_temp_patch(f->patch);
  Field_T *f_tmp = add_field("f_tmp","(3dim)",&p_tmp1,NO);
  
  if (IsAvailable_2d(f,dir1,dir2,&dir))
  {
    /* if field is ready ready */
    if (dir == -1)
    {
      remove_field(f_tmp);
      free_temp_patch(&p_tmp1);
      
      return f->v2;
    }
    else/* it only make coeffs in dir since the other one in ready */
    {  
      f_tmp->v = f->v2;
      free_attr(f_tmp);
      f_tmp->attr = dup_s(f->attr);
      f_tmp->info = dup_s(f->info);
      
      /* since f_tmp->v2 is empty we have to allocate memory for that */
      f_tmp->v2 = alloc_double(N);
      f->v2 = find_1d_coeffs_in_patch(f_tmp,(unsigned)dir);
      
    }
  }
  else
  {
    if (!f->v2)
      f->v2 = alloc_double(N);
    
    find_1d_coeffs_in_patch(f,dir1);
    f_tmp->v = f->v2;
    free(f_tmp->attr);
    f_tmp->attr = dup_s(f->attr);
    f_tmp->info = dup_s(f->info);
      
    /* since f_tmp->v2 is empty it will be made by the functions below */
    f_tmp->v2 = alloc_double(N);
    f->v2 = find_1d_coeffs_in_patch(f_tmp,dir2);
  }
  
  free_attr(f);
  free_info(f);
  f->attr = dup_s(f_tmp->attr);
  f->info = dup_s(f_tmp->info);
  free_v(f_tmp);
  f_tmp->v = 0;
  f_tmp->v2 = 0;
  remove_field(f_tmp);
  free_temp_patch(&p_tmp1);
  
  return f->v2;
}

/* values -> 3-d coeffs in specified patch.
// for more info refere to make_coeffs_1d
// ->return value: coeffs.
*/
double *make_coeffs_3d(Field_T *const f)
{
  const unsigned N = f->patch->nn;
  
  assert(strstr(f->attr,"(3dim)"));
  
  if (IsAvailable_3d(f))
  {
    return f->v2;
  }
  else
  {
    if (!f->v2)
      f->v2 = alloc_double(N);
    
    Patch_T p_tmp1 = make_temp_patch(f->patch);
    Patch_T p_tmp2 = make_temp_patch(f->patch);
    Field_T *f_tmp1 = add_field("f_tmp1","(3dim)",&p_tmp1,NO);
    Field_T *f_tmp2 = add_field("f_tmp2","(3dim)",&p_tmp2,NO);
 
    f_tmp1->v = find_1d_coeffs_in_patch(f,0);
    free_attr(f_tmp1);
    f_tmp1->attr = dup_s(f->attr);
    free_info(f_tmp1);
    f_tmp1->info = dup_s(f->info);
    /* f_tmp1->v == f->v2 
    // f_tmp1->v2 == 0
    */
    f_tmp1->v2 = alloc_double(N);
    f_tmp2->v = find_1d_coeffs_in_patch(f_tmp1,1);
    free_attr(f_tmp2);
    f_tmp2->attr = dup_s(f_tmp1->attr);
    free_info(f_tmp2);
    f_tmp2->info = dup_s(f_tmp1->info);
    /* f_tmp2->v == f_tmp1->v2
    // f_tmp2->v2 == 0
    */
    free_v2(f);/* => free(f_tmp1->v)*/
    f_tmp2->v2 = alloc_double(N);
    f->v2      = find_1d_coeffs_in_patch(f_tmp2,2);
    free_attr(f);
    f->attr = dup_s(f_tmp2->attr);
    free_info(f);
    f->info = dup_s(f_tmp2->info);
    /* f_tmp2->v == f_tmp1->v2
    // f_tmp2->v2 == f->v2
    */
    free_v2(f_tmp1);
    f_tmp1->v = 0;
    f_tmp2->v = 0;
    f_tmp2->v2 = 0;
    remove_field(f_tmp1);
    remove_field(f_tmp2);
    free_temp_patch(&p_tmp1);
    free_temp_patch(&p_tmp2);
  }
  
  return f->v2;
}

/* check if this coeffs has been already made.
// furthermore if the coeffs are not appropriate, clean them.
// ->return value: 1 if already made, 0 otherwise.
*/
static unsigned IsAvailable_1d(Field_T *const f,const unsigned dir)
{
  Patch_T *const patch = f->patch;
  const Collocation_T collocation = patch->collocation[dir];
  const Basis_T basis = patch->basis[dir];
  unsigned r = 0;
  const unsigned pn = patch->pn;
  char needle[3][MAX_STR] = {{'\0'}};
  Flag_T flg = CLEAN;
  
  sprintf(needle[0],PR_FORMAT,pn,dir,collocation,basis);
  sprintf(needle[1],"p%u,d%u",pn,(dir+1)%3);
  sprintf(needle[2],"p%u,d%u",pn,3-dir-(dir+1)%3);
  
  /* if there is no coeffs */
  if (!f->v2 || !f->info) 
  {
    flg = CLEAN;
    r = 0;
  }
  /* if other direction exist one needs to clean coeffs */
  else if (strstr(f->info,needle[1]) || strstr(f->info,needle[2]))
  {
    flg = CLEAN;
    r = 0;
  }
  else if (strstr(f->info,needle[0]))
  {
    flg = NO;
    r = 1;
  }
  else/* if some mismatches exist */
  {
    flg = CLEAN;
    r = 0;
  }
  
  if (flg == CLEAN)
  {
    free_v2(f);
    free_info(f);
  }
  
  return r;
}

/* check if this coeffs has been already made.
// if only one of them is available put the other one in dir.
// ->return value: 1 if already made, 0 otherwise.
// ->return value: if 1 and dir == -1, it means the coeffs is ready ready, 
// otherwise fill dir in required direction.
*/
static unsigned IsAvailable_2d(Field_T *const f,const unsigned dir1,const unsigned dir2,int *dir)
{
  Patch_T *const patch = f->patch;
  Collocation_T collocation[3];
  Basis_T basis[3];
  unsigned r = 0;
  unsigned pn = f->patch->pn;
  char needle[3][MAX_STR] = {{'\0'}};
  unsigned DIR[3];
  Dd_T e;
  Flag_T flg = CLEAN, pristine = YES;
  
  *dir = INT_MAX;
  collocation[dir1] = patch->collocation[dir1];
  collocation[dir2] = patch->collocation[dir2];
  basis[dir1]       = patch->basis[dir1];
  basis[dir2]       = patch->basis[dir2];
  
  /* finding the unrelated directions */
  for (e = _N0_; e <= _N2_; ++e)
  {
    if (e == dir1 || e == dir2)
    {
      DIR[e] = 0;
      sprintf(needle[e],PR_FORMAT,pn,e,collocation[e],basis[e]);
    }
    else
    {
      DIR[e] = 1;/* means unrelated direction */
      sprintf(needle[e],"p%u,d%u",pn,e);
    }
  }
  
  /* if there is no coeffs */
  if (!f->v2 || !f->info)
  {
    flg = CLEAN;
    r = 0;
  }
  else
  {
    /* check the coeffs be pristine */
    for (e = _N0_; e <= _N2_; ++e)
    {
      /* for unrelated */
      if (DIR[e])
      {
        /* if has unrelated */
        if (strstr(f->info,needle[e]))
        {
          flg = CLEAN;
          pristine = NO;
          r = 0;
          break;
        }
      }
    }/* end of for (e = 0 ; e < 3; ++e) */
    
    if (pristine == YES)
    {
      if (strstr(f->info,needle[dir1]) && 
          strstr(f->info,needle[dir2])   )
      {
        r = 1;
        flg = NO;
        *dir = -1;
      }
      else if (strstr(f->info,needle[dir1]))
      {
        r = 1;
        flg = NO;
        *dir = (int)dir2;
      }
      else if (strstr(f->info,needle[dir2]))
      {
        r = 1;
        flg = NO;
        *dir = (int)dir1;
      }
      else
      {
        flg = CLEAN;
        r = 0;
      }
      
    }/* end of if (pristine == YES) */
  }/* end of else */
  
  if (flg == CLEAN)
  {
    free_v2(f);
    free_info(f);
  }
  
  return r;
}

/* check if this coeffs has been already made.
// ->return value: 1 if already made, 0 otherwise.
*/
static unsigned IsAvailable_3d(Field_T *const f)
{
  Patch_T *const patch = f->patch;
  Collocation_T collocation[3];
  Basis_T basis[3];
  unsigned r = 0;
  unsigned pn = patch->pn;
  char needle[3][MAX_STR] = {{'\0'}};
  Flag_T flg = CLEAN;
  
  collocation[0] = patch->collocation[0];
  collocation[1] = patch->collocation[1];
  collocation[2] = patch->collocation[2];
  basis[0]       = patch->basis[0];
  basis[1]       = patch->basis[1];
  basis[2]       = patch->basis[2];
  
  sprintf(needle[0],PR_FORMAT,pn,0,collocation[0],basis[0]);
  sprintf(needle[1],PR_FORMAT,pn,1,collocation[1],basis[1]);
  sprintf(needle[2],PR_FORMAT,pn,2,collocation[2],basis[2]);
  
  /* if there is no coeffs */
  if (!f->v2 || !f->info)
  {
    flg = CLEAN;
    r = 0;
  }
  
  /* if coeffs are ready */
  else if (strstr(f->info,needle[0]) && strstr(f->info,needle[1]) 
      && strstr(f->info,needle[2]))
  {
    flg = NO;
    r = 1;
  }
  
  if (flg == CLEAN)
  {
    free_v2(f);
    free_info(f);
    r = 0;
  }
  
  return r;
}

/* finding coeffs of expansion in direction dir 
// for a given field, patch and direction.
// note: it DOESN'T ALLOCATE MEMORY.
// ->return value: coeffs
*/
static double *find_1d_coeffs_in_patch(Field_T *const f,const unsigned dir)
{
  const Collocation_T collocation = f->patch->collocation[dir];
  const Basis_T basis = f->patch->basis[dir];
  
  if(basis == Chebyshev_Tn_BASIS)
  {
    if (collocation == Chebyshev_Extrema)
    {
      coeffs_patch_Tn_Extrema_1d(f,dir);
      add_Tinfo(f,dir,Chebyshev_Extrema,Chebyshev_Tn_BASIS);
    }
    else if (collocation == Chebyshev_Nodes)
    {
      coeffs_patch_Tn_Nodes_1d(f,dir);
      add_Tinfo(f,dir,Chebyshev_Nodes,Chebyshev_Tn_BASIS);
    }
    else
      abortEr("There is no such COLLOCATION defined for this function.\n");
  }/* end of if(basis == Chebyshev_Tn_BASIS) */
  
  else if (basis == No_BASIS)
  {
    abortEr("This function is not have this part!\n");
  }
  
  else
    abortEr("There is no such BASIS defined for this function.\n");
    
  return f->v2;
}

/* finding coeffs in in patch For Tn basis 
// with Chebyshev extrema collocation.
*/
static void coeffs_patch_Tn_Extrema_1d(Field_T *const f,const unsigned dir)
{
  assert(f->v2);
  assert(f->v);
  Fourier_Transformation_1d_F *FourierTrans;
  double *const coeffs = f->v2;
  const double *const values = f->v;
  double *out, *in;
  const unsigned *n = f->patch->n;
  const unsigned B = n[(dir+1)%3]*n[3-dir-(dir+1)%3];
  unsigned l,i,j,k,s;
  out = alloc_double(n[dir]);
  in = alloc_double(n[dir]);
  
  if (strstr_i(PgetsEZ("Fourier_Transformation_Method"),"RFT"))
    FourierTrans = rft_1d_ChebyshevExtrema_coeffs;
  else
    abortEr("No such Fourier_Transformation_Method defined for this function.\n");
  
  if (dir == 0)
  {
    for (s = 0; s < B; ++s)
    {
      /* s = k+n2*j */
      j = s/n[2];
      k = s%n[2];
      
      for (i = 0; i < n[dir]; ++i)
      {
        l = L(n,i,j,k);
        in[i] = values[l];
      }
      FourierTrans(in,out,n[dir]);
      
      for (i = 0; i < n[dir]; ++i)
      {
        l = L(n,i,j,k);
        coeffs[l] = out[i];
      }
    }
  }/* end of if (dir == 0) */
  
  else if (dir == 1)
  {
    for (s = 0; s < B; ++s)
    {
      /* s = k+n2*i */
      i = s/n[2];
      k = s%n[2];
      
      for (j = 0; j < n[dir]; ++j)
      {
        l = L(n,i,j,k);
        in[j] = values[l];
      }
      FourierTrans(in,out,n[dir]);
      
      for (j = 0; j < n[dir]; ++j)
      {
        l = L(n,i,j,k);
        coeffs[l] = out[j];
      }
    }
  }/* end of if (dir == 1) */
  
  else if (dir == 2)
  {
    for (s = 0; s < B; ++s)
    {
      /* s = j+n1*i */
      i = s/n[1];
      j = s%n[1];
      
      for (k = 0; k < n[dir]; ++k)
      {
        l = L(n,i,j,k);
        in[k] = values[l];
      }
      FourierTrans(in,out,n[dir]);
      
      for (k = 0; k < n[dir]; ++k)
      {
        l = L(n,i,j,k);
        coeffs[l] = out[k];
      }
    }
  }/* end of if (dir == 2) */
  else
    abortEr("direction must be 0,1 or 2.\n");
    
  free(out);
  free(in);
}

/* finding coeffs in in patch For Tn basis 
// with Chebyshev nodes collocation.
*/
static void coeffs_patch_Tn_Nodes_1d(Field_T *const f,const unsigned dir)
{
  assert(f->v2);
  assert(f->v);
  Fourier_Transformation_1d_F *FourierTrans;
  double *const coeffs = f->v2;
  const double *const values = f->v;
  double *out, *in;
  const unsigned *n = f->patch->n;
  const unsigned B = n[(dir+1)%3]*n[3-dir-(dir+1)%3];
  unsigned l,i,j,k,s;
  out = alloc_double(n[dir]);
  in = alloc_double(n[dir]);
  
  if (strstr_i(PgetsEZ("Fourier_Transformation_Method"),"RFT"))
    FourierTrans = rft_1d_ChebyshevNodes_coeffs;
  else
    abortEr("No such Fourier_Transformation_Method defined for this function.\n");
  
  if (dir == 0)
  {
    for (s = 0; s < B; ++s)
    {
      /* s = k+n2*j */
      j = s/n[2];
      k = s%n[2];
      
      for (i = 0; i < n[dir]; ++i)
      {
        l = L(n,i,j,k);
        in[i] = values[l];
      }
      FourierTrans(in,out,n[dir]);
      
      for (i = 0; i < n[dir]; ++i)
      {
        l = L(n,i,j,k);
        coeffs[l] = out[i];
      }
    }
  }/* end of if (dir == 0) */
  
  else if (dir == 1)
  {
    for (s = 0; s < B; ++s)
    {
      /* s = k+n2*i */
      i = s/n[2];
      k = s%n[2];
      
      for (j = 0; j < n[dir]; ++j)
      {
        l = L(n,i,j,k);
        in[j] = values[l];
      }
      FourierTrans(in,out,n[dir]);
      
      for (j = 0; j < n[dir]; ++j)
      {
        l = L(n,i,j,k);
        coeffs[l] = out[j];
      }
    }
  }/* end of if (dir == 1) */
  
  else if (dir == 2)
  {
    for (s = 0; s < B; ++s)
    {
      /* s = j+n1*i */
      i = s/n[1];
      j = s%n[1];
      
      for (k = 0; k < n[dir]; ++k)
      {
        l = L(n,i,j,k);
        in[k] = values[l];
      }
      FourierTrans(in,out,n[dir]);
      
      for (k = 0; k < n[dir]; ++k)
      {
        l = L(n,i,j,k);
        coeffs[l] = out[k];
      }
    }
  }/* end of if (dir == 2) */
  else
    abortEr("direction must be 0,1 or 2.\n");
    
  free(out);
  free(in);
}


/* adding Transformation info in f->info.
// for each transformation the print format PR_FORMAT is
// (patch,direction,collocation,basis) and concatenated in row.
*/
static void add_Tinfo(Field_T *const f,const unsigned dir,const Collocation_T collocation,const Basis_T basis)
{
  char inf[MAX_STR] = {'\0'};
  unsigned l1 = 0,l2 = 0;
  const unsigned pn = f->patch->pn;
  
  if (f->info)
    l1 = (unsigned)strlen(f->info);
  /* (patch,direction,collocation,basis) */
  sprintf(inf,PR_FORMAT,pn,dir,collocation,basis);
  l2 = l1+(unsigned)strlen(inf)+1;
  f->info = realloc(f->info,l2);
  pointerEr(f->info);
  f->info[l1] = '\0';
  strcat(f->info,inf);
  
}

/* allocating 3-D fields on the whole grid to be solved 
// according to solution_man in each patch.
*/
void enable_fields(Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    char **field_name = patch->solving_man->field_name;
    unsigned nf = patch->solving_man->nf;
    unsigned i;
    
    for (i = 0; i < nf; ++i)
      add_field(field_name[i],"(3dim)",patch,YES);
  }
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

