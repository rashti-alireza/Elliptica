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
// ->return value: index of field in the pool. INT_MIN if doesn't exist.
*/
int LookUpField(const char *const name,Patch_T *const patch)
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
      ind = i;
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
  assert(dir <= 2);
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
  assert(dir1 <= 2);
  assert(dir2 <= 2);
  assert(strstr(f->attr,"(3dim)"));
    
  const unsigned N = f->patch->nn;
  int dir;
  Field_T *f_tmp = add_field("f_tmp","(3dim)",f->patch,NO);
  
  if (IsAvailable_2d(f,dir1,dir2,&dir))
  {
    /* if field is ready ready */
    if (dir == -1)
      return f->v2;
    else/* it only make coeffs in dir since the other one in ready */
    {  
      f_tmp->v = f->v2;
      /* since f_tmp->v2 is empty it will be made by the function below */
      f->v2 = find_1d_coeffs_in_patch(f_tmp,(unsigned)dir);
      free_v(f_tmp);
    }
  }
  else
  {
    if (!f->v2)
      f->v2 = alloc_double(N);
    
    find_1d_coeffs_in_patch(f,dir1);
    f_tmp->v = f->v2;
    /* since f_tmp->v2 is empty it will be made by the function below */
    f->v2 = find_1d_coeffs_in_patch(f_tmp,dir2);
    free_v(f_tmp);
  }
  
  f_tmp->v = 0;
  f_tmp->v2 = 0;
  remove_field(f_tmp);
    
  return f->v2;
}

/* values -> 3-d coeffs in specified patch.
// for more info refere to make_coeffs_1d
// ->return value: coeffs.
*/
double *make_coeffs_3d(Field_T *const f)
{
  const unsigned N = f->patch->nn;
  Collocation_T collocation[3];
  Basis_T basis[3];
  
  collocation[0] = f->patch->collocation[0];
  collocation[1] = f->patch->collocation[1];
  collocation[2] = f->patch->collocation[2];
  basis[0]       = f->patch->basis[0];
  basis[1]       = f->patch->basis[1];
  basis[2]       = f->patch->basis[2];
  
  assert(strstr(f->attr,"(3dim)"));
  
  if (IsAvailable_3d(f))
  {
    return f->v2;
  }
  else
  {
    if (!f->v2)
      f->v2 = alloc_double(N);
    
    /* when all of directions and basis are the same.
    // it gets used for high performance.
    */
    if (Is3d_fft(collocation,basis))
    {
      double *coeffs = f->v2;
      double *values = f->v;
      const unsigned *n = f->patch->n;
    
      if (f->patch->basis[0] == Chebyshev_Tn_BASIS)
        fftw_3d_ChebyshevExtrema_coeffs(values,coeffs,n);
      else
        abortEr("No such basis is defined for this function.\n");
    }
    else
    {
      Field_T *f_tmp1 = add_field("f_tmp1","(3dim)",f->patch,NO);
      Field_T *f_tmp2 = add_field("f_tmp2","(3dim)",f->patch,NO);
 
      f_tmp1->v = find_1d_coeffs_in_patch(f,0);
      /* f_tmp1->v == f->v2 
      // f_tmp1->v2 == 0
      */
      f_tmp1->v2 = alloc_double(N);
      f_tmp2->v = find_1d_coeffs_in_patch(f_tmp1,1);
      /* f_tmp2->v == f_tmp1->v2
      // f_tmp2->v2 == 0
      */
      free_v2(f);/* => free(f_tmp1->v)*/
      f_tmp2->v2 = alloc_double(N);
      f->v2      = find_1d_coeffs_in_patch(f_tmp2,2);
      /* f_tmp2->v == f_tmp1->v2
      // f_tmp2->v2 == f->v2
      */
      free_v2(f_tmp1);
      f_tmp1->v = 0;
      f_tmp2->v = 0;
      f_tmp2->v2 = 0;
      remove_field(f_tmp1);
      remove_field(f_tmp2);
    }
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
  char needle[3][MAX_STR] = {'\0'};
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
  char needle[6][MAX_STR] = {'\0'};
  Flag_T flg = CLEAN;
  
  *dir = INT_MAX;
  collocation[dir1] = patch->collocation[dir1];
  collocation[dir2] = patch->collocation[dir2];
  basis[dir1]       = patch->basis[dir1];
  basis[dir2]       = patch->basis[dir2];
  
  sprintf(needle[0],PR_FORMAT,pn,dir1,collocation[dir1],basis[dir1]);
  sprintf(needle[1],"p%u,d%u",pn,(dir1+1)%3);
  sprintf(needle[2],"p%u,d%u",pn,3-dir1-(dir1+1)%3);
  sprintf(needle[3],PR_FORMAT,pn,dir2,collocation[dir2],basis[dir2]);
  sprintf(needle[4],"p%u,d%u",pn,(dir2+1)%3);
  sprintf(needle[5],"p%u,d%u",pn,3-dir2-(dir2+1)%3);
  
  /* if there is no coeffs */
  if (!f->v2 || !f->info)
  {
    flg = CLEAN;
    r = 0;
  }
  /* if other direction exists one needs to clean coeffs */
  else if ( strstr(f->info,needle[1]) || strstr(f->info,needle[2]) || 
       strstr(f->info,needle[4]) || strstr(f->info,needle[5])   )
  {
    flg = CLEAN;
    r = 0;
  }
  /* if coeffs is ready */
  else if (strstr(f->info,needle[0]) && strstr(f->info,needle[3]))
  {
    flg = NO;
    r = 1;
    *dir = -1;
  }
  else if (strstr(f->info,needle[0]))
  {
    flg = NO;
    r = 1;
    *dir = (int)dir2;
      
  }
  else if (strstr(f->info,needle[3]))
  {
    flg = NO;
    r = 1;
    *dir = (int)dir1;
      
  }
  else
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
// ->return value: 1 if already made, 0 otherwise.
*/
static unsigned IsAvailable_3d(Field_T *const f)
{
  Patch_T *const patch = f->patch;
  Collocation_T collocation[3];
  Basis_T basis[3];
  unsigned r = 0;
  unsigned pn = patch->pn;
  char needle[3][MAX_STR] = {'\0'};
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
      add_Tinfo(f,dir,Chebyshev_Tn_BASIS,Chebyshev_Extrema);
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
  double *const coeffs = f->v2;
  const double *const values = f->v;
  double *out, *in;
  const unsigned *n = f->patch->n;
  const unsigned B = n[(dir+1)%3]*n[3-dir-(dir+1)%3];
  unsigned l,i,j,k,s;
  out = alloc_double(n[dir]);
  in = alloc_double(n[dir]);
  
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
      fftw_1d_ChebyshevExtrema_coeffs(in,out,n[dir]);
      
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
      fftw_1d_ChebyshevExtrema_coeffs(in,out,n[dir]);
      
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
      fftw_1d_ChebyshevExtrema_coeffs(in,out,n[dir]);
      
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

/* if one can performe 3d fft */
static unsigned Is3d_fft(const Collocation_T *collocation,const Basis_T *basis)
{
  unsigned r = 0;
  
  /* if all bases are the same type */
  if (basis[0] == basis[1] && basis[1] == basis[2])
  {
    /* for Chebyshev Tn */
    if (basis[0] == Chebyshev_Tn_BASIS)
    {
      /* if all collocations are the same type */
      if (collocation[0] == collocation[1] && 
          collocation[1] == collocation[2])
      {
        /* for Chebyshev extrema */
        if (collocation[0] == Chebyshev_Extrema)
          r = 1;
        else
          r = 0;
      }
    }
    else
      abortEr("No such a basis defined for this function.\n");
  }
  else
    r = 0;
  
  return r;
}
