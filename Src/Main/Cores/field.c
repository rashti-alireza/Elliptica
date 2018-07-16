/*
// Alireza Rashti
// July 2018
*/

#include "field.h"
#define MAX_STR 400

/* initializing a 3d field.
// ->return value: a ready field.
*/
Field_T *init_field_3d(const char *const name,Grid_T *const grid)
{
  Field_T *f = calloc(1,sizeof(*f));
  if (name)
    f->name = dup_s(name);
  f->grid = grid;
  f->dim = 3;
  
  return f;
}

/* initializing a 2d field.
// ->return value: a ready field.
*/
Field_T *init_field_2d(const char *const name,Grid_T *const grid)
{
  Field_T *f = calloc(1,sizeof(*f));
  if (name)
    f->name = dup_s(name);
  f->grid = grid;
  f->dim = 2;
  
  return f;
}

/* adding a field to field data base */
void add_field(Field_T *const f,Grid_T *const grid)
{
  grid->field = realloc(grid->field,(grid->nf+1)*sizeof(*grid->field));
  pointerEr(grid->field);
  grid->field[grid->nf] = f;
  grid->nf++;
}

/* getting a field when its name is given (String version) 
// return value-> pointer to found field, 0 otherwise
*/
Field_T *get_field_S(const char *const name,Grid_T *const grid)
{
  unsigned f;
  
  for (f = 0; f < grid->nf; ++f)
    if (strcmp_i(grid->field[f]->name,name))
      return grid->field[f];
      
  return 0;
}

/* values -> 1-d coeffs in specified direction and patch.
// making coeffs of a field based on basis used for expansion.
// arguments of function:
// ======================
//
// o. f is demanded field
// o. dir referes to directions which this coefficents may be found.
//	0 for a-direction, 1 for b-direction, 2 for c-direction.
// o. patch: coeffs will be made in this patch only.
//
// note1: each patch uses it own basis type and each part of coeffs
// corresponds to that basis type.
// note3: it allocates memory for coeffs if not available.
// ->return value: coeffs.
*/
double *make_coeffs_1d(Field_T *const f,const unsigned dir,const Patch_T *const patch)
{
  const unsigned N = total_nodes_grid(f->grid);
  
  assert(dir <= 2);
  assert(patch);
  assert(f->dim == 3);
  
  if (IsAvailable_1d(f,patch,dir))
    return f->coeffs;
  else
  {
    if (!f->coeffs)
      f->coeffs = alloc_double(N);
    
    find_1d_coeffs_in_patch(f,patch,dir);
  }
  
  return f->coeffs;
}

/* values -> 2-d coeffs in specified direction and patch.
// for more info refere to make_coeffs_1d
// ->return value: coeffs.
*/
double *make_coeffs_2d(Field_T *const f,const unsigned dir1,const unsigned dir2,const Patch_T *const patch)
{
  const unsigned N = total_nodes_grid(f->grid);
  int dir;
  Field_T *f_tmp = init_field_3d("f_tmp",f->grid);
  
  assert(dir2 <= 2);
  assert(dir1 <= 2);
  assert(patch);
  assert(f->dim == 3);
  
  if (IsAvailable_2d(f,patch,dir1,dir2,&dir))
  {
    /* if field is ready ready */
    if (dir == -1)
      return f->coeffs;
      
    f_tmp->values = f->coeffs;
    f->coeffs = find_1d_coeffs_in_patch(f_tmp,patch,(unsigned)dir);
    
  }
  else
  {
    if (!f->coeffs)
      f->coeffs = alloc_double(N);
    
    find_1d_coeffs_in_patch(f,patch,dir1);
    f_tmp->values = f->coeffs;
    f->coeffs = find_1d_coeffs_in_patch(f_tmp,patch,dir2);
    find_1d_coeffs_in_patch(f,patch,(unsigned)dir);
  }
  
  f_tmp->values = 0;
  f_tmp->coeffs = 0;
  free_field(f_tmp);
    
  return f->coeffs;
}

/* values -> 3-d coeffs in specified patch.
// for more info refere to make_coeffs_1d
// ->return value: coeffs.
*/
double *make_coeffs_3d(Field_T *const f,const Patch_T *const patch)
{
  Grid_T *const grid = f->grid;
  const unsigned N = total_nodes_grid(grid);
  Collocation_T collocation[3];
  Basis_T basis[3];
  
  collocation[0] = patch->collocation[0];
  collocation[1] = patch->collocation[1];
  collocation[2] = patch->collocation[2];
  basis[0]       = patch->basis[0];
  basis[1]       = patch->basis[1];
  basis[2]       = patch->basis[2];
  
  assert(patch);
  assert(f->dim == 3);
  
  if (IsAvailable_3d(f,patch))
  {
    return f->coeffs;
  }
  else
  {
    if (!f->coeffs)
      f->coeffs = alloc_double(N);
    
    /* when all of directions and basis are the same.
    // it gets used for high performance.
    */
    if (Is3d_fft(collocation,basis))
    {
      double *coeffs = &f->coeffs[patch->nc];
      double *values = &f->values[patch->nc];
      const unsigned *n = patch->n;
    
      if (patch->basis[0] == Chebyshev_Tn_BASIS)
        fftw_3d_ChebyshevExtrema_coeffs(values,coeffs,n);
      else
        abortEr("No such basis is defined for this function.\n");
    }
    else
    {
      Field_T *f_tmp1 = init_field_3d("f_tmp1",f->grid);
      Field_T *f_tmp2 = init_field_3d("f_tmp2",f->grid);
 
      f_tmp1->values = find_1d_coeffs_in_patch(f,patch,0);
      /* f_tmp1->values == f->coeffs 
      // f_tmp1->coeffs == 0
      */
      f_tmp1->coeffs = alloc_double(N);
      f_tmp2->values = find_1d_coeffs_in_patch(f_tmp1,patch,1);
      /* f_tmp2->values == f_tmp1->coeffs
      // f_tmp2->coeffs == 0
      */
      f_tmp2->coeffs = alloc_double(N);
      free(f->coeffs);/* => free(f_tmp1->values)*/
      f->coeffs      = find_1d_coeffs_in_patch(f_tmp2,patch,2);
      /* f_tmp2->values == f_tmp1->coeffs
      // f_tmp2->coeffs == f->coeffs
      */
      free(f_tmp1->coeffs);
      f_tmp1->values = 0;
      f_tmp1->coeffs = 0;
      f_tmp2->values = 0;
      f_tmp2->coeffs = 0;
      free_field(f_tmp1);
      free_field(f_tmp2);
    }
  }
  
  return f->coeffs;
}

/* values -> 1-d coeffs
// making coeffs of a field based on basis used for expansion
// on the whole grid for the specified direction.
// ->return value: coeffs.
*/
double *make_coeffs_grid_1d(Field_T *const f,const unsigned dir)
{
  Grid_T *const grid = f->grid;
  double *coeffs = 0;
  unsigned pa;
  
  FOR_ALL(pa,grid->patch)
  {
    coeffs = make_coeffs_1d(f,dir,grid->patch[pa]);
  }

  return coeffs;
}  

/* values -> 2-d coeffs
// making coeffs of a field based on basis used for expansion
// on the whole grid for the specified directions.
// ->return value: coeffs.
*/
double *make_coeffs_grid_2d(Field_T *const f,const unsigned dir1,const unsigned dir2)
{
  Grid_T *const grid = f->grid;
  double *coeffs = 0;
  unsigned pa;
  
  FOR_ALL(pa,grid->patch)
  {
    coeffs = make_coeffs_2d(f,dir1,dir2,grid->patch[pa]);
  }

  return coeffs;
}  


/* values -> 3-d coeffs.
// making coeffs of a field based on basis used for expansion
// on the whole grid.
// ->return value: coeffs.
*/
double *make_coeffs_grid_3d(Field_T *const f)
{
  Grid_T *const grid = f->grid;
  double *coeffs = 0;
  unsigned pa;
  
  FOR_ALL(pa,grid->patch)
  {
    coeffs = make_coeffs_3d(f,grid->patch[pa]);
  }

  return coeffs;
}

/* check if this coeffs has been already made.
// ->return value: 1 if already made, 0 otherwise.
*/
static unsigned IsAvailable_1d(Field_T *const f,const Patch_T *const patch,const unsigned dir)
{
  Collocation_T collocation = patch->collocation[dir];
  Basis_T basis = patch->basis[dir];
  unsigned r = 0;
  unsigned pn = patch->pn;
  char needle[3][MAX_STR] = {'\0'};
  
  /* if there is no coeffs */
  if (!f->coeffs) return 0;
  
  sprintf(needle[0],"p%u,d%u",pn,dir);
  sprintf(needle[1],"p%u,d%u",pn,(dir+1)%3);
  sprintf(needle[2],"p%u,d%u",pn,3-dir-(dir+1)%3);
  
  /* if other direction exist one needs to clean coeffs */
  if (strstr(f->info,needle[1]) || strstr(f->info,needle[2]))
  {
    free(f->coeffs);
    f->coeffs = 0;
    r = 0;
  }
  /* if coeffs is ready */
  else if (strstr(f->info,needle[0]))
  {
    /* check if the basis and collocation are also match */
    sprintf(needle[0],"p%u,d%u,c%d,b%d",pn,dir,collocation,basis);
    if (strstr(f->info,needle[0]))
      r = 1;
    else
    {
      free(f->coeffs);
      f->coeffs = 0;
      r = 0;
    }
  }
  
  return r;
}

/* check if this coeffs has been already made.
// if only one of them is available put the other one in dir.
// ->return value: 1 if already made, 0 otherwise.
// ->return value: if 1 and dir == -1, it means the field is ready ready.
*/
static unsigned IsAvailable_2d(Field_T *const f,const Patch_T *const patch,const unsigned dir1,const unsigned dir2,int *dir)
{
  Collocation_T collocation[2];
  Basis_T basis[2];
  unsigned r = 0;
  unsigned pn = patch->pn;
  char needle[6][MAX_STR] = {'\0'};
  
  collocation[0] = patch->collocation[dir1];
  collocation[1] = patch->collocation[dir2];
  basis[0]       = patch->basis[dir1];
  basis[1]       = patch->basis[dir2];
  
  /* if there is no coeffs */
  if (!f->coeffs) return 0;
  
  sprintf(needle[0],"p%u,d%u",pn,dir1);
  sprintf(needle[1],"p%u,d%u",pn,(dir1+1)%3);
  sprintf(needle[2],"p%u,d%u",pn,3-dir1-(dir1+1)%3);
  sprintf(needle[3],"p%u,d%u",pn,dir2);
  sprintf(needle[4],"p%u,d%u",pn,(dir2+1)%3);
  sprintf(needle[5],"p%u,d%u",pn,3-dir2-(dir2+1)%3);
  
  /* if other direction exists one needs to clean coeffs */
  if ( strstr(f->info,needle[1]) || strstr(f->info,needle[2]) || 
       strstr(f->info,needle[3]) || strstr(f->info,needle[4])   )
  {
    free(f->coeffs);
    f->coeffs = 0;
    r = 0;
  }
  /* if coeffs is ready */
  else if (strstr(f->info,needle[0]) && strstr(f->info,needle[3]))
  {
    /* check if the basis and collocation are also match */
    sprintf(needle[0],"p%u,d%u,c%d,b%d",pn,dir1,collocation[0],basis[0]);
    sprintf(needle[3],"p%u,d%u,c%d,b%d",pn,dir2,collocation[1],basis[1]);
    if (strstr(f->info,needle[0]) && strstr(f->info,needle[3]))
    {
      r = 1;
      *dir = -1;
    }
    else
    {
      free(f->coeffs);
      f->coeffs = 0;
      r = 0;
    }
  }
  else if (strstr(f->info,needle[0]))
  {
    /* check if the basis and collocation are also match */
    sprintf(needle[0],"p%u,d%u,c%d,b%d",pn,dir1,collocation[0],basis[0]);
    if (strstr(f->info,needle[0]))
    {
      r = 1;
      *dir = (int)dir2;
    }
    else
    {
      free(f->coeffs);
      f->coeffs = 0;
      r = 0;
    }
  }
  else if (strstr(f->info,needle[3]))
  {
    /* check if the basis and collocation are also match */
    sprintf(needle[3],"p%u,d%u,c%d,b%d",pn,dir2,collocation[1],basis[1]);
    if (strstr(f->info,needle[3]))
    {
      r = 1;
      *dir = (int)dir1;
    }
    else
    {
      free(f->coeffs);
      f->coeffs = 0;
      r = 0;
    }
  }
  
  return r;
}

/* check if this coeffs has been already made.
// ->return value: 1 if already made, 0 otherwise.
*/
static unsigned IsAvailable_3d(Field_T *const f,const Patch_T *const patch)
{
  Collocation_T collocation[3];
  Basis_T basis[3];
  unsigned r = 0;
  unsigned pn = patch->pn;
  char needle[3][MAX_STR] = {'\0'};
  
  collocation[0] = patch->collocation[0];
  collocation[1] = patch->collocation[1];
  collocation[2] = patch->collocation[2];
  basis[0]       = patch->basis[0];
  basis[1]       = patch->basis[1];
  basis[2]       = patch->basis[2];
  
  /* if there is no coeffs */
  if (!f->coeffs) return 0;
  
  sprintf(needle[0],"p%u,d%u,c%d,b%d",pn,0,collocation[0],basis[0]);
  sprintf(needle[1],"p%u,d%u,c%d,b%d",pn,1,collocation[1],basis[1]);
  sprintf(needle[2],"p%u,d%u,c%d,b%d",pn,2,collocation[2],basis[2]);
  
  /* if other direction exist one needs to clean coeffs */
  if (strstr(f->info,needle[0]) && strstr(f->info,needle[1]) 
      && strstr(f->info,needle[2]))
  {
    r = 1;
  }
  else
  {
    free(f->coeffs);
    f->coeffs = 0;
    r = 0;
  }
  
  return r;
}

/* finding NORMALIZED coeffs of expansion in direction dir 
// for a given field, patch and direction.
// ->return value: coeffs
*/
static double *find_1d_coeffs_in_patch(Field_T *const f,const Patch_T *const patch,const unsigned dir)
{
  const Collocation_T collocation = patch->collocation[dir];
  const Basis_T basis = patch->basis[dir];
  double *const coeffs = &f->coeffs[patch->nc];
  
  if(basis == Chebyshev_Tn_BASIS)
  {
    if (collocation == Chebyshev_Extrema)
    {
      coeffs_patch_Tn_Extrema_1d(f,patch,dir,coeffs);
      add_Tinfo(f,patch->pn,dir,Chebyshev_Tn_BASIS,Chebyshev_Extrema);
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
    
  return coeffs;
}

/* finding coeffs in in patch For Tn basis 
// with Chebyshev extrema collocation.
*/
static void coeffs_patch_Tn_Extrema_1d(const Field_T *const f,const Patch_T *const patch,const unsigned dir,double *const coeffs)
{
  double *coeffs_tmp, *values;
  const unsigned *n = patch->n;
  const unsigned B = n[(dir+1)%3]*n[3-dir-(dir+1)%3];
  const double Nrm = 2*(n[dir]-1);
  unsigned l,i,j,k,s;
  coeffs_tmp = alloc_double(n[dir]);
  values = alloc_double(n[dir]);
  
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
        values[i] = f->values[l];
      }
      fftw_1d_ChebyshevExtrema_coeffs(values,coeffs_tmp,n[dir]);
      
      for (i = 0; i < n[dir]; ++i)
      {
        l = L(n,i,j,k);
        coeffs[l] = coeffs_tmp[i]/Nrm;
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
        values[j] = f->values[l];
      }
      fftw_1d_ChebyshevExtrema_coeffs(values,coeffs_tmp,n[dir]);
      
      for (j = 0; j < n[dir]; ++j)
      {
        l = L(n,i,j,k);
        coeffs[l] = coeffs_tmp[j]/Nrm;
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
        values[k] = f->values[l];
      }
      fftw_1d_ChebyshevExtrema_coeffs(values,coeffs_tmp,n[dir]);
      
      for (k = 0; k < n[dir]; ++k)
      {
        l = L(n,i,j,k);
        coeffs[l] = coeffs_tmp[k]/Nrm;
      }
    }
  }/* end of if (dir == 2) */
  else
    abortEr("direction must be 0,1 or 2.\n");
    
  free(coeffs_tmp);
  free(values);
}

/* adding Transformation info in f->info.
// for each transformation we have (patch,direction,collocation,basis)
// in row.
*/
static void add_Tinfo(Field_T *const f,const unsigned pn,const unsigned dir,const Collocation_T collocation,const Basis_T basis)
{
  char inf[MAX_STR] = {'\0'};
  
  /* (patch,direction,collocation,basis) */
  sprintf(inf,"(p%u,d%u,c%d,b%d)",pn,dir,collocation,basis);
  f->info = realloc(f->info,strlen(f->info)+strlen(inf)+1);
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
    else if (basis[0] == No_BASIS)
    {
      r = 0;
    }
    else
      abortEr("No such a basis defined for this function.\n");
  }
  else
    r = 0;
  
  return r;
}
