/*
// Alireza Rashti
// July 2018
*/

#include "field.h"

/* initializing a field.
// ->return value: a ready field.
*/
Field_T *init_field(const char *const name,Grid_T *const grid)
{
  Field_T *f = alloc_field(grid);
  if (name)
    f->name = dup_s(name);
  f->grid = grid;
  
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

/* values -> coeffs.
// making coeffs of a field based on basis used for expansion.
// note: each patch uses it own basis type and each part of coeffs
// corresponds to that basis type.
// ->return value: coeffs, null no patches use basis.
*/
double *make_coeffs(Field_T *const f)
{
  Grid_T *const grid = f->grid;
  double *coeffs, *values;
  unsigned i,*n;
  unsigned pa;
  
  i = 0;
  FOR_ALL(pa,grid->patch)
  {
    coeffs = &f->coeffs[i];
    values = &f->values[i];
    n = grid->patch[pa]->n;
    
    /* basis finders come here */
    if (grid->patch[pa]->basis == Chebyshev_FirstKind_BASIS)
      fftw_3d_ChebyshevExtrema_coeffs(values,coeffs,(int *)n);
      
    i += total_nodes_patch(grid->patch[pa]);
  }
  
  return f->coeffs;
}

/* coeffs -> values
// making values of field based on basis used for expansion.
// note: each patch uses it own basis type and each part of coeffs
// corresponds to that basis type.
// ->return value: values.
*/
double *make_coeffs_inverse(Field_T *const f)
{
  Grid_T *const grid = f->grid;
  double *coeffs, *values;
  unsigned i,*n;
  unsigned pa;
  
  i = 0;
  FOR_ALL(pa,grid->patch)
  {
    coeffs = &f->coeffs[i];
    values = &f->values[i];
    n = grid->patch[pa]->n;
    
    /* basis finders come here */
    if (grid->patch[pa]->basis == Chebyshev_FirstKind_BASIS)
      fftw_3d_ChebyshevExtrema_values(values,coeffs,(int *)n);
      
    i += total_nodes_patch(grid->patch[pa]);
  }
  
  return f->values;
}
