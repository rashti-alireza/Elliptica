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

//make_coeffs(Field_T *const f)