/*
// Alireza Rashti
// July 2018
*/

#include "field.h"

/* adding a field to field data base */
void add_field(const char *const name,Grid_T *const grid)
{
  grid->field = realloc(grid->field,(grid->nf+1)*sizeof(*grid->field));
  pointerEr(grid->field);
  grid->field[grid->nf] = alloc_field(grid);
  grid->field[grid->nf]->name = dup_s(name);
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
