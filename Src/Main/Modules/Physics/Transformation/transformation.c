/*
// Alireza Rashti
// October 2019
*/

#include "transformation.h"

/* allocating memory and initializing Transformation_T structure
// ->return value: initialized structure. */
Transformation_T *initialize_transformation(void)
{
  Transformation_T *t = calloc(1,sizeof(*t));
  pointerEr(t);
  
  return t;
}

/* free transformation */
void free_transformation(Transformation_T *t)
{
  _free(t);
}
