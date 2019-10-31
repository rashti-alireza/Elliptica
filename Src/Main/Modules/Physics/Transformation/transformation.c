/*
// Alireza Rashti
// October 2019
*/

#include "transformation.h"

/* synopsis:
// =========
// 
// ** initializing **
// Transformation_T *t = initialize_transformation();
//
// ** populating **
// t->boost->Bx = Bx;
// t->boost->By = By;
// t->boost->Bz = Bz;
// t->boost->B2 = SQR(Bx)+SQR(By)+SQR(Bz); # assumed Minkowski space-time
//
// ** if you wanna inverse transformation **
// t->boost->inverse = 1;
//
// ** transforming four vector u1 to u2 **
// Lorentz_boost(t,u1,u2);
//
// ** freeing **
// free_transformation(t);
*/

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
