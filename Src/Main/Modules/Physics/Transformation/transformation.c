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
// ** B O O S T  **
// ** populating **
// t->boost->Bx = Bx;
// t->boost->By = By;
// t->boost->Bz = Bz;
// t->boost->B2 = Pow2(Bx)+Pow2(By)+Pow2(Bz); # assumed Minkowski space-time
//
// ** if you wanna inverse transformation **
// t->boost->inverse = 1;
//
// ** transforming four vector u1 to u2 **
// boost_transformation(t,u1,u2);
//
// ** R O T A T I O N **
// ** populating **
// T->rotation->Rx = phi0; => phi0 along x-axis
// T->rotation->Ry = phi1; => phi1 along y-axis
// T->rotation->Rz = phi2; => phi2 along z-axis
//
// ** transforming four vector u1 to u2 **
// ** u2 = Rz(phi2)*Ry(phi1)*Rx(phi0) u1 **
// ** note: the oreder is Rz(phi2)*Ry(phi1)*Rx(phi0) **
// rotation_transformation(t,u1,u2);
//
// ** freeing **
// free_transformation(t);
*/

/* allocating memory and initializing Transformation_T structure
// ->return value: initialized structure. */
Transformation_T *initialize_transformation(void)
{
  Transformation_T *t = calloc(1,sizeof(*t));
  IsNull(t);
  
  return t;
}

/* free transformation */
void free_transformation(Transformation_T *t)
{
  _free(t);
}
