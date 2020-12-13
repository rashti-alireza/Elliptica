#ifndef physics_transformation_LIB_H
#define physics_transformation_LIB_H
#include "elliptica_system_lib.h"


/* struct for transformation boost, rotation etc. */
typedef struct TRANSFORMATION_T
{
  /* Lorentz Boost */
  struct
  {
    double Bx;/* v_x/c in vector B */
    double By;/* v_y/c in vector B*/
    double Bz;/* v_z/c in vector B*/
    double B2;/* B.B, this needs to be given by user b/c of dot product */
    int inverse;/* if 1 it uses the INVERT of the Lorentz matrix transformation, 
                // 0 uses the Lorentz matrix transformation*/
  }boost[1];
  
  /* rotation */
  struct
  {
    double Rx;/* Rx amount along x-axis */
    double Ry;/* Ry amount along y-axis */
    double Rz;/* Rz amount along z-axis */
    int active;/* if 1 means active transformation, otherwise means
               // passive transformation */
  }rotation[1];
  
  /* spherical to Cartesian coords and reverse */
  struct
  {
    double r;/* r = (x^2+y^2+z^2)^0.5 */
    double th;/* theta */
    double ph;/* phi */
    Uint s2c:1;/* 1 means spherical to cartesian */
    Uint c2s:1;/* 1 means cartesian to spherical */
  }spheNcart[1];
  
}Transformation_T;

void rotation_transformation(Transformation_T *const t,const double *const in,double *const out);
void boost_transformation(Transformation_T *const t,const double *const in,double *const out);
void spheNcart_transformation(Transformation_T *const t,const double *const in,double *const out);
Transformation_T *init_transformation(void);
void free_transformation(Transformation_T *t);


#endif


