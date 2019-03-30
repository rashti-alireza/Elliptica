/*
// Alireza Rashti
// March 2018
*/

#include "spherical_coordinate.h"

/* making value of coords. it is a general function for Spherical type */
void make_nodes_Spherical_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const unsigned U = patch->nn;
  const unsigned *const n = patch->n;
  const Field_T *const R1_field = patch->pool[Ind("R1_radius")];
  const Field_T *const R2_field = patch->pool[Ind("R2_radius")];
  double R1,R2;
  const double *const c = patch->c;/* center of origine translated */
  unsigned i,j,k,l;
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  initialize_collocation_struct(patch,&coll_s[2],2);
  
  for (l = 0; l < U; l++)
  {
    double *X = alloc_double(3);
    double *x = patch->node[l]->x;
    double r;
    
    IJK(l,n,&i,&j,&k);
    X[0] = point_value(i,&coll_s[0]);/* r */
    X[1] = point_value(j,&coll_s[1]);/* theta */
    X[2] = point_value(k,&coll_s[2]);/* phi */
    patch->node[l]->X = X;
    
    R1 = R1_field->v[L(n,0,j,k)];
    R2 = R2_field->v[L(n,0,j,k)];
    r = X[0]*(R2-R1)+R1;
    
    x[0] = r*sin(X[1])*cos(X[2]);
    x[1] = r*sin(X[1])*sin(X[2]);
    x[2] = r*cos(X[1]);
    
    x[0]+= c[0];
    x[1]+= c[1];
    x[2]+= c[2];
    
  }
}
