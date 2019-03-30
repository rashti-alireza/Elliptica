/*
// Alireza Rashti
// March 2018
*/

#include "projective_coordinate.h"

/* making value of coords. it is a general function for Protective Hemisphere Up type */
void make_nodes_ProjectiveHemisphereUp_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const unsigned U = patch->nn;
  const unsigned *const n = patch->n;
  const Field_T *const R1_field = patch->pool[Ind("R1_ProjectiveHemisphereUp")];
  const Field_T *const R2_field = patch->pool[Ind("R2_ProjectiveHemisphereUp")];
  double R1,R2;
  const double *const c = patch->c;/* center of origine translated */
  double r2_x2_y2;
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
    X[0] = point_value(i,&coll_s[0]);
    X[1] = point_value(j,&coll_s[1]);
    X[2] = point_value(k,&coll_s[2]);
    patch->node[l]->X = X;
    
    R1 = R1_field->v[L(n,i,j,0)];
    R2 = R2_field->v[L(n,i,j,0)];
    r = 0.5*X[2]*(R2-R1)+0.5*(R2+R1);
    
    x[0] = r*X[0]*sqrt(1-0.5*SQR(X[1])); assert(!isnan(x[0]));
    x[1] = r*X[1]*sqrt(1-0.5*SQR(X[0])); assert(!isnan(x[1]));
    r2_x2_y2 = SQR(r)-SQR(x[0])-SQR(x[1]);
    if (EQL(r2_x2_y2,0))
      r2_x2_y2 = 0;/* avoiding infinitesimal negative */
    x[2] = sqrt(r2_x2_y2); assert(!isnan(x[2]));
    
    x[0]+= c[0];
    x[1]+= c[1];
    x[2]+= c[2];
    
  }
}

/* making value of coords. it is a general function for Protective Hemisphere Down type */
void make_nodes_ProjectiveHemisphereDown_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const unsigned U = patch->nn;
  const unsigned *const n = patch->n;
  const Field_T *const R1_field = patch->pool[Ind("R1_ProjectiveHemisphereDown")];
  const Field_T *const R2_field = patch->pool[Ind("R2_ProjectiveHemisphereDown")];
  double R1,R2;
  double r2_x2_y2;
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
    X[0] = point_value(i,&coll_s[0]);
    X[1] = point_value(j,&coll_s[1]);
    X[2] = point_value(k,&coll_s[2]);
    patch->node[l]->X = X;

    R1 = R1_field->v[L(n,i,j,0)];
    R2 = R2_field->v[L(n,i,j,0)];
    r = 0.5*X[2]*(R2-R1)+0.5*(R2+R1);    
    
    x[0] = r*X[0]*sqrt(1-0.5*SQR(X[1])); assert(!isnan(x[0]));
    x[1] = r*X[1]*sqrt(1-0.5*SQR(X[0])); assert(!isnan(x[1]));
    r2_x2_y2 = SQR(r)-SQR(x[0])-SQR(x[1]);
    if (EQL(r2_x2_y2,0))
      r2_x2_y2 = 0;/* avoiding infinitesimal negative */
    x[2] = -sqrt(r2_x2_y2); assert(!isnan(x[2]));
    
    x[0]+= c[0];
    x[1]+= c[1];
    x[2]+= c[2];
  }
}

/* making value of coords. it is a general function for Stereographic Sphere Left type
// projected in y = 0 plane for a sphere at center (0,-R0,0) */
void make_nodes_StereographicSphereLeft_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const unsigned U = patch->nn;
  const double R0 = fabs(patch->c[1]);
  const double R1 = patch->CoordSysInfo->R1;
  const double R2 = patch->CoordSysInfo->R2;
  const unsigned *const n = patch->n;
  unsigned i,j,k,l;
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  initialize_collocation_struct(patch,&coll_s[2],2);
  
  for (l = 0; l < U; l++)
  {
    double *X = alloc_double(3);
    double *x = patch->node[l]->x;
    double u,w;
    double r,R,A,c;
    
    IJK(l,n,&i,&j,&k);
    X[0] = point_value(i,&coll_s[0]);
    X[1] = point_value(j,&coll_s[1]);
    X[2] = point_value(k,&coll_s[2]);
    patch->node[l]->X = X;
    
    r = 0.5*X[1]*(R2-R1)+0.5*(R2+R1);
    R = sqrt(SQR(r)-SQR(R0)); assert(!isnan(R));
    u = R*X[0]*sqrt(1-0.5*SQR(X[2])); assert(!isnan(u));
    w = R*X[2]*sqrt(1-0.5*SQR(X[0])); assert(!isnan(w));
    A = SQR(u/(R0-r))+SQR(w/(R0-r))+1;
    x[1] = -r*(2/A-1)-R0;
    c = 2*r/(A*(r-R0));
    x[0] = c*u;
    x[2] = c*w;
    
  }
}

/* making value of coords. it is a general function for Stereographic Sphere Right type 
// projected in y = 0 plane for a sphere at center (0,R0,0) */
void make_nodes_StereographicSphereRight_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const unsigned U = patch->nn;
  const unsigned *const n = patch->n;
  const double R0 = patch->c[1];
  const double R1 = patch->CoordSysInfo->R1;
  const double R2 = patch->CoordSysInfo->R2;
  unsigned i,j,k,l;
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  initialize_collocation_struct(patch,&coll_s[2],2);
  
  for (l = 0; l < U; l++)
  {
    double *X = alloc_double(3);
    double *x = patch->node[l]->x;
    double u,w;
    double r,R,A,c;
    
    IJK(l,n,&i,&j,&k);
    X[0] = point_value(i,&coll_s[0]);
    X[1] = point_value(j,&coll_s[1]);
    X[2] = point_value(k,&coll_s[2]);
    patch->node[l]->X = X;
    
    r = 0.5*X[1]*(R2-R1)+0.5*(R2+R1);
    R = sqrt(SQR(r)-SQR(R0)); assert(!isnan(R));
    u = R*X[0]*sqrt(1-0.5*SQR(X[2])); assert(!isnan(u));
    w = R*X[2]*sqrt(1-0.5*SQR(X[0])); assert(!isnan(w));
    A = SQR(u/(R0-r))+SQR(w/(R0-r))+1;
    x[1] = r*(2/A-1)+R0;
    c = 2*r/(A*(r-R0));
    x[0] = c*u;
    x[2] = c*w;
    
  }
}
