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
  const Field_T *const R1_field = patch->pool[Ind("R1_ProjectiveHemisphere")];
  const Field_T *const R2_field = patch->pool[Ind("R2_ProjectiveHemisphere")];
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
  const Field_T *const R1_field = patch->pool[Ind("R1_ProjectiveHemisphere")];
  const Field_T *const R2_field = patch->pool[Ind("R2_ProjectiveHemisphere")];
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

/* making Jacobian transformation for ProjectiveHemisphere coord. */
void make_JacobianT_ProjectiveHemisphere_coord(Patch_T *const patch)
{
  patch->JacobianT->j      = JT_ProjectiveHemisphere;
  patch->JacobianT->dN0_dx = dN0_dx_ProjectiveHemisphere_patch;
  patch->JacobianT->dN0_dy = dN0_dy_ProjectiveHemisphere_patch;
  patch->JacobianT->dN0_dz = dN0_dz_ProjectiveHemisphere_patch;
  patch->JacobianT->dN1_dx = dN1_dx_ProjectiveHemisphere_patch;
  patch->JacobianT->dN1_dy = dN1_dy_ProjectiveHemisphere_patch;
  patch->JacobianT->dN1_dz = dN1_dz_ProjectiveHemisphere_patch;
  patch->JacobianT->dN2_dx = dN2_dx_ProjectiveHemisphere_patch;
  patch->JacobianT->dN2_dy = dN2_dy_ProjectiveHemisphere_patch;
  patch->JacobianT->dN2_dz = dN2_dz_ProjectiveHemisphere_patch;
}

/* making Jacobian transformation for StereographicSphere Left coord. */
void make_JacobianT_StereographicSphereLeft_coord(Patch_T *const patch)
{
  patch->JacobianT->j      = JT_StereographicSphere_Left;
  patch->JacobianT->dN0_dx = dN0_dx_StereographicSphereLeft_patch;
  patch->JacobianT->dN0_dy = dN0_dy_StereographicSphereLeft_patch;
  patch->JacobianT->dN0_dz = dN0_dz_StereographicSphereLeft_patch;
  patch->JacobianT->dN1_dx = dN1_dx_StereographicSphereLeft_patch;
  patch->JacobianT->dN1_dy = dN1_dy_StereographicSphereLeft_patch;
  patch->JacobianT->dN1_dz = dN1_dz_StereographicSphereLeft_patch;
  patch->JacobianT->dN2_dx = dN2_dx_StereographicSphereLeft_patch;
  patch->JacobianT->dN2_dy = dN2_dy_StereographicSphereLeft_patch;
  patch->JacobianT->dN2_dz = dN2_dz_StereographicSphereLeft_patch;
}

/* making Jacobian transformation for StereographicSphere Right coord. */
void make_JacobianT_StereographicSphereRight_coord(Patch_T *const patch)
{
  patch->JacobianT->j      = JT_StereographicSphere_Right;
  patch->JacobianT->dN0_dx = dN0_dx_StereographicSphereRight_patch;
  patch->JacobianT->dN0_dy = dN0_dy_StereographicSphereRight_patch;
  patch->JacobianT->dN0_dz = dN0_dz_StereographicSphereRight_patch;
  patch->JacobianT->dN1_dx = dN1_dx_StereographicSphereRight_patch;
  patch->JacobianT->dN1_dy = dN1_dy_StereographicSphereRight_patch;
  patch->JacobianT->dN1_dz = dN1_dz_StereographicSphereRight_patch;
  patch->JacobianT->dN2_dx = dN2_dx_StereographicSphereRight_patch;
  patch->JacobianT->dN2_dy = dN2_dy_StereographicSphereRight_patch;
  patch->JacobianT->dN2_dz = dN2_dz_StereographicSphereRight_patch;
}

/* Calculating dN(1/2/3)/d(x/y/z) at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN(1/2/3)/d(x/y/z) */
static double dNi_dxj_ProjectiveHemisphere(Patch_T *const patch, const Dd_T Ni, const Dd_T xj,const double *const X)
{
  double dNi_dxj = 0;
  double J = 0;/* dX?/dx? */
  double a[3];/* Cartesian value of X */
  int rx_of_X = x_of_X(a,X,patch);/* return value of function */
  const double 
         x = a[0],
         y = a[1],
         z = a[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2);
  double R1 = 0;
  double R2 = 0;
  double dR2_dx = 0;
  double dR2_dy = 0;
  double dR2_dz = 0;
  double dR1_dx = 0;
  double dR1_dy = 0;
  double dR1_dz = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  Dd_T q2_e = UNDEFINED;
  
  if (rx_of_X == 0)
    abortEr("x-coordinate could not been found.");

  if (Ni == _N0_)
    q2_e = _a_;
  else if (Ni == _N1_)
    q2_e = _b_;
  else if (Ni == _N2_)
    q2_e = _c_;
  else
    abortEr(NO_OPTION);
  
  dA_da = get_dA_da(q2_e,xj);
  
  /* preparing R1 and R2 derivatives */
  if (q2_e == _c_)
  {
    R1_R2_derivative(patch);
    Field_T *const R1_field = patch->pool[Ind("R1_ProjectiveHemisphere")];
    Field_T *const R2_field = patch->pool[Ind("R2_ProjectiveHemisphere")];
    Field_T *const dR1_dx_field = patch->pool[Ind("dR1_dx")];
    Field_T *const dR1_dy_field = patch->pool[Ind("dR1_dy")];
    Field_T *const dR2_dx_field = patch->pool[Ind("dR2_dx")];
    Field_T *const dR2_dy_field = patch->pool[Ind("dR2_dy")];
    
    R1 	   = interpolation_2d_PH(R1_field,patch,X);
    dR1_dx = interpolation_2d_PH(dR1_dx_field,patch,X);
    dR1_dy = interpolation_2d_PH(dR1_dy_field,patch,X);
    dR1_dz = 0;
    R2 	   = interpolation_2d_PH(R2_field,patch,X);
    dR2_dx = interpolation_2d_PH(dR2_dx_field,patch,X);
    dR2_dy = interpolation_2d_PH(dR2_dy_field,patch,X);
    dR2_dz = 0;
  }
  
  switch(dA_da)
  {
    case da_dx:
      J = ((Sqrt(2)*Power(x,2)*(Power(y,2) + Power(z,2)) + 
          Sqrt(2)*Power(Power(y,2) + Power(z,2),2) - 
          r1*x*(2*Power(y,2) + Power(z,2)))/
          Sqrt((-2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 
          2*Power(z,2))/r2) + 
          (Sqrt(2)*Power(x,2)*(Power(y,2) + Power(z,2)) + 
          Sqrt(2)*Power(Power(y,2) + Power(z,2),2) + 
          r1*x*(2*Power(y,2) + Power(z,2)))/
          Sqrt((2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 
          2*Power(z,2))/r2))/(2.*Power(r1,5));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      J = ((2*y*(-(Sqrt(2)*r1*x) + 2*Power(x,2) + Power(z,2)))/(Power(r2,2)*Sqrt((-2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 2*Power(z,2))/r2)) - 
          (2*y*(r1*(2*Power(x,2) + Power(z,2)) + Sqrt(2)*x*(Power(x,2) + Power(y,2) + Power(z,2))))/
          (Power(r1,5)*Sqrt((2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 2*Power(z,2))/r2)))/4.;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      J = (z*((-(Sqrt(2)*r1*x) + Power(x,2) - Power(y,2))/Sqrt((-2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 2*Power(z,2))/r2) + 
          (-(Sqrt(2)*r1*x) - Power(x,2) + Power(y,2))/Sqrt((2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 2*Power(z,2))/r2)))/(2.*Power(r2,2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dx:
      J = (x*((-(Sqrt(2)*r1*y) + 2*Power(y,2) + Power(z,2))/
          Sqrt((Power(x,2) - 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2) - 
          (Sqrt(2)*r1*y + 2*Power(y,2) + Power(z,2))/
          Sqrt((Power(x,2) + 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2)))/
          (2.*Power(r2,2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dy:
      J = ((Power(x,2)*(Sqrt(2)*r1 - 2*y) + (Sqrt(2)*r1 - y)*Power(z,2))/
          Sqrt((Power(x,2) - 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2) + 
          (Power(x,2)*(Sqrt(2)*r1 + 2*y) + (Sqrt(2)*r1 + y)*Power(z,2))/
          Sqrt((Power(x,2) + 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2))/
          (2.*Power(r2,2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dz:
      J = (z*((-Power(x,2) + y*(-(Sqrt(2)*r1) + y))/
          Sqrt((Power(x,2) - 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2) + 
          (Power(x,2) - y*(Sqrt(2)*r1 + y))/
          Sqrt((Power(x,2) + 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2)))/
          (2.*Power(r2,2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dx:
      J = (2*x/r1 + dR1_dx-dR2_dx)/(2*(R1+R2))-
          (2*r1+R1-R2)*(dR1_dx+dR2_dx)/(2*SQR(R1+R2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      J = (2*y/r1 + dR1_dy-dR2_dy)/(2*(R1+R2))-
          (2*r1+R1-R2)*(dR1_dy+dR2_dy)/(2*SQR(R1+R2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      J = (2*z/r1 + dR1_dz-dR2_dz)/(2*(R1+R2))-
          (2*r1+R1-R2)*(dR1_dz+dR2_dz)/(2*SQR(R1+R2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
      
    break;
    default:
      abortEr("No such a enum!\n");
  }
  
  dNi_dxj = dN_dX(patch,Ni,q2_e)*J;
  
  return dNi_dxj;
}

/* Jacobian transformation for ProjectiveHemisphereUp patch.
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1
*/
double JT_ProjectiveHemisphere(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
  /* dx/d? not defined */
  else if (q2_e == _x_ || q2_e == _y_|| q2_e == _z_)
    abortEr(INCOMPLETE_FUNC);
  /* dx/d? not defined */
  else if (q1_e == _a_ || q1_e == _b_|| q1_e == _c_)
    abortEr(INCOMPLETE_FUNC);
  
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double 
         x = patch->node[p]->x[0],
         y = patch->node[p]->x[1],
         z = patch->node[p]->x[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2);
  const unsigned *const n = patch->n;
  double R1 = 0;
  double R2 = 0;
  double dR2_dx = 0;
  double dR2_dy = 0;
  double dR2_dz = 0;
  double dR1_dx = 0;
  double dR1_dy = 0;
  double dR1_dz = 0;
  
  /* if dZ_d? we need dR_d? thus: */
  if (q2_e == _c_)
  {
    unsigned i,j,k;
    IJK(p,n,&i,&j,&k);
    
    /* preparing R1 and R2 derivatives */
    R1_R2_derivative(patch);
    R1 = patch->pool[Ind("R1_ProjectiveHemisphere")]->v[L(n,i,j,0)];
    R2 = patch->pool[Ind("R2_ProjectiveHemisphere")]->v[L(n,i,j,0)];
    dR2_dx = patch->pool[Ind("dR2_dx")]->v[L(n,i,j,0)];
    dR2_dy = patch->pool[Ind("dR2_dy")]->v[L(n,i,j,0)];
    dR2_dz = 0;
    dR1_dx = patch->pool[Ind("dR1_dx")]->v[L(n,i,j,0)];
    dR1_dy = patch->pool[Ind("dR1_dy")]->v[L(n,i,j,0)];
    dR1_dz = 0;
  }
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      J = ((Sqrt(2)*Power(x,2)*(Power(y,2) + Power(z,2)) + 
          Sqrt(2)*Power(Power(y,2) + Power(z,2),2) - 
          r1*x*(2*Power(y,2) + Power(z,2)))/
          Sqrt((-2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 
          2*Power(z,2))/r2) + 
          (Sqrt(2)*Power(x,2)*(Power(y,2) + Power(z,2)) + 
          Sqrt(2)*Power(Power(y,2) + Power(z,2),2) + 
          r1*x*(2*Power(y,2) + Power(z,2)))/
          Sqrt((2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 
          2*Power(z,2))/r2))/(2.*Power(r1,5));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      J = ((2*y*(-(Sqrt(2)*r1*x) + 2*Power(x,2) + Power(z,2)))/(Power(r2,2)*Sqrt((-2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 2*Power(z,2))/r2)) - 
          (2*y*(r1*(2*Power(x,2) + Power(z,2)) + Sqrt(2)*x*(Power(x,2) + Power(y,2) + Power(z,2))))/
          (Power(r1,5)*Sqrt((2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 2*Power(z,2))/r2)))/4.;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      J = (z*((-(Sqrt(2)*r1*x) + Power(x,2) - Power(y,2))/Sqrt((-2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 2*Power(z,2))/r2) + 
          (-(Sqrt(2)*r1*x) - Power(x,2) + Power(y,2))/Sqrt((2*Sqrt(2)*r1*x + 3*Power(x,2) + Power(y,2) + 2*Power(z,2))/r2)))/(2.*Power(r2,2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dx:
      J = (x*((-(Sqrt(2)*r1*y) + 2*Power(y,2) + Power(z,2))/
          Sqrt((Power(x,2) - 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2) - 
          (Sqrt(2)*r1*y + 2*Power(y,2) + Power(z,2))/
          Sqrt((Power(x,2) + 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2)))/
          (2.*Power(r2,2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dy:
      J = ((Power(x,2)*(Sqrt(2)*r1 - 2*y) + (Sqrt(2)*r1 - y)*Power(z,2))/
          Sqrt((Power(x,2) - 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2) + 
          (Power(x,2)*(Sqrt(2)*r1 + 2*y) + (Sqrt(2)*r1 + y)*Power(z,2))/
          Sqrt((Power(x,2) + 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2))/
          (2.*Power(r2,2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dz:
      J = (z*((-Power(x,2) + y*(-(Sqrt(2)*r1) + y))/
          Sqrt((Power(x,2) - 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2) + 
          (Power(x,2) - y*(Sqrt(2)*r1 + y))/
          Sqrt((Power(x,2) + 2*Sqrt(2)*r1*y + 3*Power(y,2) + 2*Power(z,2))/r2)))/
          (2.*Power(r2,2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dx:
      J = (2*x/r1 + dR1_dx-dR2_dx)/(2*(R1+R2))-
          (2*r1+R1-R2)*(dR1_dx+dR2_dx)/(2*SQR(R1+R2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      J = (2*y/r1 + dR1_dy-dR2_dy)/(2*(R1+R2))-
          (2*r1+R1-R2)*(dR1_dy+dR2_dy)/(2*SQR(R1+R2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      J = (2*z/r1 + dR1_dz-dR2_dz)/(2*(R1+R2))-
          (2*r1+R1-R2)*(dR1_dz+dR2_dz)/(2*SQR(R1+R2));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
      
    break;
    default:
      abortEr("No such a enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for StereographicSphere Left patch.
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1
*/
double JT_StereographicSphere_Left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
  /* dx/d? not defined */
  else if (q2_e == _x_ || q2_e == _y_|| q2_e == _z_)
    abortEr(INCOMPLETE_FUNC);
  /* dx/d? not defined */
  else if (q1_e == _a_ || q1_e == _b_|| q1_e == _c_)
    abortEr(INCOMPLETE_FUNC);
  
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double 
         x = patch->node[p]->x[0],
         y = patch->node[p]->x[1],
         z = patch->node[p]->x[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2),
         R0 = fabs(patch->c[1]),
         R  = sqrt(r2-SQR(R0));
  double R1 = 0;
  double R2 = 0;
  double S;
  double dS_dx = 0,
         dS_dy = 0,
         dS_dz = 0; 
  double du_dx = 0,
         du_dy = 0,
         du_dz = 0,
         dw_dx = 0,
         dw_dy = 0,
         dw_dz = 0;
  
  /* if dZ_d? we need dR_d? thus: */
  if (q2_e == _b_)
  {
    R1 = patch->CoordSysInfo->R1;
    R2 = patch->CoordSysInfo->R2;
    S = 2./(R2-R1);
  }
  else if (q1_e == _x_)
  {
    S = (R0 - r1)/(R*(R0 - r1 + y));
    dS_dx = -((x*(Power(r2,1.5) + Power(R0,2)*(r1 - y) + R0*(-2*r2 + r1*y)))/
     (r1*Power(-Power(R0,2) + r2,1.5)*Power(R0 - r1 + y,2)));
  }
  else if (q1_e == _y_)
  {
    S = (R0 - r1)/(R*(R0 - r1 + y));
    dS_dy = (Power(R0,3) + r2*(r1 - y) - (Power(R0,2)*(Power(x,2) + r1*y + Power(z,2)))/r1 - 
     R0*(Power(x,2) - 2*r1*y + 2*Power(y,2) + Power(z,2)))/
   (Power(-Power(R0,2) + r2,1.5)*Power(R0 - r1 + y,2));

  }
  else if (q1_e == _z_)
  {
    S = (R0 - r1)/(R*(R0 - r1 + y));
    dS_dz = -(((Power(r2,1.5) + Power(R0,2)*(r1 - y) + R0*(-2*r2 + r1*y))*z)/
     (r1*Power(-Power(R0,2) + r2,1.5)*Power(R0 - r1 + y,2)));
  }
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx +S;
      dw_dx = z*dS_dx;
      J = dX_du_SS(x*S,z*S)*du_dx+dX_dw_SS(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dX_du_SS(x*S,z*S)*du_dy+dX_dw_SS(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dX_du_SS(x*S,z*S)*du_dz+dX_dw_SS(x*S,z*S)*dw_dz;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dx:
      J = S*x/r1;
    break;
    case db_dy:
      J = S*y/r1;
    break;
    case db_dz:
      J = S*z/r1;
    break;
    case dc_dx:
      du_dx = x*dS_dx +S;
      dw_dx = z*dS_dx;
      J = dZ_du_SS(x*S,z*S)*du_dx+dZ_dw_SS(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dZ_du_SS(x*S,z*S)*du_dy+dZ_dw_SS(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dZ_du_SS(x*S,z*S)*du_dz+dZ_dw_SS(x*S,z*S)*dw_dz;
            
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
      
    break;
    default:
      abortEr("No such a enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for StereographicSphere Right patch.
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1
*/
double JT_StereographicSphere_Right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
  /* dx/d? not defined */
  else if (q2_e == _x_ || q2_e == _y_|| q2_e == _z_)
    abortEr(INCOMPLETE_FUNC);
  /* dx/d? not defined */
  else if (q1_e == _a_ || q1_e == _b_|| q1_e == _c_)
    abortEr(INCOMPLETE_FUNC);
  
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double 
         x = patch->node[p]->x[0],
         y = patch->node[p]->x[1],
         z = patch->node[p]->x[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2),
         R0 = fabs(patch->c[1]),
         R  = sqrt(r2-SQR(R0));
  double R1 = 0;
  double R2 = 0;
  double S;
  double dS_dx = 0,
         dS_dy = 0,
         dS_dz = 0; 
  double du_dx = 0,
         du_dy = 0,
         du_dz = 0,
         dw_dx = 0,
         dw_dy = 0,
         dw_dz = 0;
  
  /* if dZ_d? we need dR_d? thus: */
  if (q2_e == _b_)
  {
    R1 = patch->CoordSysInfo->R1;
    R2 = patch->CoordSysInfo->R2;
    S = 2./(R2-R1);
  }
  else if (q1_e == _x_)
  {
    S = (r1-R0)/(R*(r1 + y-R0));
    dS_dx = -((x*(Power(r2,1.5) + Power(R0,2)*(r1 + y) - R0*(2*r2 + r1*y)))/
     (r1*Power(-Power(R0,2) + r2,1.5)*Power(-R0 + r1 + y,2)));
  }
  else if (q1_e == _y_)
  {
    S = (r1-R0)/(R*(r1 + y-R0));
    dS_dy = (-Power(R0,3) - r2*(r1 + y) + Power(R0,2)*(r1 - y - Power(y,2)/r1) + 
     R0*(Power(x,2) + 2*r1*y + 2*Power(y,2) + Power(z,2)))/
   (Power(-Power(R0,2) + r2,1.5)*Power(-R0 + r1 + y,2));

  }
  else if (q1_e == _z_)
  {
    S = (r1-R0)/(R*(r1 + y-R0));
    dS_dz =-(((Power(r2,1.5) + Power(R0,2)*(r1 + y) - R0*(2*r2 + r1*y))*z)/
     (r1*Power(-Power(R0,2) + r2,1.5)*Power(-R0 + r1 + y,2)));
  }
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx +S;
      dw_dx = z*dS_dx;
      J = dX_du_SS(x*S,z*S)*du_dx+dX_dw_SS(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dX_du_SS(x*S,z*S)*du_dy+dX_dw_SS(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dX_du_SS(x*S,z*S)*du_dz+dX_dw_SS(x*S,z*S)*dw_dz;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dx:
      J = S*x/r1;
    break;
    case db_dy:
      J = S*y/r1;
    break;
    case db_dz:
      J = S*z/r1;
    break;
    case dc_dx:
      du_dx = x*dS_dx +S;
      dw_dx = z*dS_dx;
      J = dZ_du_SS(x*S,z*S)*du_dx+dZ_dw_SS(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dZ_du_SS(x*S,z*S)*du_dy+dZ_dw_SS(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dZ_du_SS(x*S,z*S)*du_dz+dZ_dw_SS(x*S,z*S)*dw_dz;
            
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
      
    break;
    default:
      abortEr("No such a enum!\n");
  }
  
  return J;
}

/* Calculating dN(1/2/3)/d(x/y/z) at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN(1/2/3)/d(x/y/z) */
static double dNi_dxj_StereographicSphereRight(Patch_T *const patch, const Dd_T Ni, const Dd_T xj,const double *const X)
{
  double dNi_dxj = 0;
  double J = 0;/* dX?/dx? */
  double a[3];/* Cartesian value of X */
  int rx_of_X = x_of_X(a,X,patch);/* return value of function */
  const double 
         x = a[0],
         y = a[1],
         z = a[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2),
         R0 = fabs(patch->c[1]),
         R  = sqrt(r2-SQR(R0));
  double R1 = 0;
  double R2 = 0;
  double S;
  double dS_dx = 0,
         dS_dy = 0,
         dS_dz = 0; 
  double du_dx = 0,
         du_dy = 0,
         du_dz = 0,
         dw_dx = 0,
         dw_dy = 0,
         dw_dz = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  Dd_T q2_e = UNDEFINED;
  
  if (rx_of_X == 0)
    abortEr("x-coordinate could not been found.");

  if (Ni == _N0_)
    q2_e = _a_;
  else if (Ni == _N1_)
    q2_e = _b_;
  else if (Ni == _N2_)
    q2_e = _c_;
  else
    abortEr(NO_OPTION);
  
  /* if dZ_d? we need dR_d? thus: */
  if (q2_e == _b_)
  {
    R1 = patch->CoordSysInfo->R1;
    R2 = patch->CoordSysInfo->R2;
    S = 2./(R2-R1);
  }
  else if (xj == _x_)
  {
    S = (r1-R0)/(R*(r1 + y-R0));
    dS_dx = -((x*(Power(r2,1.5) + Power(R0,2)*(r1 + y) - R0*(2*r2 + r1*y)))/
     (r1*Power(-Power(R0,2) + r2,1.5)*Power(-R0 + r1 + y,2)));
  }
  else if (xj == _y_)
  {
    S = (r1-R0)/(R*(r1 + y-R0));
    dS_dy = (-Power(R0,3) - r2*(r1 + y) + Power(R0,2)*(r1 - y - Power(y,2)/r1) + 
     R0*(Power(x,2) + 2*r1*y + 2*Power(y,2) + Power(z,2)))/
   (Power(-Power(R0,2) + r2,1.5)*Power(-R0 + r1 + y,2));
  }
  else if (xj == _z_)
  {
    S = (r1-R0)/(R*(r1 + y-R0));
    dS_dz =-(((Power(r2,1.5) + Power(R0,2)*(r1 + y) - R0*(2*r2 + r1*y))*z)/
     (r1*Power(-Power(R0,2) + r2,1.5)*Power(-R0 + r1 + y,2)));
  }
  
  dA_da = get_dA_da(q2_e,xj);
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx +S;
      dw_dx = z*dS_dx;
      J = dX_du_SS(x*S,z*S)*du_dx+dX_dw_SS(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dX_du_SS(x*S,z*S)*du_dy+dX_dw_SS(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dX_du_SS(x*S,z*S)*du_dz+dX_dw_SS(x*S,z*S)*dw_dz;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dx:
      J = S*x/r1;
    break;
    case db_dy:
      J = S*y/r1;
    break;
    case db_dz:
      J = S*z/r1;
    break;
    case dc_dx:
      du_dx = x*dS_dx +S;
      dw_dx = z*dS_dx;
      J = dZ_du_SS(x*S,z*S)*du_dx+dZ_dw_SS(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dZ_du_SS(x*S,z*S)*du_dy+dZ_dw_SS(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dZ_du_SS(x*S,z*S)*du_dz+dZ_dw_SS(x*S,z*S)*dw_dz;
            
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
      
    break;
    default:
      abortEr("No such a enum!\n");
  }
  
  dNi_dxj = dN_dX(patch,Ni,q2_e)*J;
  
  return dNi_dxj;
}

/* Calculating dN(1/2/3)/d(x/y/z) at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN(1/2/3)/d(x/y/z) */
static double dNi_dxj_StereographicSphereLeft(Patch_T *const patch, const Dd_T Ni, const Dd_T xj,const double *const X)
{
  double dNi_dxj = 0;
  double J = 0;/* dX?/dx? */
  double a[3];/* Cartesian value of X */
  int rx_of_X = x_of_X(a,X,patch);/* return value of function */
  const double 
         x = a[0],
         y = a[1],
         z = a[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2),
         R0 = fabs(patch->c[1]),
         R  = sqrt(r2-SQR(R0));
  double R1 = 0;
  double R2 = 0;
  double S;
  double dS_dx = 0,
         dS_dy = 0,
         dS_dz = 0; 
  double du_dx = 0,
         du_dy = 0,
         du_dz = 0,
         dw_dx = 0,
         dw_dy = 0,
         dw_dz = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  Dd_T q2_e = UNDEFINED;
  
  if (rx_of_X == 0)
    abortEr("x-coordinate could not been found.");

  if (Ni == _N0_)
    q2_e = _a_;
  else if (Ni == _N1_)
    q2_e = _b_;
  else if (Ni == _N2_)
    q2_e = _c_;
  else
    abortEr(NO_OPTION);
  
  /* if dZ_d? we need dR_d? thus: */
  if (q2_e == _b_)
  {
    R1 = patch->CoordSysInfo->R1;
    R2 = patch->CoordSysInfo->R2;
    S = 2./(R2-R1);
  }
  else if (xj == _x_)
  {
    S = (R0 - r1)/(R*(R0 - r1 + y));
    dS_dx = -((x*(Power(r2,1.5) + Power(R0,2)*(r1 - y) + R0*(-2*r2 + r1*y)))/
     (r1*Power(-Power(R0,2) + r2,1.5)*Power(R0 - r1 + y,2)));
  }
  else if (xj == _y_)
  {
    S = (R0 - r1)/(R*(R0 - r1 + y));
    dS_dy = (Power(R0,3) + r2*(r1 - y) - (Power(R0,2)*(Power(x,2) + r1*y + Power(z,2)))/r1 - 
     R0*(Power(x,2) - 2*r1*y + 2*Power(y,2) + Power(z,2)))/
   (Power(-Power(R0,2) + r2,1.5)*Power(R0 - r1 + y,2));

  }
  else if (xj == _z_)
  {
    S = (R0 - r1)/(R*(R0 - r1 + y));
    dS_dz = -(((Power(r2,1.5) + Power(R0,2)*(r1 - y) + R0*(-2*r2 + r1*y))*z)/
     (r1*Power(-Power(R0,2) + r2,1.5)*Power(R0 - r1 + y,2)));
  }
  
  dA_da = get_dA_da(q2_e,xj);
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx +S;
      dw_dx = z*dS_dx;
      J = dX_du_SS(x*S,z*S)*du_dx+dX_dw_SS(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dX_du_SS(x*S,z*S)*du_dy+dX_dw_SS(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dX_du_SS(x*S,z*S)*du_dz+dX_dw_SS(x*S,z*S)*dw_dz;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dx:
      J = S*x/r1;
    break;
    case db_dy:
      J = S*y/r1;
    break;
    case db_dz:
      J = S*z/r1;
    break;
    case dc_dx:
      du_dx = x*dS_dx +S;
      dw_dx = z*dS_dx;
      J = dZ_du_SS(x*S,z*S)*du_dx+dZ_dw_SS(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dZ_du_SS(x*S,z*S)*du_dy+dZ_dw_SS(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dZ_du_SS(x*S,z*S)*du_dz+dZ_dw_SS(x*S,z*S)*dw_dz;
            
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
      
    break;
    default:
      abortEr("No such a enum!\n");
  }
  
  dNi_dxj = dN_dX(patch,Ni,q2_e)*J;
  
  return dNi_dxj;
}

/* return value-> dN0_dx for ProjectiveHemisphere */
double dN0_dx_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_ProjectiveHemisphere(patch,_N0_,_x_,X);
}

/* return value-> dN0_dy for ProjectiveHemisphere */
double dN0_dy_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_ProjectiveHemisphere(patch,_N0_,_y_,X);
}
/* return value-> dN0_dz for ProjectiveHemisphere */
double dN0_dz_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_ProjectiveHemisphere(patch,_N0_,_z_,X);
}
/* return value-> dN1_dx for ProjectiveHemisphere */
double dN1_dx_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_ProjectiveHemisphere(patch,_N1_,_x_,X);
}
/* return value-> dN1_dy for ProjectiveHemisphere */
double dN1_dy_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_ProjectiveHemisphere(patch,_N1_,_y_,X);
}
/* return value-> dN1_dz for ProjectiveHemisphere */
double dN1_dz_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_ProjectiveHemisphere(patch,_N1_,_z_,X);
}
/* return value-> dN2_dx for ProjectiveHemisphere */
double dN2_dx_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_ProjectiveHemisphere(patch,_N2_,_x_,X);
}
/* return value-> dN2_dy for ProjectiveHemisphere */
double dN2_dy_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_ProjectiveHemisphere(patch,_N2_,_y_,X);
}
/* return value-> dN2_dz for ProjectiveHemisphere */
double dN2_dz_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_ProjectiveHemisphere(patch,_N2_,_z_,X);
}

/* return value-> dN0_dx for StereographicSphereLeft */
double dN0_dx_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereLeft(patch,_N0_,_x_,X);
}

/* return value-> dN0_dy for StereographicSphereLeft */
double dN0_dy_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereLeft(patch,_N0_,_y_,X);
}
/* return value-> dN0_dz for StereographicSphereLeft */
double dN0_dz_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereLeft(patch,_N0_,_z_,X);
}
/* return value-> dN1_dx for StereographicSphereLeft */
double dN1_dx_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereLeft(patch,_N1_,_x_,X);
}
/* return value-> dN1_dy for StereographicSphereLeft */
double dN1_dy_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereLeft(patch,_N1_,_y_,X);
}
/* return value-> dN1_dz for StereographicSphereLeft */
double dN1_dz_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereLeft(patch,_N1_,_z_,X);
}
/* return value-> dN2_dx for StereographicSphereLeft */
double dN2_dx_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereLeft(patch,_N2_,_x_,X);
}
/* return value-> dN2_dy for StereographicSphereLeft */
double dN2_dy_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereLeft(patch,_N2_,_y_,X);
}
/* return value-> dN2_dz for StereographicSphereLeft */
double dN2_dz_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereLeft(patch,_N2_,_z_,X);
}

/* return value-> dN0_dx for StereographicSphereRight */
double dN0_dx_StereographicSphereRight_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereRight(patch,_N0_,_x_,X);
}

/* return value-> dN0_dy for StereographicSphereRight */
double dN0_dy_StereographicSphereRight_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereRight(patch,_N0_,_y_,X);
}
/* return value-> dN0_dz for StereographicSphereRight */
double dN0_dz_StereographicSphereRight_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereRight(patch,_N0_,_z_,X);
}
/* return value-> dN1_dx for StereographicSphereRight */
double dN1_dx_StereographicSphereRight_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereRight(patch,_N1_,_x_,X);
}
/* return value-> dN1_dy for StereographicSphereRight */
double dN1_dy_StereographicSphereRight_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereRight(patch,_N1_,_y_,X);
}
/* return value-> dN1_dz for StereographicSphereRight */
double dN1_dz_StereographicSphereRight_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereRight(patch,_N1_,_z_,X);
}
/* return value-> dN2_dx for StereographicSphereRight */
double dN2_dx_StereographicSphereRight_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereRight(patch,_N2_,_x_,X);
}
/* return value-> dN2_dy for StereographicSphereRight */
double dN2_dy_StereographicSphereRight_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereRight(patch,_N2_,_y_,X);
}
/* return value-> dN2_dz for StereographicSphereRight */
double dN2_dz_StereographicSphereRight_patch(Patch_T *const patch,const double *const X)
{
  return dNi_dxj_StereographicSphereRight(patch,_N2_,_z_,X);
}

/* getting q2 and q1 coordinate, it turns them to enum_dA_da.
// ->return value: enum_dA_da
*/
enum enum_dA_da get_dA_da(const Dd_T q2_e, const Dd_T q1_e)
{
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  
  if (q2_e == _a_)
  {
   if (q1_e == _x_)
   {
     dA_da = da_dx;
   }
   else if (q1_e == _y_)
   {
     dA_da = da_dy;
   }
   else if (q1_e == _z_)
   {
     dA_da = da_dz;
   }
   else
     abortEr("Invalid entry.");
  }/* end of if (q2_e == _a_) */
  
  else if (q2_e == _b_)
  {
   if (q1_e == _x_)
   {
     dA_da = db_dx;
   }
   else if (q1_e == _y_)
   {
     dA_da = db_dy;
   }
   else if (q1_e == _z_)
   {
     dA_da = db_dz;
   }
   else
     abortEr("Invalid entry.");
  }/* end of else if (q2_e == _b_) */
  
  else if (q2_e == _c_)
  {
   if (q1_e == _x_)
   {
     dA_da = dc_dx;
   }
   else if (q1_e == _y_)
   {
     dA_da = dc_dy;
   }
   else if (q1_e == _z_)
   {
     dA_da = dc_dz;
   }
   else
     abortEr("Invalid entry.");
  }/* end of else if (q2_e == _c_) */
  
  else
    abortEr("Invalid entry.");
  
  return dA_da;  
}

/* preparing R1 and R2 derivatives of projective hemisphere coords.
// NOTE: one must remove dR_? after each updating of R. */
static void R1_R2_derivative(Patch_T *const patch)
{
  Field_T *dR1_dx = 0,
          *dR1_dy = 0,
          *dR2_dx = 0,
          *dR2_dy = 0;
  Field_T *const R1 = patch->pool[Ind("R1_ProjectiveHemisphere")];
  Field_T *const R2 = patch->pool[Ind("R2_ProjectiveHemisphere")];
          
  /* if there is no dR1_dx field */
  if (LookUpField("dR1_dx",patch) == INT_MAX)
  {
    dR1_dx = add_field("dR1_dx",0,patch,NO);
    dR1_dx->v = Partial_Derivative(R1,"x");
  }
  /* if there is no dR1_dy field */
  if (LookUpField("dR1_dy",patch) == INT_MAX)
  {
    dR1_dy = add_field("dR1_dy",0,patch,NO);
    dR1_dy->v = Partial_Derivative(R1,"y");
  }
  /* if there is no dR2_dx field */
  if (LookUpField("dR2_dx",patch) == INT_MAX)
  {
    dR2_dx = add_field("dR2_dx",0,patch,NO);
    dR2_dx->v = Partial_Derivative(R2,"x");
  }
  /* if there is no dR2_dy field */
  if (LookUpField("dR2_dy",patch) == INT_MAX)
  {
    dR2_dy = add_field("dR2_dy",0,patch,NO);
    dR2_dy->v = Partial_Derivative(R2,"y");
  }

}

/* Jacobian transformation for dN/dX?.
// ->return value: dN/dX?
*/
static double dN_dX(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e)
{
  double jN_X = 0;

  if (patch->collocation[q1_e%3] == Chebyshev_Extrema)
  {
    if (q2_e%3 == q1_e%3)
      jN_X = 2.0/(-patch->max[q1_e%3]+patch->min[q1_e%3]); 
    else
      jN_X = 0;
  }
  else if (patch->collocation[q1_e%3] == Chebyshev_Nodes)
  {
    if (q2_e%3 == q1_e%3)
      jN_X = 2.0/(-patch->max[q1_e%3]+patch->min[q1_e%3]); 
    else
      jN_X = 0;
  }
  else
  {
    abortEr(INCOMPLETE_FUNC);
  }

  return jN_X;
}

/* interpolation of 2d field for Projective Hemisphere coordinate.
// Note in Projective Hemisphere 2d fields are f(X,Y,0).
// it mostly is used for radius or derivative of radius.
// convection:
// R = interesting field
// patch = patch that has R
// X = the curvilinear coord in which we want R(X).
// ->return value: R(X)
*/
double interpolation_2d_PH(Field_T *const R, const Patch_T *const patch,const double *const X)
{
  double interp;
  Interpolation_T *interp_s = init_interpolation();
  
  interp_s->field = R;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->XY_dir_flag = 1;
  interp_s->K = 0;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  UNUSED(patch);
  
  return interp;
}

/* dX/du for Stereographic Sphere coord.
// ->return value: dX/du. */
static double dX_du_SS(const double u, const double w)
{
  return ((2*(Sqrt(2) - u))/Sqrt(2 - 2*Sqrt(2)*u + Power(u,2) - Power(w,2)) + 
     (2*(Sqrt(2) + u))/Sqrt(2 + 2*Sqrt(2)*u + Power(u,2) - Power(w,2))
     )/4.;
}

/* dX/dw for Stereographic Sphere coord.
// ->return value: dX/dw. */
static double dX_dw_SS(const double u, const double w)
{
  return (w*(1/Sqrt(2 - 2*Sqrt(2)*u + Power(u,2) - Power(w,2)) - 
       1/Sqrt(2 + 2*Sqrt(2)*u + Power(u,2) - Power(w,2))))/2.;
}

/* dZ/du for Stereographic Sphere coord.
// ->return value: dZ/du. */
static double dZ_du_SS(const double u, const double w)
{
  return (u*(1/Sqrt(2 - Power(u,2) - 2*Sqrt(2)*w + Power(w,2)) - 
       1/Sqrt(2 - Power(u,2) + 2*Sqrt(2)*w + Power(w,2))))/2.;
}

/* dZ/dw for Stereographic Sphere coord.
// ->return value: dZ/dw. */
static double dZ_dw_SS(const double u, const double w)
{
  return ((2*(Sqrt(2) - w))/Sqrt(2 - Power(u,2) - 2*Sqrt(2)*w + Power(w,2)) + 
     (2*(Sqrt(2) + w))/Sqrt(2 - Power(u,2) + 2*Sqrt(2)*w + Power(w,2))
     )/4.;
}
