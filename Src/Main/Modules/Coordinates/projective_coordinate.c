/*
// Alireza Rashti
// March 2019
*/

#include "projective_coordinate.h"

/* filling patch struct for BNS_Projective_grid */
void fill_patches_BNS_Projective_grid(Grid_T *const grid)
{
  const unsigned N_outermost_split = (unsigned) GetParameterI_E("Number_of_Outermost_Split");
  unsigned pn,i;
  
  pn = 0;
  populate_left_NS_central_box(grid,pn++);
  populate_left_NS_hemisphere_up(grid,pn++);
  populate_left_NS_hemisphere_down(grid,pn++);
  populate_left_NS_surrounding_up(grid,pn++);
  populate_left_NS_surrounding_down(grid,pn++);
  for (i = 0; i < N_outermost_split; i++)
    populate_left_outermost(grid,pn++,i);
    
  populate_right_NS_central_box(grid,pn++);
  populate_right_NS_hemisphere_up(grid,pn++);
  populate_right_NS_hemisphere_down(grid,pn++);
  populate_right_NS_surrounding_up(grid,pn++);
  populate_right_NS_surrounding_down(grid,pn++);
  for (i = 0; i < N_outermost_split; i++)
    populate_right_outermost(grid,pn++,i);
  
}

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
  const double *const c = patch->c;
  const double R0 = -c[1];
  const double R1 = patch->CoordSysInfo->ProjectiveCoord->R1;
  const double R2 = patch->CoordSysInfo->ProjectiveCoord->R2;
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
    double r,R,A,Co;
    
    IJK(l,n,&i,&j,&k);
    X[0] = point_value(i,&coll_s[0]);
    X[1] = point_value(j,&coll_s[1]);
    X[2] = point_value(k,&coll_s[2]);
    patch->node[l]->X = X;
    
    r = 0.5*X[1]*(R2-R1)+0.5*(R2+R1);
    R = sqrt(SQR(r)-SQR(R0)); assert(!isnan(R));
    u = R*X[0]*sqrt(1-0.5*SQR(X[2])); assert(!isnan(u));
    w = R*X[2]*sqrt(1-0.5*SQR(X[0])); assert(!isnan(w));
    A = SQR(u/(r-R0))+SQR(w/(r-R0))+1;
    x[1] = -r*(2/A-1);
    Co = 2*r/(A*(r-R0));
    x[0] = Co*u;
    x[2] = Co*w;
    
    x[0] += c[0];
    x[1] += c[1];
    x[2] += c[2];
  }
}

/* making value of coords. it is a general function for Stereographic Sphere Right type 
// projected in y = 0 plane for a sphere at center (0,R0,0) */
void make_nodes_StereographicSphereRight_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const unsigned U = patch->nn;
  const double *const c = patch->c;
  const double R0 = c[1];
  const double R1 = patch->CoordSysInfo->ProjectiveCoord->R1;
  const double R2 = patch->CoordSysInfo->ProjectiveCoord->R2;
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
    double r,R,A,Co;
    
    IJK(l,n,&i,&j,&k);
    X[0] = point_value(i,&coll_s[0]);
    X[1] = point_value(j,&coll_s[1]);
    X[2] = point_value(k,&coll_s[2]);
    patch->node[l]->X = X;
    
    r = 0.5*X[1]*(R2-R1)+0.5*(R2+R1);
    R = sqrt(SQR(r)-SQR(R0)); assert(!isnan(R));
    u = R*X[0]*sqrt(1-0.5*SQR(X[2])); assert(!isnan(u));
    w = R*X[2]*sqrt(1-0.5*SQR(X[0])); assert(!isnan(w));
    A = SQR(u/(r-R0))+SQR(w/(r-R0))+1;
    x[1] = r*(2/A-1);
    Co = 2*r/(A*(r-R0));
    x[0] = Co*u;
    x[2] = Co*w;
    
    x[0] += c[0];
    x[1] += c[1];
    x[2] += c[2];
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
  R1_R2_derivative(patch);
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
  const double *const c = patch->c;
  const double 
         x = a[0]-c[0],
         y = a[1]-c[1],
         z = a[2]-c[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2);
  double S,dS_dx,dS_dy,dS_dz,
         u,v,r3,
         du_dx,du_dy,du_dz,
         dv_dx,dv_dy,dv_dz,
         R1,R2,dR2_dx,dR2_dy,dR2_dz,
         dR1_dx,dR1_dy,dR1_dz;
  
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
    Field_T *const R1_field = patch->pool[Ind("R1_ProjectiveHemisphere")];
    Field_T *const R2_field = patch->pool[Ind("R2_ProjectiveHemisphere")];
    Field_T *const dR1_dx_field = patch->CoordSysInfo->ProjectiveCoord->dR1_dx;
    Field_T *const dR1_dy_field = patch->CoordSysInfo->ProjectiveCoord->dR1_dy;
    Field_T *const dR1_dz_field = patch->CoordSysInfo->ProjectiveCoord->dR1_dz;
    Field_T *const dR2_dx_field = patch->CoordSysInfo->ProjectiveCoord->dR2_dx;
    Field_T *const dR2_dy_field = patch->CoordSysInfo->ProjectiveCoord->dR2_dy;
    Field_T *const dR2_dz_field = patch->CoordSysInfo->ProjectiveCoord->dR2_dz;
    
    R1 	   = interpolation_2d_PH(R1_field,patch,X);
    dR1_dx = interpolation_2d_PH(dR1_dx_field,patch,X);
    dR1_dy = interpolation_2d_PH(dR1_dy_field,patch,X);
    dR1_dz = interpolation_2d_PH(dR1_dz_field,patch,X);
    R2 	   = interpolation_2d_PH(R2_field,patch,X);
    dR2_dx = interpolation_2d_PH(dR2_dx_field,patch,X);
    dR2_dy = interpolation_2d_PH(dR2_dy_field,patch,X);
    dR2_dz = interpolation_2d_PH(dR2_dz_field,patch,X);
  }
  else
  {
    r3 = Power3(r1);
    S = 1./r1;
    dS_dx = -x/r3;
    dS_dy = -y/r3;
    dS_dz = -z/r3;
    u = x/r1;
    v = y/r1;
  }
  
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx+S;
      dv_dx = y*dS_dx;
      J = dX_du(u,v)*du_dx+dX_dv(u,v)*dv_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dv_dy = y*dS_dy+S;
      J = dX_du(u,v)*du_dy+dX_dv(u,v)*dv_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dv_dz = y*dS_dz;
      J = dX_du(u,v)*du_dz+dX_dv(u,v)*dv_dz;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dx:
      du_dx = x*dS_dx+S;
      dv_dx = y*dS_dx;
      J = dY_du(u,v)*du_dx+dY_dv(u,v)*dv_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dy:
      du_dy = x*dS_dy;
      dv_dy = y*dS_dy+S;
      J = dY_du(u,v)*du_dy+dY_dv(u,v)*dv_dy;
            
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dz:
      du_dz = x*dS_dz;
      dv_dz = y*dS_dz;
      J = dY_du(u,v)*du_dz+dY_dv(u,v)*dv_dz;
            
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dx:
      J = 2*(x/r1/(R2-R1) - (r1*(dR2_dx-dR1_dx)-R2*dR1_dx+R1*dR2_dx)/SQR(R2-R1));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      J = 2*(y/r1/(R2-R1) - (r1*(dR2_dy-dR1_dy)-R2*dR1_dy+R1*dR2_dy)/SQR(R2-R1));
          
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      J = 2*(z/r1/(R2-R1) - (r1*(dR2_dz-dR1_dz)-R2*dR1_dz+R1*dR2_dz)/SQR(R2-R1));
      
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
  const double *const c = patch->c;
  const double 
         x = patch->node[p]->x[0]-c[0],
         y = patch->node[p]->x[1]-c[1],
         z = patch->node[p]->x[2]-c[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2);
  const unsigned *const n = patch->n;
  double S,dS_dx,dS_dy,dS_dz,
         u,v,r3,
         du_dx,du_dy,du_dz,
         dv_dx,dv_dy,dv_dz,
         R1,R2,dR2_dx,dR2_dy,dR2_dz,
         dR1_dx,dR1_dy,dR1_dz;
  
  /* if dZ_d? we need dR_d? thus: */
  if (q2_e == _c_)
  {
    unsigned i,j,k;
    IJK(p,n,&i,&j,&k);
    
    R1 = patch->CoordSysInfo->ProjectiveCoord->R1_f->v[L(n,i,j,0)];
    R2 = patch->CoordSysInfo->ProjectiveCoord->R2_f->v[L(n,i,j,0)];
    dR2_dx = patch->CoordSysInfo->ProjectiveCoord->dR2_dx->v[L(n,i,j,0)];
    dR2_dy = patch->CoordSysInfo->ProjectiveCoord->dR2_dy->v[L(n,i,j,0)];
    dR2_dz = patch->CoordSysInfo->ProjectiveCoord->dR2_dz->v[L(n,i,j,0)];
    dR1_dx = patch->CoordSysInfo->ProjectiveCoord->dR1_dx->v[L(n,i,j,0)];
    dR1_dy = patch->CoordSysInfo->ProjectiveCoord->dR1_dy->v[L(n,i,j,0)];
    dR1_dz = patch->CoordSysInfo->ProjectiveCoord->dR1_dz->v[L(n,i,j,0)];
  }
  else
  {
    r3 = Power3(r1);
    S = 1./r1;
    dS_dx = -x/r3;
    dS_dy = -y/r3;
    dS_dz = -z/r3;
    u = x/r1;
    v = y/r1;
  }
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx+S;
      dv_dx = y*dS_dx;
      J = dX_du(u,v)*du_dx+dX_dv(u,v)*dv_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dv_dy = y*dS_dy+S;
      J = dX_du(u,v)*du_dy+dX_dv(u,v)*dv_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dv_dz = y*dS_dz;
      J = dX_du(u,v)*du_dz+dX_dv(u,v)*dv_dz;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dx:
      du_dx = x*dS_dx+S;
      dv_dx = y*dS_dx;
      J = dY_du(u,v)*du_dx+dY_dv(u,v)*dv_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dy:
      du_dy = x*dS_dy;
      dv_dy = y*dS_dy+S;
      J = dY_du(u,v)*du_dy+dY_dv(u,v)*dv_dy;
            
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dz:
      du_dz = x*dS_dz;
      dv_dz = y*dS_dz;
      J = dY_du(u,v)*du_dz+dY_dv(u,v)*dv_dz;
            
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dx:
      J = 2*(x/r1/(R2-R1) - (r1*(dR2_dx-dR1_dx)-R2*dR1_dx+R1*dR2_dx)/SQR(R2-R1));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      J = 2*(y/r1/(R2-R1) - (r1*(dR2_dy-dR1_dy)-R2*dR1_dy+R1*dR2_dy)/SQR(R2-R1));
          
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      J = 2*(z/r1/(R2-R1) - (r1*(dR2_dz-dR1_dz)-R2*dR1_dz+R1*dR2_dz)/SQR(R2-R1));
      
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
  const double *const c = patch->c;
  const double 
         x = patch->node[p]->x[0]-c[0],
         y = patch->node[p]->x[1]-c[1],
         z = patch->node[p]->x[2]-c[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2),
         R0 = -c[1],
         R  = sqrt(r2-SQR(R0));
  double R1 = 0;
  double R2 = 0;
  double S;
  double dS_dx = 0,
         dS_dy = 0,
         dS_dz = 0,
         du_dx = 0,
         du_dy = 0,
         du_dz = 0,
         dw_dx = 0,
         dw_dy = 0,
         dw_dz = 0;
  
  /* if dZ_d? we need dR_d? thus: */
  if (q2_e == _b_)
  {
    R1 = patch->CoordSysInfo->ProjectiveCoord->R1;
    R2 = patch->CoordSysInfo->ProjectiveCoord->R2;
    S = 2./(R2-R1);
  }
  else if (q1_e == _x_)
  {
    S = (r1-R0)/(R*(r1 - y));
    dS_dx = -((x*(Power3(R0) + Power(r2,1.5) - 
            Power2(R0)*y + R0*(-2*r2 + r1*y)))/(r1*Power(-Power2(R0) + r2,1.5)*Power2(-r1 + y)));
  }
  else if (q1_e == _y_)
  {
    S = (r1-R0)/(R*(r1 - y));
    dS_dy = ((-R0 + r1)*(-Power2(R0) + r2 + R0*y))/(r1*Power(-Power2(R0) + r2,1.5)*(r1 - y));

  }
  else if (q1_e == _z_)
  {
    S = (r1-R0)/(R*(r1 - y));
    dS_dz = -(((Power3(R0) + Power(r2,1.5) - Power2(R0)*y + R0*(-2*r2 + r1*y))*z)/
            (r1*Power(-Power2(R0) + r2,1.5)*Power2(-r1 + y)));
  }
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx+S;
      dw_dx = z*dS_dx;
      J = dX_du(x*S,z*S)*du_dx+dX_dw(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dX_du(x*S,z*S)*du_dy+dX_dw(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dX_du(x*S,z*S)*du_dz+dX_dw(x*S,z*S)*dw_dz;
      
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
      J = dZ_du(x*S,z*S)*du_dx+dZ_dw(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dZ_du(x*S,z*S)*du_dy+dZ_dw(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dZ_du(x*S,z*S)*du_dz+dZ_dw(x*S,z*S)*dw_dz;
      
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
  const double *const c = patch->c;
  const double 
         x = patch->node[p]->x[0]-c[0],
         y = patch->node[p]->x[1]-c[1],
         z = patch->node[p]->x[2]-c[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2),
         R0 = c[1],
         R  = sqrt(r2-SQR(R0));
  double R1 = 0;
  double R2 = 0;
  double S;
  double dS_dx = 0,
         dS_dy = 0,
         dS_dz = 0,
         du_dx = 0,
         du_dy = 0,
         du_dz = 0,
         dw_dx = 0,
         dw_dy = 0,
         dw_dz = 0;
  
  /* if dZ_d? we need dR_d? thus: */
  if (q2_e == _b_)
  {
    R1 = patch->CoordSysInfo->ProjectiveCoord->R1;
    R2 = patch->CoordSysInfo->ProjectiveCoord->R2;
    S = 2./(R2-R1);
  }
  else if (q1_e == _x_)
  {
    S = (r1-R0)/(R*(r1 + y));
    dS_dx = -((x*(Power3(R0) + Power(r2,1.5) + Power2(R0)*y - R0*(2*r2 + r1*y)))/
     (r1*Power(-Power2(R0) + r2,1.5)*Power2(r1 + y)));
  }
  else if (q1_e == _y_)
  {
    S = (r1-R0)/(R*(r1 + y));
    dS_dy = -(((-R0 + r1)*(-Power2(R0) + r2 - R0*y))/
        (r1*Power(-Power2(R0) + r2,1.5)*(r1 + y)));

  }
  else if (q1_e == _z_)
  {
    S = (r1-R0)/(R*(r1 + y));
    dS_dz =-(((Power3(R0) + Power(r2,1.5) + Power2(R0)*y - R0*(2*r2 + r1*y))*z)/
     (r1*Power(-Power2(R0) + r2,1.5)*Power2(r1 + y)));
  }
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx+S;
      dw_dx = z*dS_dx;
      J = dX_du(x*S,z*S)*du_dx+dX_dw(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dX_du(x*S,z*S)*du_dy+dX_dw(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dX_du(x*S,z*S)*du_dz+dX_dw(x*S,z*S)*dw_dz;
      
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
      J = dZ_du(x*S,z*S)*du_dx+dZ_dw(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dZ_du(x*S,z*S)*du_dy+dZ_dw(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dZ_du(x*S,z*S)*du_dz+dZ_dw(x*S,z*S)*dw_dz;
      
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
  const double *const c = patch->c;
  const double 
         x = a[0]-c[0],
         y = a[1]-c[1],
         z = a[2]-c[2],
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
    R1 = patch->CoordSysInfo->ProjectiveCoord->R1;
    R2 = patch->CoordSysInfo->ProjectiveCoord->R2;
    S = 2./(R2-R1);
  }
  else if (xj == _x_)
  {
    S = (r1-R0)/(R*(r1 + y));
    dS_dx = -((x*(Power3(R0) + Power(r2,1.5) + Power2(R0)*y - R0*(2*r2 + r1*y)))/
     (r1*Power(-Power2(R0) + r2,1.5)*Power2(r1 + y)));
  }
  else if (xj == _y_)
  {
    S = (r1-R0)/(R*(r1 + y));
    dS_dy = -(((-R0 + r1)*(-Power2(R0) + r2 - R0*y))/
        (r1*Power(-Power2(R0) + r2,1.5)*(r1 + y)));

  }
  else if (xj == _z_)
  {
    S = (r1-R0)/(R*(r1 + y));
    dS_dz =-(((Power3(R0) + Power(r2,1.5) + Power2(R0)*y - R0*(2*r2 + r1*y))*z)/
     (r1*Power(-Power2(R0) + r2,1.5)*Power2(r1 + y)));
  }
  
  dA_da = get_dA_da(q2_e,xj);
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx+S;
      dw_dx = z*dS_dx;
      J = dX_du(x*S,z*S)*du_dx+dX_dw(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dX_du(x*S,z*S)*du_dy+dX_dw(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dX_du(x*S,z*S)*du_dz+dX_dw(x*S,z*S)*dw_dz;
      
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
      J = dZ_du(x*S,z*S)*du_dx+dZ_dw(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dZ_du(x*S,z*S)*du_dy+dZ_dw(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dZ_du(x*S,z*S)*du_dz+dZ_dw(x*S,z*S)*dw_dz;
      
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
  const double *const c = patch->c;
  const double 
         x = a[0]-c[0],
         y = a[1]-c[1],
         z = a[2]-c[2],
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
    R1 = patch->CoordSysInfo->ProjectiveCoord->R1;
    R2 = patch->CoordSysInfo->ProjectiveCoord->R2;
    S = 2./(R2-R1);
  }
  else if (xj == _x_)
  {
    S = (r1-R0)/(R*(r1 - y));
    dS_dx = -((x*(Power3(R0) + Power(r2,1.5) - 
            Power2(R0)*y + R0*(-2*r2 + r1*y)))/(r1*Power(-Power2(R0) + r2,1.5)*Power2(-r1 + y)));
  }
  else if (xj == _y_)
  {
    S = (r1-R0)/(R*(r1 - y));
    dS_dy = ((-R0 + r1)*(-Power2(R0) + r2 + R0*y))/(r1*Power(-Power2(R0) + r2,1.5)*(r1 - y));

  }
  else if (xj == _z_)
  {
    S = (r1-R0)/(R*(r1 - y));
    dS_dz = -(((Power3(R0) + Power(r2,1.5) - Power2(R0)*y + R0*(-2*r2 + r1*y))*z)/
            (r1*Power(-Power2(R0) + r2,1.5)*Power2(-r1 + y)));
  }
  
  dA_da = get_dA_da(q2_e,xj);
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx+S;
      dw_dx = z*dS_dx;
      J = dX_du(x*S,z*S)*du_dx+dX_dw(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dX_du(x*S,z*S)*du_dy+dX_dw(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dX_du(x*S,z*S)*du_dz+dX_dw(x*S,z*S)*dw_dz;
      
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
      J = dZ_du(x*S,z*S)*du_dx+dZ_dw(x*S,z*S)*dw_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      du_dy = x*dS_dy;
      dw_dy = z*dS_dy;
      J = dZ_du(x*S,z*S)*du_dy+dZ_dw(x*S,z*S)*dw_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      du_dz = x*dS_dz;
      dw_dz = z*dS_dz+S;
      J = dZ_du(x*S,z*S)*du_dz+dZ_dw(x*S,z*S)*dw_dz;
      
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

/* preparing R1 and R2 derivatives of projective hemisphere coords.
// NOTE: one must remove dR_? after each updating of R. */
static void R1_R2_derivative(Patch_T *const patch)
{
  Field_T *dR1_dX = add_field("dR1_dX",0,patch,NO),
          *dR1_dY = add_field("dR1_dY",0,patch,NO),
          *dR2_dX = add_field("dR2_dX",0,patch,NO),
          *dR2_dY = add_field("dR2_dY",0,patch,NO),
          *dR1_dx = add_field("dR1_dx",0,patch,YES),
          *dR1_dy = add_field("dR1_dy",0,patch,YES),
          *dR1_dz = add_field("dR1_dz",0,patch,YES),
          *dR2_dx = add_field("dR2_dx",0,patch,YES),
          *dR2_dy = add_field("dR2_dy",0,patch,YES),
          *dR2_dz = add_field("dR2_dz",0,patch,YES);
  Field_T *const R1 = patch->pool[Ind("R1_ProjectiveHemisphere")];
  Field_T *const R2 = patch->pool[Ind("R2_ProjectiveHemisphere")];
  const unsigned nn = patch->nn;
  unsigned p;
          
  dR1_dX->v = Partial_Derivative(R1,"a");
  dR1_dY->v = Partial_Derivative(R1,"b");
  dR2_dX->v = Partial_Derivative(R2,"a");
  dR2_dY->v = Partial_Derivative(R2,"b");
    
  for (p = 0; p < nn; ++p)
  {
    dR1_dx->v[p] = dR1_dX->v[p]*dq2_dq1(patch,_a_,_x_,p)+
                   dR1_dY->v[p]*dq2_dq1(patch,_b_,_x_,p);
    dR1_dy->v[p] = dR1_dX->v[p]*dq2_dq1(patch,_a_,_y_,p)+
                   dR1_dY->v[p]*dq2_dq1(patch,_b_,_y_,p);
    dR1_dz->v[p] = dR1_dX->v[p]*dq2_dq1(patch,_a_,_z_,p)+
                   dR1_dY->v[p]*dq2_dq1(patch,_b_,_z_,p);
    dR2_dx->v[p] = dR2_dX->v[p]*dq2_dq1(patch,_a_,_x_,p)+
                   dR2_dY->v[p]*dq2_dq1(patch,_b_,_x_,p);
    dR2_dy->v[p] = dR2_dX->v[p]*dq2_dq1(patch,_a_,_y_,p)+
                   dR2_dY->v[p]*dq2_dq1(patch,_b_,_y_,p);
    dR2_dz->v[p] = dR2_dX->v[p]*dq2_dq1(patch,_a_,_z_,p)+
                   dR2_dY->v[p]*dq2_dq1(patch,_b_,_z_,p);
    
  }
                      
  remove_field(dR1_dX);
  remove_field(dR1_dY);
  remove_field(dR2_dX);
  remove_field(dR2_dY);
  
  patch->CoordSysInfo->ProjectiveCoord->dR1_dx = dR1_dx;
  patch->CoordSysInfo->ProjectiveCoord->dR1_dy = dR1_dy;
  patch->CoordSysInfo->ProjectiveCoord->dR1_dz = dR1_dz;
  patch->CoordSysInfo->ProjectiveCoord->dR2_dx = dR2_dx;
  patch->CoordSysInfo->ProjectiveCoord->dR2_dy = dR2_dy;
  patch->CoordSysInfo->ProjectiveCoord->dR2_dz = dR2_dz; 
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

/* dX/du for elliptical mapping used for projective coord.
// ->return value: dX/du. */
static double dX_du(const double u, const double v)
{
  return ((2*(Sqrt2 - u))/
      Sqrt(2 - 2*Sqrt2*u + 
        Power2(u) - Power2(v)) + 
     (2*(Sqrt2 + u))/
      Sqrt(2 + 2*Sqrt2*u + 
        Power2(u) - Power2(v)))/4.;
}

/* dX/dv for elliptical mapping used for projective coord.
// ->return value: dX/dv. */
static double dX_dv(const double u, const double v)
{
  return (v*(1/Sqrt(2 - 2*Sqrt2*u + 
          Power2(u) - Power2(v)) - 
       1/
        Sqrt(2 + 2*Sqrt2*u + 
          Power2(u) - Power2(v))))/2.;
}

/* dY/du for elliptical mapping used for projective coord.
// ->return value: dY/du. */
static double dY_du(const double u, const double v)
{
  return (u*(1/Sqrt(2 - Power(u,2) - 
          2*Sqrt2*v + Power(v,2)) - 
       1/
        Sqrt(2 - Power2(u) + 
          2*Sqrt2*v + Power2(v))))/2.;
}

/* dY/du for elliptical mapping used for projective coord.
// ->return value: dY/dv. */
static double dY_dv(const double u, const double v)
{
  return ((2*(Sqrt2 - v))/
      Sqrt(2 - Power2(u) - 
        2*Sqrt2*v + Power2(v)) + 
     (2*(Sqrt2 + v))/
      Sqrt(2 - Power2(u) + 
        2*Sqrt2*v + Power2(v)))/4.;
}


/* populating properties of patch for outermost left */
void populate_left_outermost(Grid_T *const grid,const unsigned pn,const unsigned outermost_n)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_outermost%u",grid->gn,outermost_n);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"Outermost%u",outermost_n);
  make_keyword_parameter(&ret,var,"n");
  n = (unsigned)GetParameterI(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (unsigned)GetParameterI(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (unsigned)GetParameterI(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_outermost%u_R1",grid->gn,outermost_n);
  patch->CoordSysInfo->ProjectiveCoord->R1 = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_outermost%u_R2",grid->gn,outermost_n);
  patch->CoordSysInfo->ProjectiveCoord->R2 = GetParameterDoubleF_E(var);
  
  assert(GRT(patch->CoordSysInfo->ProjectiveCoord->R2,patch->CoordSysInfo->ProjectiveCoord->R1));
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = StereographicSphereLeft;
  
  /* collocation */
  patch->collocation[0] = Chebyshev_Nodes;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
  
}

/* populating properties of patch for right outermost */
void populate_right_outermost(Grid_T *const grid,const unsigned pn,const unsigned outermost_n)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_outermost%u",grid->gn,outermost_n);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"Outermost%u",outermost_n);
  make_keyword_parameter(&ret,var,"n");
  n = (unsigned)GetParameterI(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (unsigned)GetParameterI(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (unsigned)GetParameterI(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_outermost%u_R1",grid->gn,outermost_n);
  patch->CoordSysInfo->ProjectiveCoord->R1 = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_outermost%u_R2",grid->gn,outermost_n);
  patch->CoordSysInfo->ProjectiveCoord->R2 = GetParameterDoubleF_E(var);
  
  assert(GRT(patch->CoordSysInfo->ProjectiveCoord->R2,patch->CoordSysInfo->ProjectiveCoord->R1));
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = StereographicSphereRight;
  
  /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Nodes;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
  
}

/* populating properties of patch for NS left hemisphere up */
static void populate_left_NS_hemisphere_up(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double *R1_array,*R2_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_NS_hemisphere_up",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (unsigned)GetParameterI(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (unsigned)GetParameterI(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (unsigned)GetParameterI(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_NS_R1_up",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  sprintf(var,"grid%u_left_NS_R2_up",grid->gn);
  R2_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->ProjectiveCoord->R1_f = R1;
  patch->CoordSysInfo->ProjectiveCoord->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R2->v[ij0] = R2_array[ij0];
      R1->v[ij0] = R1_array[ij0];
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereUp;
  
  /* collocation */
  patch->collocation[0] = Chebyshev_Nodes;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
  
}

/* populating properties of patch for NS left hemisphere down */
static void populate_left_NS_hemisphere_down(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double *R1_array,*R2_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_NS_hemisphere_down",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (unsigned)GetParameterI(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (unsigned)GetParameterI(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (unsigned)GetParameterI(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_NS_R1_down",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  sprintf(var,"grid%u_left_NS_R2_down",grid->gn);
  R2_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->ProjectiveCoord->R1_f = R1;
  patch->CoordSysInfo->ProjectiveCoord->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R2->v[ij0] = R2_array[ij0];
      R1->v[ij0] = R1_array[ij0];
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereDown;
  
/* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Nodes;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for NS right hemisphere up */
static void populate_right_NS_hemisphere_up(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double *R1_array,*R2_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_NS_hemisphere_up",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (unsigned)GetParameterI(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (unsigned)GetParameterI(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (unsigned)GetParameterI(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_NS_R1_up",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  sprintf(var,"grid%u_right_NS_R2_up",grid->gn);
  R2_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->ProjectiveCoord->R1_f = R1;
  patch->CoordSysInfo->ProjectiveCoord->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R2->v[ij0] = R2_array[ij0];
      R1->v[ij0] = R1_array[ij0];
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereUp;
  
/* collocation */
  patch->collocation[0] = Chebyshev_Nodes;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for NS right hemisphere down */
static void populate_right_NS_hemisphere_down(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double *R1_array,*R2_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_NS_hemisphere_down",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (unsigned)GetParameterI(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (unsigned)GetParameterI(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (unsigned)GetParameterI(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_NS_R1_down",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  sprintf(var,"grid%u_right_NS_R2_down",grid->gn);
  R2_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->ProjectiveCoord->R1_f = R1;
  patch->CoordSysInfo->ProjectiveCoord->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R2->v[ij0] = R2_array[ij0];
      R1->v[ij0] = R1_array[ij0];
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereDown;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Nodes;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for left NS's surrounding up */
static void populate_left_NS_surrounding_up(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double R2_const,*R1_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_NS_surrounding_up",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (unsigned)GetParameterI(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (unsigned)GetParameterI(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (unsigned)GetParameterI(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_NS_Surrounding_R2",grid->gn);
  R2_const = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_R2_up",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->ProjectiveCoord->R1_f = R1;
  patch->CoordSysInfo->ProjectiveCoord->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R1->v[ij0] = R1_array[ij0];
      R2->v[ij0] = R2_const;
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereUp;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Nodes;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for left NS's surrounding down */
static void populate_left_NS_surrounding_down(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double R2_const,*R1_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_left_NS_surrounding_down",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"left_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (unsigned)GetParameterI(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (unsigned)GetParameterI(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (unsigned)GetParameterI(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_left_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_left_NS_Surrounding_R2",grid->gn);
  R2_const = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_left_NS_R2_down",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->ProjectiveCoord->R1_f = R1;
  patch->CoordSysInfo->ProjectiveCoord->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R1->v[ij0] = R1_array[ij0];
      R2->v[ij0] = R2_const;
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereDown;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Nodes;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for right NS's surrounding up */
static void populate_right_NS_surrounding_up(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double R2_const,*R1_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_NS_surrounding_up",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (unsigned)GetParameterI(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (unsigned)GetParameterI(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (unsigned)GetParameterI(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_NS_Surrounding_R2",grid->gn);
  R2_const = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_R2_up",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->ProjectiveCoord->R1_f = R1;
  patch->CoordSysInfo->ProjectiveCoord->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R1->v[ij0] = R1_array[ij0];
      R2->v[ij0] = R2_const;
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereUp;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Nodes;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of patch for right_NS's surrounding down */
static void populate_right_NS_surrounding_down(Grid_T *const grid,const unsigned pn)
{
  Patch_T *const patch = grid->patch[pn];
  Field_T *R1 = add_field("R1_ProjectiveHemisphere",0,patch,NO);
  Field_T *R2 = add_field("R2_ProjectiveHemisphere",0,patch,NO);
  double R2_const,*R1_array;
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  unsigned n,i,j;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,"grid%u_right_NS_surrounding_down",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  patch->n[0] = (unsigned)GetParameterI("n_a");
  patch->n[1] = (unsigned)GetParameterI("n_b");
  patch->n[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  sprintf(var,"right_NS");
  make_keyword_parameter(&ret,var,"n");
  n = (unsigned)GetParameterI(ret.s0);
  if (n != INT_MAX)	patch->n[0] = n;
  n = (unsigned)GetParameterI(ret.s1);
  if (n != INT_MAX)	patch->n[1] = n;
  n = (unsigned)GetParameterI(ret.s2);
  if (n != INT_MAX)	patch->n[2] = n;
    
  if(patch->n[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(patch->n[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(patch->n[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,"grid%u_right_NS_center_a",grid->gn);
  patch->c[0] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_b",grid->gn);
  patch->c[1] = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_center_c",grid->gn);
  patch->c[2] = GetParameterDoubleF_E(var);
  
  /* filling Rs */
  sprintf(var,"grid%u_right_NS_Surrounding_R2",grid->gn);
  R2_const = GetParameterDoubleF_E(var);
  sprintf(var,"grid%u_right_NS_R2_down",grid->gn);
  R1_array = GetParameterArrayF_E(var);
  patch->CoordSysInfo->ProjectiveCoord->R1_f = R1;
  patch->CoordSysInfo->ProjectiveCoord->R2_f = R2;
  
  R1->v = alloc_double(patch->nn);
  R2->v = alloc_double(patch->nn);
  for (i = 0; i < patch->n[0]; ++i)
    for (j = 0; j < patch->n[1]; ++j)
    {
      unsigned ij0 = L(patch->n,i,j,0);
      R1->v[ij0] = R1_array[ij0];
      R2->v[ij0] = R2_const;
    }
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = -1;
  
  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* filling flags */
  patch->coordsys = ProjectiveHemisphereDown;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Nodes;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* memory alloc patches for BNS_Projective type */
void alloc_patches_BNS_Projective_grid(Grid_T *const grid)
{
  unsigned Np = 10;/* number of patches without outermost's*/
  unsigned outermost;
  unsigned i;
  
  outermost = (unsigned) GetParameterI("Number_of_Outermost_Split");
  if (outermost != (unsigned)INT_MAX)
    Np += 2*outermost;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  pointerEr(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    pointerEr(grid->patch[i]);
  }
  
}

