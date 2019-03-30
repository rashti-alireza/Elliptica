/*
// Alireza Rashti
// June 2018
*/

#include "coordinate_systems.h"

/* making coordinates of nodes */
int make_nodes(Grid_T *const grid)
{
  int i;
  FOR_ALL(i,grid->patch)
  {
    Patch_T *patch = grid->patch[i];
    
    /* if coord is Cartesian */
    if (patch->coordsys == Cartesian)
    {
      make_nodes_Cartesian_coord(patch);
    }
    
    /* if coord is Spherical */
    else if (patch->coordsys == Spherical)
    {
      make_nodes_Spherical_coord(patch);
    }
    
    /* if coord is Protective Hemisphere Up */
    else if (patch->coordsys == ProjectiveHemisphereUp)
    {
      make_nodes_ProjectiveHemisphereUp_coord(patch);
    }
    
    /* if coord is Protective Hemisphere Down */
    else if (patch->coordsys == ProjectiveHemisphereDown)
    {
      make_nodes_ProjectiveHemisphereDown_coord(patch);
    }
    
    /* if coord is Stereographic Sphere Left */
    else if (patch->coordsys == StereographicSphereLeft)
    {
      make_nodes_StereographicSphereLeft_coord(patch);
    }

    /* if coord is Stereographic Sphere Right */
    else if (patch->coordsys == StereographicSphereRight)
    {
      make_nodes_StereographicSphereRight_coord(patch);
    }

    else
      abortEr("No action for such coordinate.\n");
  }
  
  return EXIT_SUCCESS;
}

/* making Jacobian Transformation of coords.
// ->return value: EXIT_SUCCESS.
*/
int make_JacobianT(Grid_T *const grid)
{
  int i;
  FOR_ALL(i,grid->patch)
  {
    Patch_T *patch = grid->patch[i];
    
    /* if coord is Cartesian */
    if (patch->coordsys == Cartesian)
    {
      patch->JacobianT = calloc(1,sizeof(*patch->JacobianT));
      pointerEr(patch->JacobianT);
      make_JacobianT_Cartesian_coord(patch);
    }
    
    else
      abortEr("No job for such coordinate.\n");
  }
  
  return EXIT_SUCCESS;
}

/* making Jacobian transformation for Cartesian coord. */
static void make_JacobianT_Cartesian_coord(Patch_T *const patch)
{
  patch->JacobianT->j      = JT_Cartesian_patch;
  patch->JacobianT->dN0_dx = dN0_dx_Cartesian_patch;
  patch->JacobianT->dN0_dy = dN0_dy_Cartesian_patch;
  patch->JacobianT->dN0_dz = dN0_dz_Cartesian_patch;
  patch->JacobianT->dN1_dx = dN1_dx_Cartesian_patch;
  patch->JacobianT->dN1_dy = dN1_dy_Cartesian_patch;
  patch->JacobianT->dN1_dz = dN1_dz_Cartesian_patch;
  patch->JacobianT->dN2_dx = dN2_dx_Cartesian_patch;
  patch->JacobianT->dN2_dy = dN2_dy_Cartesian_patch;
  patch->JacobianT->dN2_dz = dN2_dz_Cartesian_patch;
}

/* making value of coords. it is a general function for Cartesian type */
static void make_nodes_Cartesian_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const unsigned U = patch->nn;
  const unsigned *const n = patch->n;
  unsigned i,j,k,l;
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  initialize_collocation_struct(patch,&coll_s[2],2);
  
  for (l = 0; l < U; l++)
  {
    double *x = patch->node[l]->x;
    
    IJK(l,n,&i,&j,&k);
    x[0] = point(i,&coll_s[0]);
    x[1] = point(j,&coll_s[1]);
    x[2] = point(k,&coll_s[2]);
    
    /* since X and x are the same we have: */
    patch->node[l]->X = x;
  }
}

/* making value of coords. it is a general function for Spherical type */
static void make_nodes_Spherical_coord(Patch_T *const patch)
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
    X[0] = point(i,&coll_s[0]);/* r */
    X[1] = point(j,&coll_s[1]);/* theta */
    X[2] = point(k,&coll_s[2]);/* phi */
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

/* making value of coords. it is a general function for Protective Hemisphere Up type */
static void make_nodes_ProjectiveHemisphereUp_coord(Patch_T *const patch)
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
    X[0] = point(i,&coll_s[0]);
    X[1] = point(j,&coll_s[1]);
    X[2] = point(k,&coll_s[2]);
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
static void make_nodes_ProjectiveHemisphereDown_coord(Patch_T *const patch)
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
    X[0] = point(i,&coll_s[0]);
    X[1] = point(j,&coll_s[1]);
    X[2] = point(k,&coll_s[2]);
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
static void make_nodes_StereographicSphereLeft_coord(Patch_T *const patch)
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
    X[0] = point(i,&coll_s[0]);
    X[1] = point(j,&coll_s[1]);
    X[2] = point(k,&coll_s[2]);
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
static void make_nodes_StereographicSphereRight_coord(Patch_T *const patch)
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
    X[0] = point(i,&coll_s[0]);
    X[1] = point(j,&coll_s[1]);
    X[2] = point(k,&coll_s[2]);
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

/* initializing collocation struct 
// for making coords based on type of collocation.
//
// transformations:
// ========================
//
// note: min and max refer to the minimum and maximum of a line respectively,
// and n refers to number of points on that line.
//
// o. EquiSpaced:
//      dividing a line in equal sizes
//	-> grid space = (max-min)/(n-1)
//
// o. Chebyshev_Extrema:
//	mapping a line [min,max] to [1,-1] and then using Chebyshev extrema
//	in [0,Pi] to find the nodes on the line. the following map is used:
// 	x = a*N+b = a*cos(t)+b, t = i*Pi/(n-1) where i = 0,...,n-1; i.e.
//	x = 0.5*(-max+min)*N + 0.5*(max+min) where x in [min,max] and N in [-1,1].
// 	N = cos(t). Note that the order of x and N are in reverse. the reason is that
//	it x MUST increase by i. it's crucial for interface realization.
//
// o. Chebyshev_Node:
//	mapping a line [min,max] to [1,-1] and then using Chebyshev nodes
//	in [0,Pi] to find the nodes on the line. the following map is used:
// 	x = a*N+b = a*cos(t)+b, t = (2i+1)*Pi/(2n) where i = 0,...,n-1; i.e.
//	x = 0.5*(-max+min)*N + 0.5*(max+min) where x in (min,max) and N in (-1,1).
// 	N = cos(t). Note that the order of x and N are in reverse. the reason is that
//	it x MUST increase by i. it's crucial for interface realization.
//	
*/
static void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const coll_s,const unsigned dir)
{
  /* some assertions */
  assert(dir < 3);
  assert(LSS(patch->min[dir],patch->max[dir]));
  
  coll_s->c = patch->collocation[dir];
  coll_s->n = patch->n[dir];
  coll_s->min = patch->min[dir];
  coll_s->max = patch->max[dir];
  assert(coll_s->n-1 != 0);

  if (coll_s->c == EquiSpaced)
  {
    coll_s->stp = (coll_s->max-coll_s->min)/(coll_s->n-1);
  }
  else if (coll_s->c == Chebyshev_Extrema)
  {
    coll_s->a = 0.5*(-coll_s->max+coll_s->min);
    coll_s->b = 0.5*(coll_s->max+coll_s->min);
  }
  else if (coll_s->c == Chebyshev_Nodes)
  {
    coll_s->a = 0.5*(-coll_s->max+coll_s->min);
    coll_s->b = 0.5*(coll_s->max+coll_s->min);
  }
  else
    abortEr("There is no such COLLOCATION.\n");
}

/* point value based on collocation 
// ->return value: collocation point on i-th location.
*/
static double point(const unsigned i, const struct Collocation_s *const coll_s)
{
  double x = DBL_MAX;
  
  /* x = min + i*grid_size */ 
  if (coll_s->c == EquiSpaced)
  {
    x = coll_s->min + coll_s->stp*i;
  }
  /* x =  a*N+b => x = a*cos(t)+b */
  else if (coll_s->c == Chebyshev_Extrema)
  {
    double t = i*M_PI/(coll_s->n-1);
    
    x = coll_s->a*cos(t)+coll_s->b;
  }
  /* x = a*N+b => x = a*cos(t)+b */
  else if (coll_s->c == Chebyshev_Nodes)
  {
    double t = (2*i+1)*M_PI/coll_s->n/2.;
    
    x = coll_s->a*cos(t)+coll_s->b;
  }
  else
    abortEr("There is no such collocation.\n");
  
  return x;
}

/* dq2/dq1 Jacobian transformaion for coordinates,
// where q? could be any N0,N1,N2,x,y,z,a,b,c.
// list of q?_e:
// =============
//
// _N0_ for normalized 0-coord [-1,1] 
// _N1_ for normalized 1-coord [-1,1] 
// _N2_ for normalized 2-coord [-1,1] 
// _x_  for Carteisian 0-coord 
// _y_  for Carteisian 1-coord 
// _z_  for Carteisian 2-coord 
// _a_  for Curvilinear 0-coord
// _b_  for Curvilinear 1-coord
// _c_  for Curvilinear 2-coord
//
// o. q2_e is for q2 direction
// o. q1_e is for q1 direction
// o. p is the interested point.
//
// ->return value = dq2/dq1.
*/
double dq2_dq1(const Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  double j = 0;
  
  /* quick answer */
  if (q2_e == q1_e)
    return 1;
    
  else if (q2_e == _N0_ || q2_e == _N1_ || q2_e == _N2_ )
  {
    j = dN_dq(patch,q2_e,q1_e,p);
  }
  else if (q1_e == _N0_ || q1_e == _N1_ || q1_e == _N2_ )
  {
    abortEr(INCOMPLETE_FUNC);
  }
  /* this part means q2_e and q1_e are from x,y,z or a,b,c */
  else
    return patch->JacobianT->j(patch,q2_e,q1_e,p);
  
  return j;
}

/* Jacobian transformation for dN/dX?.
// ->return value: dN/dX?
*/
static double dN_dX(const Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  double jN_X = 0;
  
  if (patch->collocation[q1_e%3] == Chebyshev_Extrema)
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
  
  UNUSED(p);
  
  return jN_X;
}

/* Jacobian transformation for dN/dq?.
// ->return value: dN/dq?
*/
static double dN_dq(const Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  double jN_X = 0;
  
  if (q2_e%3 == q1_e%3)/* e.g. _N2_%3 = 1 and _c_%3 = 1 */
  {
    /* means dN?/dx? = dN?/da*da/dx? + dN?/db*db/dx? + dN?/dc*dc/dx? */
    if (q1_e == _x_ || q1_e == _y_ || q1_e == _z_ )
    {
      jN_X = dN_dX(patch,q2_e,_a_,p)*dq2_dq1(patch,_a_,q1_e,p)+
             dN_dX(patch,q2_e,_b_,p)*dq2_dq1(patch,_b_,q1_e,p)+
             dN_dX(patch,q2_e,_c_,p)*dq2_dq1(patch,_c_,q1_e,p);
              
    }
    else /* means q1_e is between _a_, _b_ or _c_*/
    {
      return dN_dX(patch,q2_e,q1_e,p);
    }
  }
  else
    return 0;
  
  return jN_X;
}


/* Jacobian transformation for Cartesian patch.
// ->return value: dq2/dq1
*/
double JT_Cartesian_patch(const Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  double j;
  
  if (q2_e%3 == q1_e%3)/* e.g. _y_%3 = 1 and _b_%3 = 1 */
    j = 1;
  else
    j = 0;
    
  UNUSED(patch);
  UNUSED(p);
  
  return j;
}

/* Calculating dN0/dx at arbitrary curvilinear point point X.
// used for interpolation.
// x = 0.5*(min-max)*N0 +0.5*(min+max)
// ->return value: dN0/dx */
double dN0_dx_Cartesian_patch(const Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  return 2./(-patch->max[0]+patch->min[0]);
}

/* Calculating dN0/dy at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN0/dy */
double dN0_dy_Cartesian_patch(const Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN0/dz at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN0/dz */
double dN0_dz_Cartesian_patch(const Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN1/dx at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN1/dx */
double dN1_dx_Cartesian_patch(const Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN1/dy at arbitrary curvilinear point point X.
// used for interpolation.
// y = 0.5*(min-max)*N1 +0.5*(min+max)
// ->return value: dN1/dy */
double dN1_dy_Cartesian_patch(const Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  return 2./(-patch->max[1]+patch->min[1]);
}
/* Calculating dN1/dz at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN1/dz */
double dN1_dz_Cartesian_patch(const Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN0/dx at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN0/dx */
double dN2_dx_Cartesian_patch(const Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN0/dx at arbitrary curvilinear point point X.
// used for interpolation.
// ->return value: dN0/dx */
double dN2_dy_Cartesian_patch(const Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  UNUSED(patch);
  return 0;
}
/* Calculating dN0/dx at arbitrary curvilinear point point X.
// used for interpolation.
// z = 0.5*(min-max)*N2 +0.5*(min+max)
// ->return value: dN0/dx */
double dN2_dz_Cartesian_patch(const Patch_T *const patch,const double *const X)
{
  UNUSED(X);
  return 2./(-patch->max[2]+patch->min[2]);
}

/* given patch, general coord of a point and its direction,
// it will change the coord to Chebyshev Extrema format.
// note: X = 0.5*(min-max)*N +0.5*(min+max)
// ->return value: Chebyshe Extrema point.
*/
double General2ChebyshevExtrema(const double X,const unsigned dir,const Patch_T *const patch)
{
  if (patch->basis[dir]       != Chebyshev_Tn_BASIS && 
      patch->collocation[dir] != Chebyshev_Extrema    )
    abortEr("This direction at the specified patch doesn't use first kind Chebyshev,\n"
      "and Chebyshev Extrema collocation which this function needs.\n");
    
  const double a = -patch->max[dir]+patch->min[dir];
  const double b =  patch->max[dir]+patch->min[dir];
  
  return (X*2-b)/a;
}
