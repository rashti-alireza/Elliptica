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
      patch->JacobianT = calloc(1,sizeof(patch->JacobianT));
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
  patch->JacobianT->j = JT_Cartesian_patch;
}

/* making value of coords. it is a general function for Cartesian type */
static void make_nodes_Cartesian_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const unsigned U = countf(patch->node);
  unsigned i,j,k,l,*n;
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  initialize_collocation_struct(patch,&coll_s[2],2);
  n = patch->n;
  
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
// 	x = a*N+b = a*cos(t)+b i.e.
//	x = 0.5*(-max+min)*N + 0.5*(max+min) where x in [min,max] and N in [-1,1].
// 	N = cos(t). Note that the order of x and N are in reverse. the reason is that
//	it x MUST increase by i. it's crucial for interface realization.
//
//	
*/
static void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const coll_s,const unsigned dir)
{
  /* some assertions */
  assert(dir < 3);
  assert(patch->min[dir] < patch->max[dir]);
  
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
