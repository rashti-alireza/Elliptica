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
    
    /* if coord is Cubed Spherical */
    else if (patch->coordsys == CubedSpherical)
    {
      make_nodes_CubedSpherical_coord(patch);
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
    /* if coord is Cubed Spherical */
    else if (patch->coordsys == CubedSpherical)
    {
      patch->JacobianT = calloc(1,sizeof(*patch->JacobianT));
      pointerEr(patch->JacobianT);
      make_JacobianT_CubedSpherical_coord(patch);
    }
    /* if coord is ProjectiveHemisphereUp */
    else if (patch->coordsys == ProjectiveHemisphereUp)
    {
      patch->JacobianT = calloc(1,sizeof(*patch->JacobianT));
      pointerEr(patch->JacobianT);
      make_JacobianT_ProjectiveHemisphere_coord(patch);
    }
    
    /* if coord is ProjectiveHemisphereDown */
    else if (patch->coordsys == ProjectiveHemisphereDown)
    {
      patch->JacobianT = calloc(1,sizeof(*patch->JacobianT));
      pointerEr(patch->JacobianT);
      make_JacobianT_ProjectiveHemisphere_coord(patch);
    }
    /* if coord is StereographicSphereLeft */
    else if (patch->coordsys == StereographicSphereLeft)
    {
      patch->JacobianT = calloc(1,sizeof(*patch->JacobianT));
      pointerEr(patch->JacobianT);
      make_JacobianT_StereographicSphereLeft_coord(patch);
    }
    
    /* if coord is StereographicSphereRight */
    else if (patch->coordsys == StereographicSphereRight)
    {
      patch->JacobianT = calloc(1,sizeof(*patch->JacobianT));
      pointerEr(patch->JacobianT);
      make_JacobianT_StereographicSphereRight_coord(patch);
    }
    
    else
      abortEr("No job for such coordinate.\n");
  }
  
  return EXIT_SUCCESS;
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
void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const coll_s,const unsigned dir)
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
double point_value(const unsigned i, const struct Collocation_s *const coll_s)
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
    double t = (i+0.5)*M_PI/coll_s->n;
    
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
double dq2_dq1(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
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
static double dN_dX(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
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
  
  UNUSED(p);
  
  return jN_X;
}

/* Jacobian transformation for dN/dq?.
// ->return value: dN/dq?
*/
static double dN_dq(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  double jN_X = 0;
  
  /* means dN?/dx? = dN?/da*da/dx? + dN?/db*db/dx? + dN?/dc*dc/dx? */
  if (q1_e == _x_ || q1_e == _y_ || q1_e == _z_ )
  {
    /* it means N only depends on one of a,b or c */
    if (patch->basis[q2_e%3] == Chebyshev_Tn_BASIS)
    {
      switch(q2_e)
      {
        case _N0_:
          jN_X = dN_dX(patch,q2_e,_a_,p)*dq2_dq1(patch,_a_,q1_e,p);
        break;
        case _N1_:
          jN_X = dN_dX(patch,q2_e,_b_,p)*dq2_dq1(patch,_b_,q1_e,p);
        break;
        case _N2_:
          jN_X = dN_dX(patch,q2_e,_c_,p)*dq2_dq1(patch,_c_,q1_e,p);
        break;
        default:
          abortEr("It should not reach here!");
      }
    }
    else
      jN_X = dN_dX(patch,q2_e,_a_,p)*dq2_dq1(patch,_a_,q1_e,p)+
             dN_dX(patch,q2_e,_b_,p)*dq2_dq1(patch,_b_,q1_e,p)+
             dN_dX(patch,q2_e,_c_,p)*dq2_dq1(patch,_c_,q1_e,p);
            
  }
  else /* means q1_e is between _a_, _b_ or _c_*/
  {
    return dN_dX(patch,q2_e,q1_e,p);
  }
  
  return jN_X;
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

/* making keyword for searching of parameter value */
void make_keyword_parameter(struct Ret_S *const ret,const char *const box,const char *const needle)
{
  /* for box?_n_? */
  if (strcmp_i(needle,"n"))
  {
    sprintf(ret->s0,"%s_n_a",box);
    sprintf(ret->s1,"%s_n_b",box);
    sprintf(ret->s2,"%s_n_c",box);
  }
  
  /* for box?_center_? */
  else if (strcmp_i(needle,"center"))
  {
    sprintf(ret->s0,"%s_center_a",box);
    sprintf(ret->s1,"%s_center_b",box);
    sprintf(ret->s2,"%s_center_c",box);
  }
  
  /* for box?_size_? */
  else if (strcmp_i(needle,"size"))
  {
    sprintf(ret->s0,"%s_size_a",box);
    sprintf(ret->s1,"%s_size_b",box);
    sprintf(ret->s2,"%s_size_c",box);
  }
  /* for box?_collocation_? */
  else if (strcmp_i(needle,"collocation"))
  {
    sprintf(ret->s0,"%s_collocation_a",box);
    sprintf(ret->s1,"%s_collocation_b",box);
    sprintf(ret->s2,"%s_collocation_c",box);
  }
  /* for box?_basis_? */
  else if (strcmp_i(needle,"basis"))
  {
    sprintf(ret->s0,"%s_basis_a",box);
    sprintf(ret->s1,"%s_basis_b",box);
    sprintf(ret->s2,"%s_basis_c",box);
  }
  else
  {
    abortEr_s("There is no such %s.\n",needle);
  }
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

/* calculating the main characteristics of grid needed for setting up
// the patches, like lengths, radii and etc for tutorial purposes. */
void grid_characteristics_example(Grid_T *const grid)
{
  const char *kind;
  
  /* finding the kind of grid */
  kind = GetParameterS_E("grid_kind");
  grid->kind = dup_s(kind);
  
  if (strcmp_i(grid->kind,"Cartesian_grid"))
    characteristics_Cartesian_grid(grid);
  
  else if (strcmp_i(grid->kind,"BNS_Projective_grid"))
    characteristics_BNS_Projective_grid(grid); 
    
  else if (strcmp_i(grid->kind,"BNS_Spherical_grid"))
    characteristics_BNS_Spherical_grid(grid); 
    
  else if (strcmp_i(grid->kind,"BNS_CubedSpherical_grid"))
    characteristics_BNS_CubedSpherical_grid(grid); 
   
  else
    abortEr_s("There is no such %s grid kind.\n",grid->kind);
  
}

/* calculating the main characteristic of grid for Cartesian grid */
static void characteristics_Cartesian_grid(Grid_T *const grid)
{
  /* this type of grid is so simple; nothing to calculate. */
  UNUSED(grid);
}

/* calculating the main characteristic of grid for BNS_CubedSpherical grid */
static void characteristics_BNS_CubedSpherical_grid(Grid_T *const grid)
{
  const unsigned gn   = grid->gn;
  const double CONST  = 5.0;
  const double C      = GetParameterD_E("BNS_Distance");
  const double R_NS_l = GetParameterD_E("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = GetParameterD_E("right_NS_radius");/* assuming perfect sphere */
  double box_size_l,box_size_r;
  const unsigned N_Outermost_Split = (unsigned)GetParameterI_E("Number_of_Outermost_Split"); 
  double *R_outermost = calloc(N_Outermost_Split,sizeof(*R_outermost));
  unsigned n_box_l, n_box_r,n_r,n_l;
  char var[100] = {'\0'};
  char par[100] = {'\0'};
  char val[100] = {'\0'};
  unsigned i;
  
  assert(GRT(C,0));
  assert(GRT(R_NS_l,0));
  assert(GRT(R_NS_r,0));
  assert(LSS(2*R_NS_l,C));
  assert(LSS(2*R_NS_r,C));
  
  n_l = (unsigned)GetParameterI("n_c");
  i   = (unsigned)GetParameterI("left_NS_n_c");
  if (i != INT_MAX) 	n_l = i;
  if (n_l == INT_MAX)   abortEr("n_l could not be set.");
  assert(n_l > 2);
  
  n_r = (unsigned)GetParameterI("n_c");
  i   = (unsigned)GetParameterI("right_NS_n_c");
  if (i != INT_MAX) 	n_r = i;
  if (n_r == INT_MAX)   abortEr("n_r could not be set.");
  assert(n_r > 2);
  
  /* making NS's surface function */
  NS_surface_BNS_CubedSpherical_grid(grid);
  
  box_size_l = 2*R_NS_l/n_l;
  box_size_r = 2*R_NS_r/n_r;
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R_outermost[i] = GetParameterD_E(var);
    
    if (LSS(R_outermost[i],2*C))
      abortEr("the radius of outermost patches must be greater than twice of BNS distance.");
    
    if (i > 0)
      if (LSSEQL(R_outermost[i],R_outermost[i-1]))
        abortEr("The radius of outermost must be increasing.");
    
  }
  
  n_box_l = (unsigned)(3 + n_l/CONST);
  n_box_r = (unsigned)(3 + n_r/CONST);
  
  /* adding the results to the parameter data base */
  
  /* n_a, n_b, n_c */
  sprintf(par,"grid%u_left_centeral_box_n_abc",gn);
  sprintf(val,"%u",n_box_l);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_centeral_box_n_abc",gn);
  sprintf(val,"%u",n_box_r);
  add_parameter_string(par,val);
  
  /* size a,b,c */
  sprintf(par,"grid%u_left_centeral_box_size_a",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_b",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_c",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_right_centeral_box_size_a",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_right_centeral_box_size_b",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_right_centeral_box_size_c",gn);
  add_parameter_double(par,box_size_r);
  
  /* surrounding box length */
  sprintf(par,"grid%u_surrounding_box_length",gn);
  add_parameter_double(par,C);
  
  /* R1 and R2 outermost */
  sprintf(par,"grid%u_outermost%u_R2",gn,0);
  add_parameter_double(par,R_outermost[0]);
    
  for (i = 1; i < N_Outermost_Split; i++)
  {
    /* R1: */
    sprintf(par,"grid%u_outermost%u_R1",gn,i);
    add_parameter_double(par,R_outermost[i-1]);
    
    /* R2: */
    sprintf(par,"grid%u_outermost%u_R2",gn,i);
    add_parameter_double(par,R_outermost[i]);
    
  }
  
  /* assuming the center of left NS at (0,-C/2,0) */
  sprintf(par,"grid%u_left_NS_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_left_NS_center_b",gn);
  add_parameter_double(par,-C/2);
  
  sprintf(par,"grid%u_left_NS_center_c",gn);
  add_parameter_double(par,0.0);
  
  /* assuming the center of right NS at (0,C/2,0) */
  sprintf(par,"grid%u_right_NS_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_right_NS_center_b",gn);
  add_parameter_double(par,C/2);
  
  sprintf(par,"grid%u_right_NS_center_c",gn);
  add_parameter_double(par,0.0);
  
  free(R_outermost);
}

/* calculating the main characteristic of grid for BNS_Projective grid */
static void characteristics_BNS_Projective_grid(Grid_T *const grid)
{
  const unsigned gn   = grid->gn;
  const double CONST  = 5.0;
  const double C      = GetParameterD_E("BNS_Distance");
  const double R_NS_l = GetParameterD_E("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = GetParameterD_E("right_NS_radius");/* assuming perfect sphere */
  struct R_max r_max[1] = {0};
  double R_max_l,R_max_r;/* maximum distance from the center of each star */
  const unsigned N_Outermost_Split = (unsigned)GetParameterI_E("Number_of_Outermost_Split"); 
  double O,O_l,O_r,
         R_Surr_l,R_Surr_r,
         box_size_l,box_size_r,
         *R_outmost_l = alloc_double(N_Outermost_Split),
         *R_outmost_r = alloc_double(N_Outermost_Split),
         *R0 = alloc_double(N_Outermost_Split);
  unsigned n_box_l, n_box_r,n_r,n_l;
  double M,s;
  char var[100] = {'\0'};
  char par[100] = {'\0'};
  char val[100] = {'\0'};
  unsigned i;
  
  assert(GRT(C,0));
  assert(GRT(R_NS_l,0));
  assert(GRT(R_NS_r,0));
  
  n_l = (unsigned)GetParameterI("n_c");
  i   = (unsigned)GetParameterI("left_NS_n_c");
  if (i != INT_MAX) 	n_l = i;
  if (n_l == INT_MAX)   abortEr("n_l could not be set.");
  assert(n_l > 2);
  
  n_r = (unsigned)GetParameterI("n_c");
  i   = (unsigned)GetParameterI("right_NS_n_c");
  if (i != INT_MAX) 	n_r = i;
  if (n_r == INT_MAX)   abortEr("n_r could not be set.");
  assert(n_r > 2);
  
  /* making field of NS's radius and 
  // finding the max distance from the center of the star */
  NS_radii_BNS_Projective_grid(grid,r_max);
  R_max_l = r_max->R_max_l;
  R_max_r = r_max->R_max_r;
  
  box_size_l = 2*R_max_l/(n_l-1);
  box_size_r = 2*R_max_r/(n_r-1);
  
  O = C-R_NS_l-R_NS_r;
  if(!GRT(O,0))
    abortEr("The centers of neutron stars are too close.");
    
  O_l = O/2+R_NS_l;
  O_r = O/2+R_NS_r;
  
  M = GRT(O_l,O_r) ? O_l : O_r;
  s = O/M/2;
  R_Surr_l = sqrt(SQR(O_l)+s);
  R_Surr_r = sqrt(SQR(O_r)+s);
  assert(LSS(R_Surr_l-O_l,O/2));
  assert(LSS(R_Surr_r-O_r,O/2));
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R0[i] = GetParameterD_E(var);
    
    assert(GRT(R0[i],s));/* => R2_outermost <= R2_Surr */
    
    if (LSS(R0[i],C))
      abortEr("the radius of outermost patches must be greater than BNS distance.");
    
    if (i > 0)
      if (LSSEQL(R0[i],R0[i-1]))
        abortEr("The radius of outermost must be increasing.");
        
    R_outmost_l[i] = sqrt(SQR(O_l)+SQR(R0[i]));
    R_outmost_r[i] = sqrt(SQR(O_r)+SQR(R0[i]));
    
  }
  
  n_box_l = (unsigned)(3 + n_l/CONST);
  n_box_r = (unsigned)(3 + n_r/CONST);
  
  /* adding the results to the parameter data base */
  
  /* n_a, n_b, n_c */
  sprintf(par,"grid%u_left_centeral_box_n_abc",gn);
  sprintf(val,"%u",n_box_l);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_centeral_box_n_abc",gn);
  sprintf(val,"%u",n_box_r);
  add_parameter_string(par,val);
  
  /* size a,b,c */
  sprintf(par,"grid%u_left_centeral_box_size_a",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_b",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_c",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_right_centeral_box_size_a",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_right_centeral_box_size_b",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_right_centeral_box_size_c",gn);
  add_parameter_double(par,box_size_r);
  
  /* R2 surroundings */
  sprintf(par,"grid%u_left_NS_surrounding_R2",gn);
  add_parameter_double(par,R_Surr_l);
  
  sprintf(par,"grid%u_right_NS_surrounding_R2",gn);
  add_parameter_double(par,R_Surr_r);
  
  /* R1 and R2 outermost */
  for (i = 0; i < N_Outermost_Split; i++)
  {
    /* R1: */
    if (i == 0)
    {
      sprintf(par,"grid%u_left_outermost%u_R1",gn,i);
      add_parameter_double(par,R_Surr_l);
      
      sprintf(par,"grid%u_right_outermost%u_R1",gn,i);
      add_parameter_double(par,R_Surr_r);
    }
    else
    {
      sprintf(par,"grid%u_left_outermost%u_R1",gn,i);
      add_parameter_double(par,R_outmost_l[i-1]);
      
      sprintf(par,"grid%u_right_outermost%u_R1",gn,i);
      add_parameter_double(par,R_outmost_r[i-1]);
    }
    
    /* R2: */
    sprintf(par,"grid%u_left_outermost%u_R2",gn,i);
    add_parameter_double(par,R_outmost_l[i]);
    
    sprintf(par,"grid%u_right_outermost%u_R2",gn,i);
    add_parameter_double(par,R_outmost_r[i]);
  }
  
  /* assuming the center of left NS at (0,-O_l,0) */
  sprintf(par,"grid%u_left_NS_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_left_NS_center_b",gn);
  add_parameter_double(par,-O_l);
  
  sprintf(par,"grid%u_left_NS_center_c",gn);
  add_parameter_double(par,0.0);
  
  /* assuming the center of right NS at (0,O_r,0) */
  sprintf(par,"grid%u_right_NS_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_right_NS_center_b",gn);
  add_parameter_double(par,O_r);
  
  sprintf(par,"grid%u_right_NS_center_c",gn);
  add_parameter_double(par,0.0);
  
  free(R0);
  free(R_outmost_r);
  free(R_outmost_l);
}

/* calculating the main characteristic of grid for BNS Spherical grid */
static void characteristics_BNS_Spherical_grid(Grid_T *const grid)
{
  const unsigned gn   = grid->gn;
  const double C      = GetParameterD_E("BNS_Distance");
  const double R_NS_l = GetParameterD_E("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = GetParameterD_E("right_NS_radius");/* assuming perfect sphere */
  const unsigned N_Outermost_Split = (unsigned)GetParameterI_E("Number_of_Outermost_Split"); 
  double O,O_l,O_r,
         R_Surr_l,R_Surr_r,
         *R_outmost_l = alloc_double(N_Outermost_Split),
         *R_outmost_r = alloc_double(N_Outermost_Split),
         *R0 = alloc_double(N_Outermost_Split);
  double M,s;
  char var[100] = {'\0'};
  char par[100] = {'\0'};
  unsigned i;
  
  assert(GRT(C,0));
  assert(GRT(R_NS_l,0));
  assert(GRT(R_NS_r,0));
  
  /* making field of NS's radius and 
  // finding the max distance from the center of the star */
  NS_radii_BNS_Spherical_grid(grid,0);
  
  O = C-R_NS_l-R_NS_r;
  if(!GRT(O,0))
    abortEr("The centers of neutron stars are too close.");
    
  O_l = O/2+R_NS_l;
  O_r = O/2+R_NS_r;
  
  M = GRT(O_l,O_r) ? O_l : O_r;
  s = O/M/2;
  R_Surr_l = sqrt(SQR(O_l)+s);
  R_Surr_r = sqrt(SQR(O_r)+s);
  assert(LSS(R_Surr_l-O_l,O/2));
  assert(LSS(R_Surr_r-O_r,O/2));
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R0[i] = GetParameterD_E(var);
    
    assert(GRT(R0[i],s));/* => R2_outermost <= R2_Surr */
    
    if (LSS(R0[i],C))
      abortEr("the radius of outermost patches must be greater than BNS distance.");
    
    if (i > 0)
      if (LSSEQL(R0[i],R0[i-1]))
        abortEr("The radius of outermost must be increasing.");
        
    R_outmost_l[i] = sqrt(SQR(O_l)+SQR(R0[i]));
    R_outmost_r[i] = sqrt(SQR(O_r)+SQR(R0[i]));
    
  }
  
  /* adding the results to the parameter data base: */
  
  /* R2 surroundings */
  sprintf(par,"grid%u_left_NS_Surrounding_R2",gn);
  add_parameter_double(par,R_Surr_l);
  
  sprintf(par,"grid%u_right_NS_Surrounding_R2",gn);
  add_parameter_double(par,R_Surr_r);
  
  /* R1 and R2 outermost */
  for (i = 0; i < N_Outermost_Split; i++)
  {
    /* R1: */
    if (i == 0)
    {
      sprintf(par,"grid%u_left_outermost%u_R1",gn,i);
      add_parameter_double(par,R_Surr_l);
      
      sprintf(par,"grid%u_right_outermost%u_R1",gn,i);
      add_parameter_double(par,R_Surr_r);
    }
    else
    {
      sprintf(par,"grid%u_left_outermost%u_R1",gn,i);
      add_parameter_double(par,R_outmost_l[i-1]);
      
      sprintf(par,"grid%u_right_outermost%u_R1",gn,i);
      add_parameter_double(par,R_outmost_r[i-1]);
    }
    
    /* R2: */
    sprintf(par,"grid%u_left_outermost%u_R2",gn,i);
    add_parameter_double(par,R_outmost_l[i]);
    
    sprintf(par,"grid%u_right_outermost%u_R2",gn,i);
    add_parameter_double(par,R_outmost_r[i]);
  }
  
  /* assuming the center of left NS at (0,-O_l,0) */
  sprintf(par,"grid%u_left_NS_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_left_NS_center_b",gn);
  add_parameter_double(par,-O_l);
  
  sprintf(par,"grid%u_left_NS_center_c",gn);
  add_parameter_double(par,0.0);
  
  /* assuming the center of right NS at (0,O_r,0) */
  sprintf(par,"grid%u_right_NS_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_right_NS_center_b",gn);
  add_parameter_double(par,O_r);
  
  sprintf(par,"grid%u_right_NS_center_c",gn);
  add_parameter_double(par,0.0);
  
  free(R0);
  free(R_outmost_r);
  free(R_outmost_l);
}

/* making field of NS's inside (R1) and surface (R2) radius and 
// finding the max distance from the center of the star */
static void NS_radii_BNS_Projective_grid(Grid_T *const grid,void *vp)
{
  const double R_NS_l = GetParameterD_E("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = GetParameterD_E("right_NS_radius");/* assuming perfect sphere */
  double *R1_down,*R2_down;
  double *R1_up,*R2_up;
  struct R_max *const r_max = vp;
  char par[100] = {'\0'};
  unsigned N[3],n,i,j,N_total;
  
  /* left NS */
  r_max->R_max_l = R_NS_l;/* since this star is perfect sphere */
  
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("left_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("left_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("left_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface: */
  
  /* up side */
  R2_up = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R2_up[L(N,i,j,0)] = R_NS_l;
      
  sprintf(par,"grid%u_left_NS_R2_up",grid->gn);
  add_parameter_array(par,R2_up,N_total);
  /* end of up side */
  /* down side */
  R2_down = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R2_down[L(N,i,j,0)] = R_NS_l;
      
  sprintf(par,"grid%u_left_NS_R2_down",grid->gn);
  add_parameter_array(par,R2_down,N_total);
  /* end of down side */
  
  /* inside: */
  
  /* up side */
  R1_up = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R1_up[L(N,i,j,0)] = R2_up[L(N,i,j,0)]/N[2];
      
  sprintf(par,"grid%u_left_NS_R1_up",grid->gn);
  add_parameter_array(par,R1_up,N_total);
  free(R1_up);
  free(R2_up);
  /* end of up side */
  /* down side */
  R1_down = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R1_down[L(N,i,j,0)] = R2_down[L(N,i,j,0)]/N[2];
      
  sprintf(par,"grid%u_left_NS_R1_down",grid->gn);
  add_parameter_array(par,R1_down,N_total);
  free(R1_down);
  free(R2_down);
  /* end of down side */
  
  /* right NS */
  r_max->R_max_r = R_NS_r;/* since this star is perfect sphere */
  
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("right_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("right_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("right_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface: */
  
  /* up side */
  R2_up = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R2_up[L(N,i,j,0)] = R_NS_r;
      
  sprintf(par,"grid%u_right_NS_R2_up",grid->gn);
  add_parameter_array(par,R2_up,N_total);
  /* end of up side */
  /* down side */
  R2_down = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R2_down[L(N,i,j,0)] = R_NS_r;
      
  sprintf(par,"grid%u_right_NS_R2_down",grid->gn);
  add_parameter_array(par,R2_down,N_total);
  /* end of down side */
  
  /* inside: */
  
  /* up side */
  R1_up = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R1_up[L(N,i,j,0)] = R2_up[L(N,i,j,0)]/N[2];
      
  sprintf(par,"grid%u_right_NS_R1_up",grid->gn);
  add_parameter_array(par,R1_up,N_total);
  free(R1_up);
  free(R2_up);
  /* end of up side */
  /* down side */
  R1_down = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R1_down[L(N,i,j,0)] = R2_down[L(N,i,j,0)]/N[2];
      
  sprintf(par,"grid%u_right_NS_R1_down",grid->gn);
  add_parameter_array(par,R1_down,N_total);
  free(R1_down);
  free(R2_down);
  /* end of down side */
  
}

/* making field of NS's inside (R1) and surface (R2) radius and 
// finding the max distance from the center of the star */
static void NS_radii_BNS_Spherical_grid(Grid_T *const grid,void *vp)
{
  const double R_NS_l = GetParameterD_E("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = GetParameterD_E("right_NS_radius");/* assuming perfect sphere */
  double *R1,*R2;
  char par[100] = {'\0'};
  unsigned N[3],n,j,k,N_total;
  
  UNUSED(vp);
  
  /* left NS */
  
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("left_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("left_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("left_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface: R2 */
  R2 = alloc_double(N_total);
  for (j = 0; j < N[1]; ++j)
    for (k = 0; k < N[2]; ++k)
      R2[L(N,0,j,k)] = R_NS_l;
      
  sprintf(par,"grid%u_left_NS_R2",grid->gn);
  add_parameter_array(par,R2,N_total);
  
  /* inside: R1 */
  R1 = alloc_double(N_total);
  for (j = 0; j < N[1]; ++j)
    for (k = 0; k < N[2]; ++k)
      R1[L(N,0,j,k)] = 0;
      
  sprintf(par,"grid%u_left_NS_R1",grid->gn);
  add_parameter_array(par,R1,N_total);
  free(R1);
  free(R2);
  
  /* right NS */
  
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("right_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("right_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("right_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface: R2 */
  R2 = alloc_double(N_total);
  for (j = 0; j < N[1]; ++j)
    for (k = 0; k < N[2]; ++k)
      R2[L(N,0,j,k)] = R_NS_r;
      
  sprintf(par,"grid%u_right_NS_R2",grid->gn);
  add_parameter_array(par,R2,N_total);
  
  /* inside: R1*/
  R1 = alloc_double(N_total);
  for (j = 0; j < N[1]; ++j)
    for (k = 0; k < N[2]; ++k)
      R1[L(N,0,j,k)] = 0;

  sprintf(par,"grid%u_right_NS_R1",grid->gn);
  add_parameter_array(par,R1,N_total);
  free(R1);
  free(R2);
}

/* making  NS's surface function */
static void NS_surface_BNS_CubedSpherical_grid(Grid_T *const grid)
{
  const double R_NS_l = GetParameterD_E("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = GetParameterD_E("right_NS_radius");/* assuming perfect sphere */
  double *R;
  char par[100] = {'\0'};
  unsigned N[3],n,i,j,N_total;
  
  /* left NS */
  
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("left_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("left_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("left_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface */
  R = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R[L(N,i,j,0)] = R_NS_l;
      
  sprintf(par,"grid%u_left_NS_surface_function_up",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_left_NS_surface_function_down",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_left_NS_surface_function_back",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_left_NS_surface_function_front",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_left_NS_surface_function_left",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_left_NS_surface_function_right",grid->gn);
  add_parameter_array(par,R,N_total);
  
  free(R);
  
  /* right NS */
  
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("right_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("right_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("right_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface */
  R = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R[L(N,i,j,0)] = R_NS_r;
      
  sprintf(par,"grid%u_right_NS_surface_function_up",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_right_NS_surface_function_down",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_right_NS_surface_function_back",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_right_NS_surface_function_front",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_right_NS_surface_function_left",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_right_NS_surface_function_right",grid->gn);
  add_parameter_array(par,R,N_total);
  
  free(R);
}

