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
    
    else
      Error0("No action for such coordinate.\n");
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
      IsNull(patch->JacobianT);
      make_JacobianT_Cartesian_coord(patch);
    }
    /* if coord is Cubed Spherical */
    else if (patch->coordsys == CubedSpherical)
    {
      patch->JacobianT = calloc(1,sizeof(*patch->JacobianT));
      IsNull(patch->JacobianT);
      make_JacobianT_CubedSpherical_coord(patch);
    }
    else
      Error0("No job for such coordinate.\n");
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
    Error0("There is no such COLLOCATION.\n");
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
    Error0("There is no such collocation.\n");
  
  return x;
}

/* dq2/dq1 Jacobian transformaion for coordinates,
// where q? could be any N0,N1,N2,x,y,z,a,b,c.
// list of q?_e:
// =============
//
// _N0_ for normalized 0-coord like [-1,1] that is used in an expansion
// _N1_ for normalized 1-coord like [-1,1] that is used in an expansion
// _N2_ for normalized 2-coord like [-1,1] that is used in an expansion
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
    j = dq_dN(patch,q2_e,q1_e,p);
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
    Error0(INCOMPLETE_FUNC);
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
          Error0("It should not reach here!");
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

/* Jacobian transformation for dq/dN?. since we know dN/dq, 
// we can use the inverse property of J^-1 = adj(J)/det(J) to find dq/dN.
// ->return value: dq/dN? */
static double dq_dN(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  double j = 0;
  
  if (q2_e == _x_ || q2_e == _y_ || q2_e == _z_)
  {
    const double a00 = dq2_dq1(patch,_N0_,_x_,p);
    const double a01 = dq2_dq1(patch,_N0_,_y_,p);
    const double a02 = dq2_dq1(patch,_N0_,_z_,p);
    const double a10 = dq2_dq1(patch,_N1_,_x_,p);
    const double a11 = dq2_dq1(patch,_N1_,_y_,p);
    const double a12 = dq2_dq1(patch,_N1_,_z_,p);
    const double a20 = dq2_dq1(patch,_N2_,_x_,p);
    const double a21 = dq2_dq1(patch,_N2_,_y_,p);
    const double a22 = dq2_dq1(patch,_N2_,_z_,p);
    const double det = a00 *a11 *a22  - a00 *a12 *a21  -
                       a01 *a10 *a22  + a01 *a12 *a20  +
                       a02 *a10 *a21  - a02 *a11 *a20  ;/* determinate */
    const unsigned a = q2_e%3;
    const unsigned b = q1_e%3;
    assert(!EQL(det,0));
    
    if (a == 0)
    {
      switch(b)
      {
        case 0:
          j = ADJ_00/det;
        break;
        case 1:
          j = ADJ_01/det;
        break;
        case 2:
          j = ADJ_02/det;
        break;
        default:
          Error0(NO_OPTION);
      }
    }
    else if(a == 1)
    {
      switch(b)
      {
        case 0:
          j = ADJ_10/det;
        break;
        case 1:
          j = ADJ_11/det;
        break;
        case 2:
          j = ADJ_12/det;
        break;
        default:
          Error0(NO_OPTION);
      }
    }
    else if (a == 2)
    {
      switch(b)
      {
        case 0:
          j = ADJ_20/det;
        break;
        case 1:
          j = ADJ_21/det;
        break;
        case 2:
          j = ADJ_22/det;
        break;
        default:
          Error0(NO_OPTION);
      }
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(INCOMPLETE_FUNC);
  
  return j;
}

/* testing dq_dN function */
void test_dq_dN(Grid_T *const grid)
{
  double a00,a01,a02,
         a10,a11,a12,
         a20,a21,a22,
         b00,b01,b02,
         b10,b11,b12,
         b20,b21,b22,
         c00,c01,c02,
         c10,c11,c12,
         c20,c21,c22;
  unsigned p,ijk,nn;
  
  /* testing dx/dN */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *const patch = grid->patch[p];
    
    nn = patch->nn;
    for (ijk = 0; ijk < nn; ++ijk)
    {
      a00 = dq2_dq1(patch,_N0_,_x_,ijk);
      a01 = dq2_dq1(patch,_N0_,_y_,ijk);
      a02 = dq2_dq1(patch,_N0_,_z_,ijk);
      a10 = dq2_dq1(patch,_N1_,_x_,ijk);
      a11 = dq2_dq1(patch,_N1_,_y_,ijk);
      a12 = dq2_dq1(patch,_N1_,_z_,ijk);
      a20 = dq2_dq1(patch,_N2_,_x_,ijk);
      a21 = dq2_dq1(patch,_N2_,_y_,ijk);
      a22 = dq2_dq1(patch,_N2_,_z_,ijk);
      
      b00 = dq2_dq1(patch,_x_,_N0_,ijk);
      b01 = dq2_dq1(patch,_x_,_N1_,ijk);
      b02 = dq2_dq1(patch,_x_,_N2_,ijk);
      b10 = dq2_dq1(patch,_y_,_N0_,ijk);
      b11 = dq2_dq1(patch,_y_,_N1_,ijk);
      b12 = dq2_dq1(patch,_y_,_N2_,ijk);
      b20 = dq2_dq1(patch,_z_,_N0_,ijk);
      b21 = dq2_dq1(patch,_z_,_N1_,ijk);
      b22 = dq2_dq1(patch,_z_,_N2_,ijk);
      
      c00 = a00*b00 + a01*b10 + a02*b20;
      c01 = a00*b01 + a01*b11 + a02*b21;
      c02 = a00*b02 + a01*b12 + a02*b22;
      c10 = a10*b00 + a11*b10 + a12*b20;
      c11 = a10*b01 + a11*b11 + a12*b21;
      c12 = a10*b02 + a11*b12 + a12*b22;
      c20 = a20*b00 + a21*b10 + a22*b20;
      c21 = a20*b01 + a21*b11 + a22*b21;
      c22 = a20*b02 + a21*b12 + a22*b22;
      
      if(!EQL(c00,1))  Error0("dx_dN is not correct!\n");
      if(!EQL(c01,0))  Error0("dx_dN is not correct!\n");
      if(!EQL(c02,0))  Error0("dx_dN is not correct!\n");
      if(!EQL(c10,0))  Error0("dx_dN is not correct!\n");
      if(!EQL(c11,1))  Error0("dx_dN is not correct!\n");
      if(!EQL(c12,0))  Error0("dx_dN is not correct!\n");
      if(!EQL(c20,0))  Error0("dx_dN is not correct!\n");
      if(!EQL(c21,0))  Error0("dx_dN is not correct!\n");
      if(!EQL(c22,1))  Error0("dx_dN is not correct!\n");
    }
  }
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
    Error0("This direction at the specified patch doesn't use first kind Chebyshev,\n"
      "and Chebyshev Extrema collocation which this function needs.\n");
    
  const double a = -patch->max[dir]+patch->min[dir];
  const double b =  patch->max[dir]+patch->min[dir];
  double cheb = (X*2-b)/a;
  
  if (EQL(cheb,1.0))
    cheb = 1;
  else if (EQL(cheb,-1.0))
    cheb = -1;
    
  return cheb;
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
    Error1("There is no such %s.\n",needle);
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
     Error0("Invalid entry.");
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
     Error0("Invalid entry.");
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
     Error0("Invalid entry.");
  }/* end of else if (q2_e == _c_) */
  
  else
    Error0("Invalid entry.");
  
  return dA_da;  
}

/* calculating the main characteristics of grid needed for setting up
// the patches, like lengths, radii and etc for tutorial purposes. */
void grid_characteristics_example(Grid_T *const grid)
{
  const char *kind;
  
  /* finding the kind of grid */
  kind = Pgets("grid_kind");
  
  if (strcmp_i(kind,"Cartesian_grid"))
    characteristics_Cartesian_grid_eg(grid);
  
  else if (strcmp_i(kind,"BNS_Spherical_grid"))
    characteristics_BNS_Spherical_grid_eg(grid); 
    
  else if (strcmp_i(kind,"BNS_CubedSpherical_grid"))
    characteristics_BNS_CS_grid_eg(grid); 
    
  else if (strcmp_i(kind,"BBN_CubedSpherical_grid"))
    characteristics_BBN_CS_grid_eg(grid); 
   
  else if (strcmp_i(kind,"BBN_Split_CubedSpherical_grid"))
    characteristics_BBN_Split_CS_grid_eg(grid); 
    
  else
    Error1("There is no such %s grid kind.\n",grid->kind);
  
}

/* calculating the main characteristic of grid for Cartesian grid */
static void characteristics_Cartesian_grid_eg(Grid_T *const grid)
{
  /* this type of grid is so simple; nothing to calculate. */
  grid->kind = dup_s(Pgets("grid_kind"));
}

/* calculating the main characteristic of grid for BBN_CubedSpherical grid */
static void characteristics_BBN_CS_grid_eg(Grid_T *const grid)
{
  /* calculate the characteristics of this grid */
  const unsigned gn   = grid->gn;
  const double C      = Pgetd("BH_NS_separation");
  const double R_NS_l = Pgetd("NS_radius"),
               bh_m   = Pgetd("BH_mass"),
               bh_chi = Pgetd("BH_dimensionless_spin"),
               R_BH_r = bh_m*(1+sqrt(1-Pow2(bh_chi)));
               
  double box_size_l;
  const unsigned N_Outermost_Split = (unsigned)Pgeti("Number_of_Outermost_Split"); 
  double *R_outermost = calloc(N_Outermost_Split,sizeof(*R_outermost));
  unsigned nlb[3]/*left box*/,n;
  char var[100] = {'\0'};
  char par[100] = {'\0'};
  char val[100] = {'\0'};
  const char *kind;
  unsigned i;
  
  /* finding the kind of grid */
  kind = Pgets("grid_kind");
  grid->kind = dup_s(kind);
  
  assert(GRT(C,0));
  assert(GRT(R_NS_l,0));
  assert(GRT(R_BH_r,0));
  assert(LSS(2*R_NS_l,C));
  assert(LSS(2*R_BH_r,C));
  
  /* making NS and BH surfaces function */
  NS_BH_surface_CS_grid_eg(grid,R_NS_l,R_BH_r,bh_m*bh_chi);
  
  box_size_l = Pgetd("left_central_box_length_ratio")*R_NS_l;
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R_outermost[i] = Pgetd(var);
    
    if (LSS(R_outermost[i],2*C))
      Error0("the radius of outermost patches must be greater than twice of BNS distance.");
    
    if (i > 0)
      if (LSSEQL(R_outermost[i],R_outermost[i-1]))
        Error0("The radius of outermost must be increasing.");
    
  }
  
  /* filling n */
  
  /* left box */
  nlb[0] = (unsigned)PgetiEZ("n_a");
  nlb[1] = (unsigned)PgetiEZ("n_b");
  nlb[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("left_NS_n_a");
  if (n != INT_MAX)   nlb[0] = n;
  n = (unsigned)PgetiEZ("left_NS_n_b");
  if (n != INT_MAX)   nlb[1] = n;
  n = (unsigned)PgetiEZ("left_NS_n_c");
  if (n != INT_MAX)   nlb[2] = n;
    
  if(nlb[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(nlb[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(nlb[2] == INT_MAX)
    Error0("n_c could not be set.\n");
  
  /* adding the results to the parameter data base */
  
  /* n_a, n_b, n_c */
  /* left box */
  sprintf(par,"grid%u_left_central_box_n_a",gn);
  sprintf(val,"%u",nlb[0]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_left_central_box_n_b",gn);
  sprintf(val,"%u",nlb[1]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_left_central_box_n_c",gn);
  sprintf(val,"%u",nlb[2]);
  add_parameter_string(par,val);
  
  /* size a,b,c */
  sprintf(par,"grid%u_left_central_box_size_a",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_central_box_size_b",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_central_box_size_c",gn);
  add_parameter_double(par,box_size_l);
  
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
  
  /* assuming the center of right BH at (0,C/2,0) */
  sprintf(par,"grid%u_right_BH_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_right_BH_center_b",gn);
  add_parameter_double(par,C/2);
  
  sprintf(par,"grid%u_right_BH_center_c",gn);
  add_parameter_double(par,0.0);
  
  free(R_outermost);
}

/* example of calculating the main characteristic of grid for 
// BBN_Split_CubedSpherical grid. 
// note: NS and BH sufaces are required to be given in Ylm expansions. */
static void characteristics_BBN_Split_CS_grid_eg(Grid_T *const grid)
{
  /* calculate the characteristics of this grid */
  Grid_Char_T grid_char[1] = {0};/* grid characteristics */
  const unsigned ns = 0, bh = 1;
  const unsigned lmax   = 5;
  const unsigned Ntheta = Ntheta_Ylm(lmax);
  const unsigned Nphi   = Nphi_Ylm(lmax);
  const unsigned Ntot   = Ntotal_Ylm(lmax);
  const double C      = Pgetd("BH_NS_separation");
  const double R_NS_l = Pgetd("NS_radius"),
               bh_m   = Pgetd("BH_mass"),
               bh_chi = Pgetd("BH_dimensionless_spin"),
               R_BH_r = bh_m*(1+sqrt(1-Pow2(bh_chi)));
  double *rns = alloc_double(Ntot);/* surface function r = r(th,ph). */
  double *rbh = alloc_double(Ntot);/* surface function r = r(th,ph). */
  double *reClm_rns = alloc_ClmYlm(lmax),
         *imClm_rns = alloc_ClmYlm(lmax);
  double *reClm_rbh = alloc_ClmYlm(lmax),
         *imClm_rbh = alloc_ClmYlm(lmax);
  double box_size_l,box_size_r;
  const char *kind;
  unsigned ij;
  
  /* set surface functions (required in Ylm) */
  /* initialize tables */
  init_Legendre_root_function();
  for (ij = 0; ij < Ntot; ++ij)
  {
    rns[ij] = R_NS_l;
    rbh[ij] = R_BH_r;
  }
  /* calculating coeffs */
  get_Ylm_coeffs(reClm_rns,imClm_rns,rns,Ntheta,Nphi,lmax);
  get_Ylm_coeffs(reClm_rbh,imClm_rbh,rbh,Ntheta,Nphi,lmax);
  
  /* finding the kind of grid */
  kind = Pgets("grid_kind");
  grid->kind = dup_s(kind);
  
  assert(C > 0);
  assert(R_NS_l > 0);
  assert(R_BH_r > 0);
  assert(LSS(2*R_NS_l,C));
  assert(LSS(2*R_BH_r,C));
  
  box_size_r = box_size_l = Pgetd("left_central_box_length_ratio")*R_NS_l;
  
  /* set char of grid */
  grid_char->grid = grid;
  grid_char->S    = C;
  /* NS */
  grid_char->params[ns]->obj    = "NS";
  grid_char->params[ns]->dir    = "left";
  grid_char->params[ns]->relClm = reClm_rns;
  grid_char->params[ns]->imgClm = imClm_rns;
  grid_char->params[ns]->lmax   = lmax;
  grid_char->params[ns]->r_min  = R_NS_l;
  grid_char->params[ns]->r_max  = R_NS_l;
  grid_char->params[ns]->l      = box_size_l;
  grid_char->params[ns]->w      = box_size_l;
  grid_char->params[ns]->h      = box_size_l;
  /* BH */
  grid_char->params[bh]->obj    = "BH";
  grid_char->params[bh]->dir    = "right";
  grid_char->params[bh]->relClm = reClm_rbh;
  grid_char->params[bh]->imgClm = imClm_rbh;
  grid_char->params[bh]->lmax   = lmax;
  grid_char->params[bh]->r_min  = R_BH_r;
  grid_char->params[bh]->r_max  = R_BH_r;
  grid_char->params[bh]->l      = box_size_r;
  grid_char->params[bh]->w      = box_size_r;
  grid_char->params[bh]->h      = box_size_r;
  
  /* set number of splits, points in each directions,
  // surface functions etc. */
  set_params_split_CS(grid_char);
  
  /* free */
  _free(rns);
  _free(rbh);
  _free(reClm_rns);
  _free(imClm_rns);
  _free(reClm_rbh);
  _free(imClm_rbh);
  
}

/* making  NS and BH surfaces function */
static void NS_BH_surface_CS_grid_eg(Grid_T *const grid,const double R_NS_l,const double R_BH_r,const double a_BH)
{
  double *R;
  char par[100] = {'\0'};
  unsigned N[3],n,i,j,k,N_total;
  Patch_T patch[1] = {0};
  struct Collocation_s coll_s[2] = {0};
  double X[2],r;
  
  /* left NS */
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = 0;

  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;

  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
  
  /* left NS */
  
  /* filling N */
  N[0] = (unsigned)PgetiEZ("n_a");
  N[1] = (unsigned)PgetiEZ("n_b");
  N[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("left_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)PgetiEZ("left_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)PgetiEZ("left_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    Error0("n_c could not be set.\n");
  
  patch->n[0] = N[0];
  patch->n[1] = N[1];
  patch->n[2] = N[2];
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
    
  N_total = N[0]*N[1]*N[2];
  
  /* NS surface */
  
  R = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = R_NS_l;
        
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
  
  /* right BH */
  
  /* filling N */
  N[0] = (unsigned)PgetiEZ("n_a");
  N[1] = (unsigned)PgetiEZ("n_b");
  N[2] = (unsigned)PgetiEZ("n_c");
  
  /* check for override */
  n = (unsigned)PgetiEZ("right_BH_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)PgetiEZ("right_BH_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)PgetiEZ("right_BH_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    Error0("n_c could not be set.\n");
   
  patch->n[0] = N[0];
  patch->n[1] = N[1];
  patch->n[2] = N[2];
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  
  N_total = N[0]*N[1]*N[2];
  
  R = alloc_double(N_total);
  
  /* surface up and down */
  for (i = 0; i < N[0]; ++i)
  {
    X[0] = point_value(i,&coll_s[0]);
    for (j = 0; j < N[1]; ++j)
    {
      X[1] = point_value(j,&coll_s[1]);
      r = sqrt(
               (1+Pow2(X[0])+Pow2(X[1]))/
               ((Pow2(X[0])+Pow2(X[1]))/(Pow2(R_BH_r)+Pow2(a_BH)) + 1/Pow2(R_BH_r))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_right_BH_surface_function_up",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_right_BH_surface_function_down",grid->gn);
  add_parameter_array(par,R,N_total);
  
  /* surface back */
  for (i = 0; i < N[0]; ++i)
  {
    X[0] = point_value(i,&coll_s[0]);/* a = z/x */
    for (j = 0; j < N[1]; ++j)
    {
      X[1] = point_value(j,&coll_s[1]);/* b = y/x */
      r = sqrt(
               (1+Pow2(X[0])+Pow2(X[1]))/
               (((1+Pow2(X[1])))/(Pow2(R_BH_r)+Pow2(a_BH)) + Pow2(X[0])/Pow2(R_BH_r))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_right_BH_surface_function_back",grid->gn);
  add_parameter_array(par,R,N_total);
  
  /* surface front */
  for (i = 0; i < N[0]; ++i)
  {
    X[0] = point_value(i,&coll_s[0]);/* a = y/x */
    for (j = 0; j < N[1]; ++j)
    {
      X[1] = point_value(j,&coll_s[1]);/* b = z/x */
      r = sqrt(
               (1+Pow2(X[0])+Pow2(X[1]))/
               (((1+Pow2(X[0])))/(Pow2(R_BH_r)+Pow2(a_BH)) + Pow2(X[1])/Pow2(R_BH_r))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_right_BH_surface_function_front",grid->gn);
  add_parameter_array(par,R,N_total);
  
  /* surface left */
  for (i = 0; i < N[0]; ++i)
  {
    X[0] = point_value(i,&coll_s[0]);/* a = x/y */
    for (j = 0; j < N[1]; ++j)
    {
      X[1] = point_value(j,&coll_s[1]);/* b = z/y */
      r = sqrt(
               (1+Pow2(X[0])+Pow2(X[1]))/
               (((1+Pow2(X[0])))/(Pow2(R_BH_r)+Pow2(a_BH)) + Pow2(X[1])/Pow2(R_BH_r))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_right_BH_surface_function_left",grid->gn);
  add_parameter_array(par,R,N_total);
  
  /* surface right */
  for (i = 0; i < N[0]; ++i)
  {
    X[0] = point_value(i,&coll_s[0]);/* a = z/y */
    for (j = 0; j < N[1]; ++j)
    {
      X[1] = point_value(j,&coll_s[1]);/* b = x/y */
      r = sqrt(
               (1+Pow2(X[0])+Pow2(X[1]))/
               (((1+Pow2(X[1])))/(Pow2(R_BH_r)+Pow2(a_BH)) + Pow2(X[0])/Pow2(R_BH_r))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_right_BH_surface_function_right",grid->gn);
  add_parameter_array(par,R,N_total);
  
  free(R);
}

/* calculating the main characteristic of grid for BNS_CubedSpherical grid */
static void characteristics_BNS_CS_grid_eg(Grid_T *const grid)
{
  const unsigned gn   = grid->gn;
  const double C      = Pgetd("BNS_Distance");
  const double R_NS_l = Pgetd("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = Pgetd("right_NS_radius");/* assuming perfect sphere */
  double box_size_l,box_size_r;
  const unsigned N_Outermost_Split = (unsigned)Pgeti("Number_of_Outermost_Split"); 
  double *R_outermost = calloc(N_Outermost_Split,sizeof(*R_outermost));
  unsigned nlb[3]/*left box*/, nrb[3]/*right box*/,n;
  char var[100] = {'\0'};
  char par[100] = {'\0'};
  char val[100] = {'\0'};
  unsigned i;
  
  assert(GRT(C,0));
  assert(GRT(R_NS_l,0));
  assert(GRT(R_NS_r,0));
  assert(LSS(2*R_NS_l,C));
  assert(LSS(2*R_NS_r,C));
  
  /* making NS's surface function */
  NS_surface_BNS_CS_grid_eg(grid);
  
  box_size_l = Pgetd("left_central_box_length_ratio")*R_NS_l;
  box_size_r = Pgetd("right_central_box_length_ratio")*R_NS_r;
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R_outermost[i] = Pgetd(var);
    
    if (LSS(R_outermost[i],2*C))
      Error0("the radius of outermost patches must be greater than twice of BNS distance.");
    
    if (i > 0)
      if (LSSEQL(R_outermost[i],R_outermost[i-1]))
        Error0("The radius of outermost must be increasing.");
    
  }
  
  /* filling n */
  
  /* left box */
  nlb[0] = (unsigned)PgetiEZ("n_a");
  nlb[1] = (unsigned)PgetiEZ("n_b");
  nlb[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("left_NS_n_a");
  if (n != INT_MAX)   nlb[0] = n;
  n = (unsigned)PgetiEZ("left_NS_n_b");
  if (n != INT_MAX)   nlb[1] = n;
  n = (unsigned)PgetiEZ("left_NS_n_c");
  if (n != INT_MAX)   nlb[2] = n;
    
  if(nlb[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(nlb[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(nlb[2] == INT_MAX)
    Error0("n_c could not be set.\n");
  
  /* right box */
  nrb[0] = (unsigned)PgetiEZ("n_a");
  nrb[1] = (unsigned)PgetiEZ("n_b");
  nrb[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("right_NS_n_a");
  if (n != INT_MAX)   nrb[0] = n;
  n = (unsigned)PgetiEZ("right_NS_n_b");
  if (n != INT_MAX)   nrb[1] = n;
  n = (unsigned)PgetiEZ("right_NS_n_c");
  if (n != INT_MAX)   nrb[2] = n;
    
  if(nrb[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(nrb[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(nrb[2] == INT_MAX)
    Error0("n_c could not be set.\n");
  
  /* adding the results to the parameter data base */
  
  /* n_a, n_b, n_c */
  /* left box */
  sprintf(par,"grid%u_left_central_box_n_a",nlb[0]);
  sprintf(val,"%u",nlb[0]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_left_central_box_n_b",nlb[1]);
  sprintf(val,"%u",nlb[1]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_left_central_box_n_c",nlb[2]);
  sprintf(val,"%u",nlb[2]);
  add_parameter_string(par,val);
  
  /* right box */
  sprintf(par,"grid%u_right_central_box_n_a",nrb[0]);
  sprintf(val,"%u",nrb[0]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_central_box_n_b",nrb[1]);
  sprintf(val,"%u",nrb[1]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_central_box_n_c",nrb[2]);
  sprintf(val,"%u",nrb[2]);
  add_parameter_string(par,val);
  
  /* size a,b,c */
  sprintf(par,"grid%u_left_central_box_size_a",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_central_box_size_b",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_central_box_size_c",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_right_central_box_size_a",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_right_central_box_size_b",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_right_central_box_size_c",gn);
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

/* calculating the main characteristic of grid for BNS Spherical grid */
static void characteristics_BNS_Spherical_grid_eg(Grid_T *const grid)
{
  const unsigned gn   = grid->gn;
  const double C      = Pgetd("BNS_Distance");
  const double R_NS_l = Pgetd("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = Pgetd("right_NS_radius");/* assuming perfect sphere */
  const unsigned N_Outermost_Split = (unsigned)Pgeti("Number_of_Outermost_Split"); 
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
  NS_radii_BNS_Spherical_grid_eg(grid,0);
  
  O = C-R_NS_l-R_NS_r;
  if(!GRT(O,0))
    Error0("The centers of neutron stars are too close.");
    
  O_l = O/2+R_NS_l;
  O_r = O/2+R_NS_r;
  
  M = GRT(O_l,O_r) ? O_l : O_r;
  s = O/M/2;
  R_Surr_l = sqrt(Pow2(O_l)+s);
  R_Surr_r = sqrt(Pow2(O_r)+s);
  assert(LSS(R_Surr_l-O_l,O/2));
  assert(LSS(R_Surr_r-O_r,O/2));
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R0[i] = Pgetd(var);
    
    assert(GRT(R0[i],s));/* => R2_outermost <= R2_Surr */
    
    if (LSS(R0[i],C))
      Error0("the radius of outermost patches must be greater than BNS distance.");
    
    if (i > 0)
      if (LSSEQL(R0[i],R0[i-1]))
        Error0("The radius of outermost must be increasing.");
        
    R_outmost_l[i] = sqrt(Pow2(O_l)+Pow2(R0[i]));
    R_outmost_r[i] = sqrt(Pow2(O_r)+Pow2(R0[i]));
    
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
static void NS_radii_BNS_Spherical_grid_eg(Grid_T *const grid,void *vp)
{
  const double R_NS_l = Pgetd("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = Pgetd("right_NS_radius");/* assuming perfect sphere */
  double *R1,*R2;
  char par[100] = {'\0'};
  unsigned N[3],n,ijk,N_total;
  
  UNUSED(vp);
  
  /* left NS */
  
  /* filling N */
  N[0] = (unsigned)PgetiEZ("n_a");
  N[1] = (unsigned)PgetiEZ("n_b");
  N[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("left_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)PgetiEZ("left_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)PgetiEZ("left_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    Error0("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface: R2 */
  R2 = alloc_double(N_total);
  for (ijk = 0; ijk < N_total; ++ijk)
      R2[ijk] = R_NS_l;
      
  sprintf(par,"grid%u_left_NS_R2",grid->gn);
  add_parameter_array(par,R2,N_total);
  
  /* inside: R1 */
  R1 = alloc_double(N_total);
  for (ijk = 0; ijk < N_total; ++ijk)
      R1[ijk] = 0;
      
  sprintf(par,"grid%u_left_NS_R1",grid->gn);
  add_parameter_array(par,R1,N_total);
  free(R1);
  free(R2);
  
  /* right NS */
  
  /* filling N */
  N[0] = (unsigned)PgetiEZ("n_a");
  N[1] = (unsigned)PgetiEZ("n_b");
  N[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("right_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)PgetiEZ("right_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)PgetiEZ("right_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    Error0("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface: R2 */
  R2 = alloc_double(N_total);
  for (ijk = 0; ijk < N_total; ++ijk)
    R2[ijk] = R_NS_r;
      
  sprintf(par,"grid%u_right_NS_R2",grid->gn);
  add_parameter_array(par,R2,N_total);
  
  /* inside: R1*/
  R1 = alloc_double(N_total);
  for (ijk = 0; ijk < N_total; ++ijk)
      R1[ijk] = 0;

  sprintf(par,"grid%u_right_NS_R1",grid->gn);
  add_parameter_array(par,R1,N_total);
  free(R1);
  free(R2);
}

/* making  NS's surface function */
static void NS_surface_BNS_CS_grid_eg(Grid_T *const grid)
{
  const double R_NS_l = Pgetd("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = Pgetd("right_NS_radius");/* assuming perfect sphere */
  double *R;
  char par[100] = {'\0'};
  unsigned N[3],n,i,j,k,N_total;
  
  /* left NS */
  
  /* filling N */
  N[0] = (unsigned)PgetiEZ("n_a");
  N[1] = (unsigned)PgetiEZ("n_b");
  N[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("left_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)PgetiEZ("left_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)PgetiEZ("left_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    Error0("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface */
  R = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = R_NS_l;
      
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
  N[0] = (unsigned)PgetiEZ("n_a");
  N[1] = (unsigned)PgetiEZ("n_b");
  N[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("right_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)PgetiEZ("right_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)PgetiEZ("right_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    Error0("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface */
  R = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = R_NS_r;
      
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

