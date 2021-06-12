/*
// Alireza Rashti
// June 2018
*/

#include "coordinate_utilities.h"


/* find the point in patches given by needle and fill
// the needle with the answers.
// for more information about needle look at the typede_data.h in 
// Core folder. for an example to how make needle look at grid.c and
// search for point_finder.
// note: don't forget to free memory at the end of the day for needle
// using free_needle.
*/
void point_finder(Needle_T *const needle)
{
  Flag_T flg = NO;
  Uint i,j;
  
  /* look inside guess patches */
  if (needle->Ng != 0)
  {
    find(needle,GUESS);
    
    /* if guess didn't work seek in the rest patches */
    if (needle->Nans == 0)
    {
      free(needle->in);
      needle->Nin = 0;
      
      FOR_ALL(i,needle->grid->patch)
      {
        flg = NO;
        /* look if this patch has been already found to exclude it */
        for (j = 0; j < needle->Ng; j++)
          if (needle->guess[j] == needle->grid->patch[i]->pn)
          {
            flg = YES;
            break;
          }
        if (flg == NO)
          needle_in(needle,needle->grid->patch[i]);
      }
      
      find(needle,FORCE_IN);
    }/* end of if (needle->Nans == 0) */
    
  }/* end of if (needle->Ng != 0) */
  
  /* look inside include patches */
  else if (needle->Nin != 0)
  {
    find(needle,FORCE_IN);
  }/* end of if (needle->Nin > 0) */
  /* find in all patched exluding needle->ex */
  else if (needle->Nex != 0)
  {
    FOR_ALL(i,needle->grid->patch)
    {
      flg = NO;
      for (j = 0; j < needle->Nex; j++)
      {
        if (needle->ex[j] == needle->grid->patch[i]->pn)
        {
          flg = YES;
          break;
        }
      }
      if (flg == NO)
        needle_in(needle,needle->grid->patch[i]);
    }
    
    find(needle,FORCE_IN);
  }/* end of if (needle->Nex > 0) */
  
}

/* adding a patch to needle->ex */
void needle_ex(Needle_T *const needle,const Patch_T *const patch)
{
  assert(needle);
  Uint i;
  
  i = 0;
  while(i < needle->Nex)
  {
    if (needle->ex[i] == patch->pn)  return;
    i++;
  }
  
  needle->ex = 
    realloc(needle->ex,(needle->Nex+1)*sizeof(*needle->ex));
  IsNull(needle->ex);
  
  needle->ex[needle->Nex] = patch->pn;
  needle->Nex++;
}

/* adding a patch to needle->in */
void needle_in(Needle_T *const needle,const Patch_T *const patch)
{
  assert(needle);
  Uint i;
  
  i = 0;
  while(i < needle->Nin)
  {
    if (needle->in[i] == patch->pn)  return;
    i++;
  }
  
  needle->in = 
    realloc(needle->in,(needle->Nin+1)*sizeof(*needle->in));
  IsNull(needle->in);
  
  needle->in[needle->Nin] = patch->pn;
  needle->Nin++;
}

/* adding a patch to needle->guess */
void needle_guess(Needle_T *const needle,const Patch_T *const patch)
{
  assert(needle);
  Uint i;
  
  i = 0;
  while(i < needle->Ng)
  {
    if (needle->guess[i] == patch->pn)  return;
    i++;
  }
  
  needle->guess = 
    realloc(needle->guess,(needle->Ng+1)*sizeof(*needle->guess));
  IsNull(needle->guess);
  
  needle->guess[needle->Ng] = patch->pn;
  needle->Ng++;
}

/* adding a patch to needle->ans */
void needle_ans(Needle_T *const needle,const Patch_T *const patch)
{
  assert(needle);
  Uint i;
  
  i = 0;
  while(i < needle->Nans)
  {
    if (needle->ans[i] == patch->pn)
      Error0("This point has been found twice in a same patch.\n"
      "Apparently, some part of point_finder is wrong or the needle"
      "has not been initialized correctly.\n");
    i++;
  }
  
  needle->ans = 
    realloc(needle->ans,(needle->Nans+1)*sizeof(*needle->ans));
  IsNull(needle->ans);
  
  needle->ans[needle->Nans] = patch->pn;
  needle->Nans++;
}

/* find for point in designated patches inside needle based on mode */
static void find(Needle_T *const needle,Mode_T mode)
{
  Uint *p = 0, np = UINT_MAX;
  Uint i;
  
  if (mode == GUESS)
  {
    p = needle->guess;
    np = needle->Ng;
  }
  else if (mode == FORCE_IN)
  {
    p = needle->in;
    np = needle->Nin;
  }
  else
    Error0("There is no such mode.\n");
  
  for (i = 0; i < np; i++)
  {
    double X[3];
    Patch_T *patch = needle->grid->patch[p[i]];
    int a;
    
    a = X_of_x_precision(X,needle->x,patch,needle->precision_factor);
    
    if (a)
      needle_ans(needle,patch);
      
  }
}

/* find point X(general coords) correspond to x(Cartesian coords) 
// in the given patch. one can increase precision by the precision_factor.
// ->return value 1 if it is successful, otherwise 0. */
int X_of_x_precision(double *const X,const double *const x,const Patch_T *const patch,const double precision_factor)
{
  int r = 0;
  
  if (patch->coordsys == Cartesian)
    r = X_of_x_Cartesian_coord(X,x,patch);
  else if (patch->coordsys == CubedSpherical)
    r = X_of_x_CS_coord(X,x,patch,precision_factor,1);
  else
    Error0("No finder for this coordinate.\n");
 
  return r;
}

/* find x in cartesian coord correspond to X (general coords) 
// in the given patch.
// ->return value 1 if it is successful, otherwise 0.
*/
int x_of_X_precision(double *const x,const double *const X,const Patch_T *const patch,const double precision_factor)
{
  int ret = 0;
  
  if (patch->coordsys == Cartesian)
    ret = x_of_X_Cartesian_coord(x,X,patch);
  else if (patch->coordsys == CubedSpherical)
    ret = x_of_X_CS_coord(x,X,patch,precision_factor,1);
  else
    Error0(NO_JOB);
 
  return ret;
}

/* find x in cartesian coord correspond to X (general coords) 
// for Cartesian coord. Note: x reported with respect to the origin (0,0,0)
// ->return value 1 if it is successful, otherwise 0. */
static int x_of_X_Cartesian_coord(double *const x,const double *const X,const Patch_T *const patch)
{
  x[0] = X[0];
  x[1] = X[1];
  x[2] = X[2];
  
  /* test if this is a valid answer */
  if (IsInside(x,patch->min,patch->max,EPS_coord_general))
    return 1;
 
  return 0;
}

/* find x in cartesian coord correspond to X (general coords) 
// for Cubed Spherical. Note: x reported with respect to the origin (0,0,0)
// if check_flg = 1, it checks the solution.
// ->return value 1 if it is successful, otherwise 0. */
static int x_of_X_CS_coord(double *const x,const double *const X,const Patch_T *const patch,const double precision_factor,const int check_flg)
{
  const Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
  const Flag_T type = patch->CoordSysInfo->CubedSphericalCoord->type;
  double S;/* sign */
  Uint a,b,c;/* permuted indices */
  Field_T *R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f,
          *R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
  const double xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1,
               xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
  double R1 = patch->CoordSysInfo->CubedSphericalCoord->R1,
         R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
  const double *const C = patch->c;/* center of origine translated */
  double x1,x2,d,ratio;
  double x_test[3],X_test[3],dX;
  
  SignAndIndex_permutation_CubedSphere(side,&a,&b,&c,&S);

  switch (type)
  {
    case OB_T_SCS:
      d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1   = S*(xc1 == DBL_MAX ? R_interpolation_CS(R1_f,X)/d : xc1);
      x2   = S*(xc2 == DBL_MAX ? R_interpolation_CS(R2_f,X)/d : xc2);
      
      x[c] = x1+(x2-x1)*X[2];
      x[a] = S*x[c]*X[0];
      x[b] = S*x[c]*X[1];
      
      x[a]+= C[a];
      x[b]+= C[b];
      x[c]+= C[c];
      
      break;
    case OT_T_SCS:
      R1 = R_interpolation_CS(R1_f,X);
      R2 = R_interpolation_CS(R2_f,X);
      ratio = 1.-R1/R2;
      d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = S*(xc1 == DBL_MAX ? R1/d : xc1);
      x[c] = x1/(1.-ratio*X[2]);
      x[a] = S*x[c]*X[0];
      x[b] = S*x[c]*X[1];
      
      x[a]+= C[a];
      x[b]+= C[b];
      x[c]+= C[c];
      
    break;
    case NS_T_CS:
      d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = xc1;
      x2 = S*R_interpolation_CS(R2_f,X)/d;
      
      x[c] = x1+(x2-x1)*X[2];
      x[a] = S*x[c]*X[0];
      x[b] = S*x[c]*X[1];
      
      x[a]+= C[a];
      x[b]+= C[b];
      x[c]+= C[c];
    break;
    case SR_T_CS:
      d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x2 = xc2;
      x1 = S*R_interpolation_CS(R1_f,X)/d;
      
      x[c] = x1+(x2-x1)*X[2];
      x[a] = S*x[c]*X[0];
      x[b] = S*x[c]*X[1];
      
      x[a]+= C[a];
      x[b]+= C[b];
      x[c]+= C[c];
    break;
    case OT_T1_CS:
      d  = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = xc1;
      ratio  = 1.-S*d*xc1/R2;
      
      x[c] = x1/(1-ratio*X[2]);
      x[a] = S*x[c]*X[0];
      x[b] = S*x[c]*X[1];
      
      x[a]+= C[a];
      x[b]+= C[b];
      x[c]+= C[c];
    break;
    case OT_T2_CS:
      ratio = 1.-R1/R2;
      d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = S*R1/d;
      
      x[c] = x1/(1.-ratio*X[2]);
      x[a] = S*x[c]*X[0];
      x[b] = S*x[c]*X[1];
      
      x[a]+= C[a];
      x[b]+= C[b];
      x[c]+= C[c];
    break;
    default:
      Error0(NO_OPTION);
  }
  
  /* test the solution */
  if (check_flg)
  {
    x_test[0] = x[0];
    x_test[1] = x[1];
    x_test[2] = x[2];
    X_of_x_CS_coord(X_test,x_test,patch,precision_factor,0);
    dX = root_square(3,X,X_test);
    
    if (!EQL(dX,0))
      return 0;
  }
  
  return 1;
}

/* find point X correspond to x for patch with Cartesian coord.
// ->return value: 1 if it is successful, otherwise 0.
*/
static int X_of_x_Cartesian_coord(double *const X,const double *const x,const Patch_T *const patch)
{
  X[0] = x[0];
  X[1] = x[1];
  X[2] = x[2];
  
  /* test if this is a valid answer */
  if (IsInside(X,patch->min,patch->max,EPS_coord_general))
  {
    /* adjusting boundary values to avoid some unexpeted behavior
    // for example at interpolation. */
    #if 0
    if (EQL_coord(X[0],patch->max[0],EPS_coord_general))
      X[0] = patch->max[0];
    if (EQL_coord(X[0],patch->min[0],EPS_coord_general))
      X[0] = patch->min[0];
    if (EQL_coord(X[1],patch->max[1],EPS_coord_general))
      X[1] = patch->max[1];
    if (EQL_coord(X[1],patch->min[1],EPS_coord_general))
      X[1] = patch->min[1];
    if (EQL_coord(X[2],patch->max[2],EPS_coord_general))
      X[2] = patch->max[2];
    if (EQL_coord(X[2],patch->min[2],EPS_coord_general))
      X[2] = patch->min[2];
    #endif
      
    /* adjust X */
    X[0] = (EQL(X[0],patch->max[0]) ? patch->max[0] : X[0]);
    X[0] = (EQL(X[0],patch->min[0]) ? patch->min[0] : X[0]);

    X[1] = (EQL(X[1],patch->max[1]) ? patch->max[1] : X[1]);
    X[1] = (EQL(X[1],patch->min[1]) ? patch->min[1] : X[1]);

    X[2] = (EQL(X[2],patch->max[2]) ? patch->max[2] : X[2]);
    X[2] = (EQL(X[2],patch->min[2]) ? patch->min[2] : X[2]);
    
    return 1;
  }
  
  return 0;
}

/* find point X correspond to cart-coord for patch with cubed spherical coord.
// it's a general algorithm and for even if the point is not collocated.
// if check_flg = 1, it checks the solution.
// ->return value: 1 if it is successful, otherwise 0. */
static int X_of_x_CS_coord(double *const X,
                           const double *const cart,
                           const Patch_T *const patch,
                           const double precision_factor,
                           const int check_flg)
{
  const double *const C = patch->c;/* center of origine translated */
  const Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
  const Flag_T type = patch->CoordSysInfo->CubedSphericalCoord->type;
  const double x[3]= {cart[0]-C[0],
                      cart[1]-C[1],
                      cart[2]-C[2]};
  double S; /* sign */
  Uint i,j,k;/* permuted indices */
  Field_T *R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f,
          *R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
  const double xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1,
               xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
  double R1 = patch->CoordSysInfo->CubedSphericalCoord->R1,
         R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
  double x1,x2,d,ratio;
  double x_test[3],X_test[3],dx;
  double eps = EPS_coord_general*precision_factor;
  
  SignAndIndex_permutation_CubedSphere(side,&i,&j,&k,&S);
  
  X[0] = S*x[i]/x[k];
  X[1] = S*x[j]/x[k];
  
  switch (type)
  {
    case OB_T_SCS:
      d    = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1   = S*(xc1 == DBL_MAX ? R_interpolation_CS(R1_f,X)/d : xc1);
      x2   = S*(xc2 == DBL_MAX ? R_interpolation_CS(R2_f,X)/d : xc2);
      X[2] = (x[k]-x1)/(x2-x1);
      
      #if 0
      const Uint *n;
      /*  for interpolation error */
      n = patch->n;
      if (patch->nsplit[2] == 1)
      {
        if (n[2] < LOW_n)
          eps = EPS_coord_LOW_n1/(n[0]*n[1]*n[2]);
        else
          eps = EPS_coord_OB_SCS1/(n[0]*n[1]*n[2]);
      }
      else
      {
        if (n[2] < LOW_n)
          eps = EPS_coord_LOW_n2/(n[0]*n[1]*n[2]);
        else
          eps = EPS_coord_OB_SCS2/(n[0]*n[1]*n[2]);
      }
      #endif
      
    break;
    case OT_T_SCS:
      R1 = R_interpolation_CS(R1_f,X);
      R2 = R_interpolation_CS(R2_f,X);
      ratio = 1.-R1/R2;
      d  = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = S*(xc1 == DBL_MAX ? R1/d : xc1);
      X[2] = (1-x1/x[k])/ratio;

      #if 0
      const Uint *n;
      /*  for interpolation error */
      n = patch->n;
      if (patch->nsplit[2] == 1)
      {
        if (n[2] < LOW_n)
          eps = EPS_coord_LOW_n1/(n[0]*n[1]*n[2]);
        else
          eps = EPS_coord_OT_SCS1/(n[0]*n[1]*n[2]);
      }
      else
      {
        if (n[2] < LOW_n)
          eps = EPS_coord_LOW_n2/(n[0]*n[1]*n[2]);
        else
          eps = EPS_coord_OT_SCS2/(n[0]*n[1]*n[2]);
      }
      #endif
      
    break;
    case NS_T_CS:
      d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = xc1;
      x2 = S*R_interpolation_CS(R2_f,X)/d;
      X[2] = (x[k]-x1)/(x2-x1);
    break;
    case SR_T_CS:
      d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x2 = xc2;
      x1 = S*R_interpolation_CS(R1_f,X)/d;
      X[2] = (x[k]-x1)/(x2-x1);
    break;
    case OT_T1_CS:
      d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      ratio = 1.-S*d*xc1/R2;
      X[2] = (1-xc1/x[k])/ratio;
    break;
    case OT_T2_CS:
      ratio = 1.-R1/R2;
      d  = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = S*R1/d;
      X[2] = (1-x1/x[k])/ratio;
    break;
    default:
      Error0(NO_OPTION);
  }
  
  /* adjust X */
  X[0] = (EQL(X[0],patch->max[0]) ? patch->max[0] : X[0]);
  X[0] = (EQL(X[0],patch->min[0]) ? patch->min[0] : X[0]);

  X[1] = (EQL(X[1],patch->max[1]) ? patch->max[1] : X[1]);
  X[1] = (EQL(X[1],patch->min[1]) ? patch->min[1] : X[1]);

  X[2] = (EQL(X[2],patch->max[2]) ? patch->max[2] : X[2]);
  X[2] = (EQL(X[2],patch->min[2]) ? patch->min[2] : X[2]);  
  
  /* test the solution, don't modify X. */
  if (check_flg)
  {
    if (!IsInside(X,patch->min,patch->max,eps))
      return 0;
      
    /* test if it gives you the same x coords */
    X_test[0] = X[0];
    X_test[1] = X[1];
    X_test[2] = X[2];
    
    x_of_X_CS_coord(x_test,X_test,patch,precision_factor,0);
    dx = root_square(3,cart,x_test);
    double scale = MaxMag_d(root_square(3,cart,0),root_square(3,x_test,0));
    scale = scale < 1 ? 1 : scale;
    if (!EQL(dx/scale,0.))
      return 0;
  }
  
  /* extra care for patches with difficult X[2] like outermost patches. */
  X[2] = (X[2] > patch->max[2] ? patch->max[2] : X[2]);
  X[2] = (X[2] < patch->min[2] ? patch->min[2] : X[2]);
  
  return 1;
}

/* given point and patch find if the is any node collocated 
// to that point and then return its index.
// ->return value: found index, and put flg = FOUND, otherwise, flg = NONE.
*/
Uint find_node(const double *const x, const Patch_T *const patch,Flag_T *const flg)
{
  Uint v = UINT_MAX;
  double res = EPS_collocation*root_square(3,x,0);/* resolution */
  Uint i;
  double *y, nrm;

  res = GRT(res,EPS_collocation) ? res : EPS_collocation;
  *flg = NONE;
    
  FOR_ALL(i,patch->node)
  {
    y = patch->node[i]->x;
    nrm = root_square(3,x,y);
    if (LSSEQL_coord(nrm,res,EPS_coord_general))
    {
      v = i;
      res = nrm;
      *flg = FOUND;
    }
  }

  return v;
}

/* ->return value: x coord in specified patch */
double x_coord(const Uint i,const Patch_T *const patch)
{
  return patch->node[i]->x[0];
}

/* ->return value: y coord in specified patch */
double y_coord(const Uint i,const Patch_T *const patch)
{
  return patch->node[i]->x[1];
}

/* ->return value: z coord in specified patch */
double z_coord(const Uint i,const Patch_T *const patch)
{
  return patch->node[i]->x[2];
}

/* ->return value: X coord in specified patch */
double X_coord(const Uint i,const Patch_T *const patch)
{
  return patch->node[i]->X[0];
}

/* ->return value: Y coord in specified patch */
double Y_coord(const Uint i,const Patch_T *const patch)
{
  return patch->node[i]->X[1];
}

/* ->return value: Z coord in specified patch */
double Z_coord(const Uint i,const Patch_T *const patch)
{
  return patch->node[i]->X[2];
}

/* make an empty needle 
// ->return value: pointer to new needle
*/
void *alloc_needle(void)
{
  Needle_T *needle;
  
  needle = calloc(1,sizeof(*needle));
  IsNull(needle);

  needle->precision_factor = 1.;
  return needle;
}

/* free needle */
void free_needle(Needle_T *needle)
{
  if (needle == 0) return;
  
  else
  {
    if (needle->Nin  != 0) free (needle->in);
    if (needle->Nex  != 0) free (needle->ex);
    if (needle->Ng   != 0) free (needle->guess);
    if (needle->Nans != 0) free (needle->ans);
    
  }
  
  free(needle);
}

/* given (X,Y) in a specified surface of Z from NS/BH in cubed spherical coords
// it finds the associated polar and azimuthal angels on the surface.
// note: it does not depend on Z */
void theta_phi_of_XY_CS(double *const theta,double *const phi,const double *const X,const Flag_T side)
{
  const double a = X[0];
  const double b = X[1];
  const double d = sqrt(1+Pow2(a)+Pow2(b));
  
  switch (side)
  {
    case UP:
      *phi   = arctan(b,a);
      *theta = acos(1/d);
    break;
    case DOWN:
      *phi   = arctan(a,b);
      *theta = acos(-1/d);
    break;
    case LEFT:
      *phi   = arctan(-1,a);
      *theta = acos(b/d);
    break;
    case RIGHT:
      *phi   = arctan(1,b);
      *theta = acos(a/d);
    break;
    case BACK:
      *phi   = arctan(b,-1);
      *theta = acos(a/d);
    break;
    case FRONT:
      *phi   = arctan(a,1);
      *theta = acos(b/d);
    break;
    default:
      Error0(NO_OPTION);
  }
  
  /* more test */
  if(0)
  {
    double th = *theta;
    double ph = *phi;
    switch (side)
    {
      case UP:
        if (!EQL(X[0],tan(th)*cos(ph)))
          Error0("Wrong transformation");
        if (!EQL(X[1],tan(th)*sin(ph)))
          Error0("Wrong transformation");
      break;
      case DOWN:
        if (!EQL(X[0],-tan(th)*sin(ph)))
          Error0("Wrong transformation");
        if (!EQL(X[1],-tan(th)*cos(ph)))
          Error0("Wrong transformation");
      break;
      case LEFT:
        if (!EQL(X[0],-1./tan(ph)))
          Error0("Wrong transformation");
        if (!EQL(X[1],-1./(tan(th)*sin(ph))))
          Error0("Wrong transformation");
      break;
      case RIGHT:
        if (!EQL(X[0],1./(tan(th)*sin(ph))))
          Error0("Wrong transformation");
        if (!EQL(X[1],1./tan(ph)))
          Error0("Wrong transformation");
      break;
      case BACK:
        if (!EQL(X[0],-1./(tan(th)*cos(ph))))
          Error0("Wrong transformation");
        if (!EQL(X[1],-tan(ph)))
          Error0("Wrong transformation");
      break;
      case FRONT:
        if (!EQL(X[0],tan(ph)))
          Error0("Wrong transformation");
        if (!EQL(X[1],1/(tan(th)*cos(ph))))
          Error0("Wrong transformation");
      break;
      default:
        Error0(NO_OPTION);
    }
  }
}

/* ->: collected patches which cover the region in 
// patch->CoordSysInfo->region and number of patches Np.
// one can collect an assortment of patches separated with comma, eg:
//
// Patch_T **patches = collect_patches(grid,"NS1_around,NS1",&np); 
// which cover NS1 and NS1_around regardless of direction. 
// note: it's blind with respect to repetition in the specified region.
// NOTE: if region = ".*" it collects all of the grid patches. */
Patch_T **
collect_patches
  (
  Grid_T *const grid,/* the grid */
  const char *const region,/* see the list in IsItCovering function */
  Uint *const Np/* number of patches found */
  )
{
  Patch_T **patches = 0;
  const int IsMatchAll = !strcmp(region,".*");
  Uint np,p;
  
  /* init */
  *Np = 0;
  
  np = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if (IsMatchAll || IsItCovering(patch,region))
    {
      patches = realloc(patches,(np+2)*sizeof(*patches));
      IsNull(patches);
      patches[np]   = patch;
      patches[np+1] = 0;
      ++np;
    }
  }
  
  /* check if there is no such region */
  if (np == 0)
    Errors("No such '%s'!",region);
  
  *Np = np;
  return patches;
}

/* ->: Is this patch covering this region? yes = 1, no = 0. 
// list of regions (? can be 1 2 or nothing):
// 
// "NS?" == the whole NS patches including central box (if any)
// "BH?" == the whole BH patches including central box (if any)
// "NS?_around" == the whole NS around patches
// "BH?_around" == the whole NS around patches
// "outermost" == the whole outermost patches
// "filling_box" == only patches cover the filling box
// "central_box" == only patches cover the central box
// NOTE: one can add suffix _OB or IB if interested in outer boundary and
//       inner boundary of patches (the surfaces of the patches)
//       some examples below:
//
// "NS?_OB" == only patches include the Outer Boundary i.e NS surface from inside.
// "BH?_OB" == only patches include the Outer Boundary i.e BH surface from inside.
// "NS?_around_IB" == only patches include the Inner Boundary i.e NS surface from outside
// "BH?_around_IB" == only patches include the Inner Boundary i.e BH surface from outside
// "NS?_around_OB" == only patches include the Outer Boundary i.e farthest surface
// "BH?_around_OB" == only patches include the Outer Boundary i.e farthest surface
//
// ex:
// ===
// IsItCovering(patch,"outermost");     => outemost patch?
// IsItCovering(patch,"NS_OB");         => NS_surface patch?
// IsItCovering(patch,"NS2_OB");        => NS2_surface patch?
// IsItCovering(patch,"NS1,NS2_OB");    => NS1 or NS2_surface patch?
// IsItCovering(patch,"NS2_around_IB"); => patches cover NS2 surface from around patches
*/
int 
IsItCovering
  (
  const Patch_T *const patch,/* the patch */
  const char *const region/* BH/NS etc. see the list above */
  )
{
  Grid_T *const grid = patch->grid;
  char **reg = read_separated_items_in_string(region,',');
  char s[999] = {'\0'};
  Uint i = 0;
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_BHBH ||
      grid->kind == Grid_SplitCubedSpherical_SNS  ||
      grid->kind == Grid_SplitCubedSpherical_SBH
     )
  {
    i = 0;
    while (reg[i])
    {
      /* if this is a patch name request. */
      if (strstr(reg[i],PATCH_NAME_P_))
      {
        /* does match? */
        if (!strcmp(patch->name,reg[i]))
        {
          free_2d(reg);
          return 1;
        } 
      }
      else
      {
        sprintf(s,"(%s)",reg[i]);
        /* check the spell of region (for debug purposes) */
        if (SPELL_CHECK)
        {
          if (!strstr(Dictionary_SCS,s))
            Errors("Please spell check '%s'!",reg[i]);
        }
        
        /* does cover? */
        if (strstr(patch->CoordSysInfo->region,s))
        {
          free_2d(reg);
          return 1;
        }
      }
      i++;
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
  free_2d(reg);
  
  return 0;
}

/* ->: initialize a grid character. */
Grid_Char_T *init_grid_char(Grid_T *const new_grid)
{
  Grid_Char_T *g = calloc(1,sizeof(*g));IsNull(g);
  g->grid = new_grid;
  
  return g;
}

/* free all entries of Grid_Char_T */
void free_grid_char(Grid_Char_T *g)
{
  Uint i;
  
  if (!g)
    return;
  
  for (i = 0; i < NPARAMS_GRID_CHAR; ++i)
  {
    if (g->params[i]->occupied)
    {
    Free(g->params[i]->relClm);
    Free(g->params[i]->imgClm);
    }
  }
  Free(g);
}

Grid_Kind_T set_grid_kind(const char *const grid_kind)
{
  Grid_Kind_T ret = Grid_UNDEFINED;
  
  if (strcmp_i(grid_kind,"SplitCubedSpherical(BH+NS)"))
    ret = Grid_SplitCubedSpherical_BHNS;
    
  else if (strcmp_i(grid_kind,"SplitCubedSpherical(NS+NS)"))
    ret = Grid_SplitCubedSpherical_NSNS;
    
  else if (strcmp_i(grid_kind,"SplitCubedSpherical(BH+BH)"))
    ret = Grid_SplitCubedSpherical_BHBH;
    
  else if (strcmp_i(grid_kind,"SplitCubedSpherical(NS)"))
    ret = Grid_SplitCubedSpherical_SNS;
    
  else if (strcmp_i(grid_kind,"SplitCubedSpherical(BH)"))
    ret = Grid_SplitCubedSpherical_SBH;
    
  else if (strcmp_i(grid_kind,"CubedSpherical(BH+NS)"))
   ret = Grid_CubedSpherical_BHNS;
   
  else if (strcmp_i(grid_kind,"CubedSpherical(NS+NS)"))
   ret = Grid_CubedSpherical_NSNS;
   
  else if (strcmp_i(grid_kind,"box"))
   ret = Grid_Box; 
   
  else
    Errors("There is no such %s grid kind.\n",grid_kind);
  
  return ret;
}

/* given (X,Y,Z) in cubed spherical or split cubed spherical coords
// it finds the associated polar and azimuthal angels for given X.
// works for both split and normal cubed spherical. */
void find_theta_phi_of_XYZ_CS(double *const theta,double *const phi,
                              const double *const X,const Flag_T side)
{
  const double a = X[0];
  const double b = X[1];
  const double d = sqrt(1+Pow2(a)+Pow2(b));
  
  switch (side)
  {
    case UP:
      *phi   = arctan(b,a);
      *theta = acos(1/d);
    break;
    case DOWN:
      *phi   = arctan(a,b);
      *theta = acos(-1/d);
    break;
    case LEFT:
      *phi   = arctan(-1,a);
      *theta = acos(b/d);
    break;
    case RIGHT:
      *phi   = arctan(1,b);
      *theta = acos(a/d);
    break;
    case BACK:
      *phi   = arctan(b,-1);
      *theta = acos(a/d);
    break;
    case FRONT:
      *phi   = arctan(a,1);
      *theta = acos(b/d);
    break;
    default:
      Error0(NO_OPTION);
  }
  
  /* more test */
  if(DEBUGGING)
  {
    double th = *theta;
    double ph = *phi;
    switch (side)
    {
      case UP:
        if (!EQL(X[0],tan(th)*cos(ph)))
          Error0("Wrong transformation");
        if (!EQL(X[1],tan(th)*sin(ph)))
          Error0("Wrong transformation");
      break;
      case DOWN:
        if (!EQL(X[0],-tan(th)*sin(ph)))
          Error0("Wrong transformation");
        if (!EQL(X[1],-tan(th)*cos(ph)))
          Error0("Wrong transformation");
      break;
      case LEFT:
        if (!EQL(X[0],-1./tan(ph)))
          Error0("Wrong transformation");
        if (!EQL(X[1],-1./(tan(th)*sin(ph))))
          Error0("Wrong transformation");
      break;
      case RIGHT:
        if (!EQL(X[0],1./(tan(th)*sin(ph))))
          Error0("Wrong transformation");
        if (!EQL(X[1],1./tan(ph)))
          Error0("Wrong transformation");
      break;
      case BACK:
        if (!EQL(X[0],-1./(tan(th)*cos(ph))))
          Error0("Wrong transformation");
        if (!EQL(X[1],-tan(ph)))
          Error0("Wrong transformation");
      break;
      case FRONT:
        if (!EQL(X[0],tan(ph)))
          Error0("Wrong transformation");
        if (!EQL(X[1],1/(tan(th)*cos(ph))))
          Error0("Wrong transformation");
      break;
      default:
        Error0(NO_OPTION);
    }
  }
}


/* given theta, phi and knowing the fact that they are on a surface, 
// it finds the corresponding patch and X,Y,Z coordinate. 
// works for both split and normal cubed spherical. */
void 
find_XYZ_and_patch_of_theta_phi_CS
 (
 double *const X/* found X,Y,Z Note: X[2] must be filled 
                // to determine the surface */,
 Patch_T **const ppatch,/* found patch */
 const double *const center/* center of S2 in general is not patch->c */,
 const double theta/* given theta */,
 const double phi/* given phi */,
 Patch_T **const patches,/* search among these patches */
 const Uint Np/* number of patches */
 )
{
  const double tan_phi    = tan(phi);
  const double cos_theta  = cos(theta);
  const double tan_phi2   = Pow2(tan_phi);
  const double cos_theta2 = Pow2(cos_theta);
  Flag_T found_flg = NO;
  Uint p;
  
  /* check all of the given patches in which (x,y,z) and 
  // (X,Y,Z) and (theta,phi) are consistent */
  for (p = 0; p < Np; ++p)
  {
    Patch_T *patch = patches[p];
    
    Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
    double a = 0, b = 0;
    double a_sign = 0,b_sign = 0,c_sign = 0;
    double x[3],phi2,theta2,r;
    int is_inside;
    
    /* first calculate the magnetitude of a and b 
    // which are related to X[0] and X[1] with a sign */
    switch (side)
    {
      case UP:
        a = sqrt((1 - cos_theta2)/(cos_theta2 + cos_theta2*tan_phi2));
        b = tan_phi*sqrt((1 - cos_theta2)/(cos_theta2*(1 + tan_phi2)));
      break;
      case DOWN:
        b = sqrt((1 - cos_theta2)/(cos_theta2 + cos_theta2*tan_phi2));
        a = tan_phi*sqrt((1 - cos_theta2)/(cos_theta2*(1 + tan_phi2)));
      break;
      case LEFT:
        a = 1/tan_phi;
        b = sqrt((cos_theta2 + cos_theta2*tan_phi2)/((1. - cos_theta2)*tan_phi2));
      break;
      case RIGHT:
        b = 1/tan_phi;
        a = sqrt((cos_theta2 + cos_theta2*tan_phi2)/((1. - cos_theta2)*tan_phi2));
      break;
      case BACK:
        a = sqrt((cos_theta2 + cos_theta2*tan_phi2)/(1 - cos_theta2));
        b = tan_phi;
      break;
      case FRONT:
        b = sqrt((cos_theta2 + cos_theta2*tan_phi2)/(1 - cos_theta2));
        a = tan_phi;
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* having found the magnitude of a and b, we need to find out the sign of them.
    // this is done by paying attention to side, signum(cos_theta) and range of tanphi */
    switch (side)
    {
      case UP:
        arctan_argument_signum(&b_sign,&a_sign,phi);
      break;
      case DOWN:
        arctan_argument_signum(&a_sign,&b_sign,phi);
      break;
      case LEFT:
        arctan_argument_signum(&c_sign,&a_sign,phi);
        if (cos_theta > 0) b_sign = 1;
        else		   b_sign = -1;
      break;
      case RIGHT:
        arctan_argument_signum(&c_sign,&b_sign,phi);
        if (cos_theta > 0) a_sign = 1;
        else		   a_sign = -1;
      break;
      case BACK:
        arctan_argument_signum(&b_sign,&c_sign,phi);
        if (cos_theta > 0) a_sign = 1;
        else		   a_sign = -1;
      break;
      case FRONT:
        arctan_argument_signum(&a_sign,&c_sign,phi);
        if (cos_theta > 0) b_sign = 1;
        else		   b_sign = -1;
      break;
      default:
        Error0(NO_OPTION);
    }
    
    X[0] = fabs(a)*a_sign;
    X[1] = fabs(b)*b_sign;
    
    /* check if x of X really gives you the correct angles */
    is_inside = x_of_X(x,X,patch);
    x[0] -= center[0];
    x[1] -= center[1];
    x[2] -= center[2];
    r = root_square(3,x,0);
    theta2 = acos(x[2]/r);
    phi2   = arctan(x[1],x[0]);
    if (EQL(theta2,theta) && EQL(phi2,phi) && is_inside)
    {
      found_flg = YES;
      *ppatch = patch;
      break;
    }
  }
  if (found_flg == NO)
    Error0("(X,Y,Z) or patch could not be found.\n");
}

/* ->: first patch has the Cartesian point x or null if no patch has this.
// get a Cartesian point x and collection of patches,
// it returns the first patch has this point.
// Np is the number of patches. */
Patch_T *x_in_which_patch(const double x[3],Patch_T **const patches,
                          const Uint Np)
{
  double X[3] = {0};
  Uint p;
  
  for (p = 0; p < Np; ++p)
  {
    Patch_T *patch = patches[p];
    
    if (X_of_x(X,x,patch))
      return patch;
  }
  
  return 0;
}

/* ->: find forcefully the patch that has the Cartesian point x;
// ->: null if failes. also finds the corresponding X coordinates.
// get a Cartesian point x and collection of patches,
// it returns the closest patch that has this point.
// Np is the number of patches. 
// Note: sometimes x_in_which_patch may fail and we have to use
// this function; so ONLY use this if x_in_which_patch failes. it,
// by and large, used for outermost patches when there is a big jump 
// in resolution like 10->16. */
Patch_T *x_in_which_patch_force(const double x[3],Patch_T **const patches,
                                const Uint Np,double *const X)
{
  if (!Np)
    return 0;
  
  const double Eps   = EPS_coord_general*10.;/* from experiment */
  Patch_T *patch     = 0;
  double min = DBL_MAX;
  Uint pmin  = UINT_MAX;/* patch with min dx */
  Uint p;
  
  for (p = 0; p < Np; ++p)
  {
    patch        = patches[p];
    double Xp[3] = {0};
    double xp[3] = {0};
    double dx;
    X_of_x(Xp,x,patch);
    
    /* only if it is inside count it */
    if (IsInside(Xp,patch->min,patch->max,Eps))
    {
      /* adjust */
      Xp[0] = (Xp[0] > patch->max[0] ? patch->max[0] : Xp[0]);
      Xp[0] = (Xp[0] < patch->min[0] ? patch->min[0] : Xp[0]);
      Xp[1] = (Xp[1] > patch->max[1] ? patch->max[1] : Xp[1]);
      Xp[1] = (Xp[1] < patch->min[1] ? patch->min[1] : Xp[1]);
      Xp[2] = (Xp[2] > patch->max[2] ? patch->max[2] : Xp[2]);
      Xp[2] = (Xp[2] < patch->min[2] ? patch->min[2] : Xp[2]);
      
      x_of_X(xp,Xp,patch);
      dx = L2_norm(3,x,xp);
      if (dx < min)
      {
        pmin = p;
        min  = dx;
      }
    }
  }
  
  if (pmin == UINT_MAX)
    return 0;
  
  patch = patches[pmin];
  X_of_x(X,x,patch);
  /* since we gonna pick this patch anyway let's adjust boundaries 
  // to avoid some unexpeted behavior. for example at interpolation. */
  X[0] = (X[0] > patch->max[0] ? patch->max[0] : X[0]);
  X[0] = (X[0] < patch->min[0] ? patch->min[0] : X[0]);
  X[1] = (X[1] > patch->max[1] ? patch->max[1] : X[1]);
  X[1] = (X[1] < patch->min[1] ? patch->min[1] : X[1]);
  X[2] = (X[2] > patch->max[2] ? patch->max[2] : X[2]);
  X[2] = (X[2] < patch->min[2] ? patch->min[2] : X[2]);  
  
  /* test it */
  if(0)
  {
    double test_x[3] = {0};
    double test_dx;
    assert(x_of_X(test_x,X,patch));
    test_dx = L2_norm(3,test_x,x);
    printf(Pretty0"'%s'|displacement = %e\n",patch->name,test_dx);
  }
  
  return patch;
}



/* ->: first patch has patch coords. X or null if no patch has this.
// get a patch coords. X and collection of patches,
// it returns the first patch has this point.
// Np is the number of patches. 
// some times, you know the X coords but not the patch, for instance, 
// in mass_shedding function this becomes handy.
// note: the collection of patch must be refine enough to
// avoid potential wrong return since X is not uniqe. */
Patch_T *X_in_which_patch(const double X[3],Patch_T **const patches,
                          const Uint Np)
{
  double x[3] = {0};
  Uint p;
  
  for (p = 0; p < Np; ++p)
  {
    Patch_T *patch = patches[p];
    
    if (x_of_X(x,X,patch))
      return patch;
  }
  
  return 0;
}

/* ->: collected patches which have the given regex match in patch->name.
// ex:
// Patch_T **patches = collect_patches_regex(grid,".+(left|right)_NS.?_around.+",&np); 
// which in a NSNS system covers NS1_around and NS2_around regardless of direction. */
Patch_T **
collect_patches_regex
  (
  Grid_T *const grid,/* the grid */
  const char *const regex,/* regex */
  Uint *const Np/* number of patches found */
  )
{
  Patch_T **patches = 0;
  Uint np,p;
  
  /* init */
  *Np = 0;
  
  np = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if (regex_search(regex,patch->name))
    {
      patches = realloc(patches,(np+2)*sizeof(*patches));
      IsNull(patches);
      patches[np]   = patch;
      patches[np+1] = 0;
      ++np;
    }
  }
  /* check */
  if (np == 0)
    Errors("No name was matched for '%s'!",regex);
  
  *Np = np;
  return patches;
}

