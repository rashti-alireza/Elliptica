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
  unsigned i,j;
  
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
  unsigned i;
  
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
  unsigned i;
  
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
  unsigned i;
  
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
  unsigned i;
  
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
  unsigned *p = 0, np = UINT_MAX;
  unsigned i;
  
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
    
    a = X_of_x(X,needle->x,patch);
    
    if (a)
      needle_ans(needle,patch);
      
  }
}

/* find point X(general coords) correspond to x(Cartesian coords) 
// in the given patch.
// ->return value 1 if it is successful, otherwise 0.
*/
int X_of_x(double *const X,const double *const x,const Patch_T *const patch)
{
  int r = 0;
  
  if (patch->coordsys == Cartesian)
    r = X_of_x_Cartesian_coord(X,x,patch);
  else if (patch->coordsys == CubedSpherical)
    r = X_of_x_CS_coord(X,x,patch,1);
  else
      Error0("No finder for this coordinate.\n");
 
  return r;
}

/* find x in cartesian coord correspond to X (general coords) 
// in the given patch.
// ->return value 1 if it is successful, otherwise 0.
*/
int x_of_X(double *const x,const double *const X,const Patch_T *const patch)
{
  int ret = 0;
  
  if (patch->coordsys == Cartesian)
    ret = x_of_X_Cartesian_coord(x,X,patch);
  else if (patch->coordsys == CubedSpherical)
    ret = x_of_X_CS_coord(x,X,patch,1);
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
static int x_of_X_CS_coord(double *const x,const double *const X,const Patch_T *const patch,const int check_flg)
{
  const Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
  const Flag_T type = patch->CoordSysInfo->CubedSphericalCoord->type;
  double S;/* sign */
  unsigned a,b,c;/* permuted indices */
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
    case OJ_T_SCS:
      d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = S*R_interpolation_CS(R1_f,X)/d;
      x2 = S*R_interpolation_CS(R2_f,X)/d;
      
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
      x1 = S*R1/d;
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
    X_of_x_CS_coord(X_test,x_test,patch,0);
    dX = root_square(3,X,X_test);
    
    if (!EQL_coord(dX,0,EPS_coord_general))
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
    return 1;
  
  return 0;
}

/* find point X correspond to cart-coord for patch with cubed spherical coord.
// it's a general algorithm and for even if the point is not collocated.
// if check_flg = 1, it checks the solution.
// ->return value: 1 if it is successful, otherwise 0. */
static int X_of_x_CS_coord(double *const X,const double *const cart,const Patch_T *const patch,const int check_flg)
{
  const double *const C = patch->c;/* center of origine translated */
  const Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
  const Flag_T type = patch->CoordSysInfo->CubedSphericalCoord->type;
  const double x[3]= {cart[0]-C[0],
                      cart[1]-C[1],
                      cart[2]-C[2]};
  double S; /* sign */
  unsigned i,j,k;/* permuted indices */
  Field_T *R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f,
          *R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
  const double xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1,
               xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
  double R1 = patch->CoordSysInfo->CubedSphericalCoord->R1,
         R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
  double x1,x2,d,ratio;
  double x_test[3],X_test[3],dx;
  double eps = EPS_coord_general;
  const unsigned *n;
  
  SignAndIndex_permutation_CubedSphere(side,&i,&j,&k,&S);
  
  X[0] = S*x[i]/x[k];
  X[1] = S*x[j]/x[k];
  
  switch (type)
  {
    case OJ_T_SCS:
      d    = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1   = S*R_interpolation_CS(R1_f,X)/d;
      x2   = S*R_interpolation_CS(R2_f,X)/d;
      X[2] = (x[k]-x1)/(x2-x1);
      
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
      
    break;
    case OT_T_SCS:
      R1 = R_interpolation_CS(R1_f,X);
      R2 = R_interpolation_CS(R2_f,X);
      ratio = 1.-R1/R2;
      d  = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = S*R1/d;
      X[2] = (1-x1/x[k])/ratio;
      
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
  
  /* adujusting boundary number to avoid some unexpeted behavior
  // due to interpolation error. */
  if (EQL_coord(X[0],patch->max[0],eps))  X[0] = patch->max[0];
  if (EQL_coord(X[0],patch->min[0],eps))  X[0] = patch->min[0];
  if (EQL_coord(X[1],patch->max[1],eps))  X[1] = patch->max[1];
  if (EQL_coord(X[1],patch->min[1],eps))  X[1] = patch->min[1];
  if (EQL_coord(X[2],patch->max[2],eps))  X[2] = patch->max[2];
  if (EQL_coord(X[2],patch->min[2],eps))  X[2] = patch->min[2];  
  
  /* test the solution */
  if (check_flg)
  {
    unsigned interval_test = 0;
    
    if (IsInside(X,patch->min,patch->max,eps))
       interval_test = 1;
    
    if (!interval_test)
      return 0;
      
    /* test if it gives you the same x coords */
    X_test[0] = X[0];
    X_test[1] = X[1];
    X_test[2] = X[2];
    x_of_X_CS_coord(x_test,X_test,patch,0);
    dx = root_square(3,cart,x_test);
    double scale = MaxMag_d(root_square(3,cart,0),root_square(3,x_test,0));
    scale = scale < 1 ? 1 : scale;
    if (!EQL_coord(dx/scale,0,eps))
      return 0;
  }
  
  return 1;
}

/* given point and patch find if the is any node collocated 
// to that point and then return its index.
// ->return value: found index, and put flg = FOUND, otherwise, flg = NONE.
*/
unsigned find_node(const double *const x, const Patch_T *const patch,Flag_T *const flg)
{
  unsigned v = UINT_MAX;
  double res = EPS_collocation*root_square(3,x,0);/* resolution */
  unsigned i;
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
double x_coord(const unsigned i,const Patch_T *const patch)
{
  return patch->node[i]->x[0];
}

/* ->return value: y coord in specified patch */
double y_coord(const unsigned i,const Patch_T *const patch)
{
  return patch->node[i]->x[1];
}

/* ->return value: z coord in specified patch */
double z_coord(const unsigned i,const Patch_T *const patch)
{
  return patch->node[i]->x[2];
}

/* ->return value: X coord in specified patch */
double X_coord(const unsigned i,const Patch_T *const patch)
{
  return patch->node[i]->X[0];
}

/* ->return value: Y coord in specified patch */
double Y_coord(const unsigned i,const Patch_T *const patch)
{
  return patch->node[i]->X[1];
}

/* ->return value: Z coord in specified patch */
double Z_coord(const unsigned i,const Patch_T *const patch)
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

/* ->: collected patches which cover the object and number of patches Np
// which cover the specified region. */
Patch_T **
collect_patches
  (
  Grid_T *const grid,/* the grid */
  const char *const region,/* see the list in IsItCovering function */
  const Flag_T side,/* LEFT or RIGHT or CENTER or NONE */
  unsigned *const Np/* number of patches found */
  )
{
  Patch_T **patches = 0;
  unsigned np,p;
  
  np = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if (IsItCovering(patch,region,side))
    {
      patches = realloc(patches,(np+1)*sizeof(*patches));
      IsNull(patches);
      patches[np] = patch;
      ++np;
    }
  }
  
  /* check if there is no such region */
  if (np == 0)
    Error0("No such region!");
  
  *Np = np;
  return patches;
}

/* ->: Is this patch covering this region? yes = 1, no = 0. 
// list of regions (? can be 1 2 or nothing):
// 
// "NS?" == the whole NS patches including central box (if any)
// "BH?" == the whole BH patches including central box (if any)
// "NS?_surface" == only patches include the NS surface from inside
// "BH?_surface" == only patches include the BH surface from inside
// "NS?_surrounding_surface" == only patches include the NS surface from outside
// "BH?_surrounding_surface" == only patches include the BH surface from outside
// "NS?_surrounding" == the whole NS surrounding patches
// "BH?_surrounding" == the whole NS surrounding patches
// "outermost" == the whole outermost patches
// "filling_box" == only patches cover the filling box
// "central_box" == only patches cover the central box
//
// ex:
// ===
// IsItCovering(patch,"outermost",NONE);   => outemost patch?
// IsItCovering(patch,"NS_surface",LEFT);  => NS_surface patch?
// IsItCovering(patch,"NS2_surface",LEFT); => NS2_surface patch?
*/
int 
IsItCovering
  (
  const Patch_T *const patch,/* the patch */
  const char *const region,/* BH/NS etc. see the list above */
  const Flag_T Fside/* LEFT or RIGHT or CENTER or NONE (side of region, if any) */
  )
{
  Grid_T *const grid = patch->grid;
  int ret = 0;  
  const char *side = 0;
  char s[999] = {'\0'};
   
  if(Fside == LEFT || Fside == RIGHT || Fside == CENTER)
    side = StrSide[Fside];
  else
    side = 0;
    
  if (strcmp_i(grid->kind,"SplitCubedSpherical(BH+NS)") ||
      strcmp_i(grid->kind,"SplitCubedSpherical(NS+NS)") ||
      strcmp_i(grid->kind,"SplitCubedSpherical(BH+BH)") ||
      strcmp_i(grid->kind,"SplitCubedSpherical(NS)")    ||
      strcmp_i(grid->kind,"SplitCubedSpherical(BH)")
     )
  {
    if (side)
      sprintf(s,"(%s_%s)",side,region);
    else
      sprintf(s,"(%s)",region);
    
    /* if the request is obvious */  
    if (strstr_i(patch->CoordSysInfo->region,s))
      return 1;
  }
  else
  {
    Error0(NO_OPTION);
  }
  
  return ret;
}


