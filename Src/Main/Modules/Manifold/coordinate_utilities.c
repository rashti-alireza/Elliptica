/*
// Alireza Rashti
// June 2018
*/

#include "coordinate_utilities.h"

#define RES_EPS 1E-11

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
  pointerEr(needle->ex);
  
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
  pointerEr(needle->in);
  
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
  pointerEr(needle->guess);
  
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
      abortEr("This point has been found twice in a same patch.\n"
      "Apparently, some part of point_finder is wrong or the needle"
      "has not been initialized correctly.\n");
    i++;
  }
  
  needle->ans = 
    realloc(needle->ans,(needle->Nans+1)*sizeof(*needle->ans));
  pointerEr(needle->ans);
  
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
    abortEr("There is no such mode.\n");
  
  for (i = 0; i < np; i++)
  {
    double X[3];
    double lim[TOT_Limit];
    Patch_T *patch = needle->grid->patch[p[i]];
    int a;
    
    a = X_of_x(X,needle->x,patch);
    
    fill_limits(lim,patch);
    
    if (a && IsInside(X,lim))
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
    r = X_of_x_Cartesian_coord(X,x);
  else if (patch->coordsys == CubedSpherical)
    r = X_of_x_CS_coord(X,x,patch,1);
  else
      abortEr("No finder for this coordinate.\n");
 
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
      abortEr(NO_JOB);
 
  return ret;
}

/* find x in cartesian coord correspond to X (general coords) 
// for Cartesian coord. Note: x reported with respect to the origin (0,0,0)
// ->return value 1 if it is successful, otherwise 0. */
static int x_of_X_Cartesian_coord(double *const x,const double *const X,const Patch_T *const patch)
{
  UNUSED(patch);  
  
  x[0] = X[0];
  x[1] = X[1];
  x[2] = X[2];
  
  return 1;
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
               xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2,
                R1 = patch->CoordSysInfo->CubedSphericalCoord->R1,
                R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
  const double *const C = patch->c;/* center of origine translated */
  double x1,x2,d,L;
  double x_test[3],X_test[3],dX;
  
  SignAndIndex_permutation_CubedSphere(side,&a,&b,&c,&S);

  switch (type)
  {
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
      L  = 1.-S*d*xc1/R2;
      
      x[c] = x1/(1-L*X[2]);
      x[a] = S*x[c]*X[0];
      x[b] = S*x[c]*X[1];
      
      x[a]+= C[a];
      x[b]+= C[b];
      x[c]+= C[c];
    break;
    case OT_T2_CS:
      L = 1.-R1/R2;
      d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = S*R1/d;
      
      x[c] = x1/(1.-L*X[2]);
      x[a] = S*x[c]*X[0];
      x[b] = S*x[c]*X[1];
      
      x[a]+= C[a];
      x[b]+= C[b];
      x[c]+= C[c];
    break;
    default:
      abortEr(NO_OPTION);
  }
  
  /* test the solution */
  if (check_flg)
  {
    x_test[0] = x[0];
    x_test[1] = x[1];
    x_test[2] = x[2];
    X_of_x_CS_coord(X_test,x_test,patch,0);
    dX = rms(3,X,X_test);
    if (!EQL(dX,0))
      return 0;
  }
  
  return 1;
}

/* find point X correspond to x for patch with Cartesian coord.
// ->return value: 1 if it is successful, otherwise 0.
*/
static int X_of_x_Cartesian_coord(double *const X,const double *const x)
{
  X[0] = x[0];
  X[1] = x[1];
  X[2] = x[2];
  
  return 1;
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
               xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2,
                R1 = patch->CoordSysInfo->CubedSphericalCoord->R1,
                R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
  double x1,x2,d,L;
  double x_test[3],X_test[3],dx;
 
  SignAndIndex_permutation_CubedSphere(side,&i,&j,&k,&S);
  
  X[0] = S*x[i]/x[k];
  X[1] = S*x[j]/x[k];
  
  switch (type)
  {
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
      L = 1.-S*d*xc1/R2;
      X[2] = (1-xc1/x[k])/L;
    break;
    case OT_T2_CS:
      L = 1.-R1/R2;
      d  = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      x1 = S*R1/d;
      X[2] = (1-x1/x[k])/L;
    break;
    default:
      abortEr(NO_OPTION);
  }
  
  /* adujusting boundary number to avoid some unexpeted behavior
  // due to round off error. */
  if (EQL(X[0],1))  X[0] = 1;
  if (EQL(X[0],-1)) X[0] = -1;
  
  if (EQL(X[1],1))  X[1] = 1;
  if (EQL(X[1],-1)) X[1] = -1;
  
  if (EQL(X[2],1))  X[2] = 1;
  if (EQL(X[2],0))  X[2] = 0;
  
  /* test the solution */
  if (check_flg)
  {
    X_test[0] = X[0];
    X_test[1] = X[1];
    X_test[2] = X[2];
    x_of_X_CS_coord(x_test,X_test,patch,0);
    dx = rms(3,cart,x_test);
    double scale = MaxMag_d(rms(3,cart,0),rms(3,x_test,0));
    scale = scale < 1 ? 1 : scale;
    if (!EQL(dx/scale,0))
      return 0;
  }
  
  return 1;
}

/* fill limits based on patch boundary */
static void fill_limits(double *const lim, const Patch_T *const patch)
{
  lim[MIN0] = patch->min[0];
  lim[MIN1] = patch->min[1];
  lim[MIN2] = patch->min[2];
  lim[MAX0] = patch->max[0];
  lim[MAX1] = patch->max[1];
  lim[MAX2] = patch->max[2];
}

/* if x occurs inside the limits (boundary of a patch).
// ->return value : 1 if yes, 0 otherwise.
*/
static int IsInside(const double *const x,const double *const lim)
{
  int v = 0;
  
  if (
      LSSEQL(x[0],lim[MAX0]) && GRTEQL(x[0],lim[MIN0]) &&
      LSSEQL(x[1],lim[MAX1]) && GRTEQL(x[1],lim[MIN1]) &&
      LSSEQL(x[2],lim[MAX2]) && GRTEQL(x[2],lim[MIN2])
      )
    v = 1;
  else
    v = 0;
    
  return v;
}

/* given point and patch find if the is any node collocated 
// to that point and then return its index.
// ->return value: found index, and put flg = FOUND, otherwise, flg = NONE.
*/
unsigned find_node(const double *const x, const Patch_T *const patch,Flag_T *const flg)
{
  unsigned v = UINT_MAX;
  double res = RES_EPS*rms(3,x,0);/* resolution */
  unsigned i;
  double *y, nrm;

  res = GRT(res,RES_EPS) ? res : RES_EPS;
  *flg = NONE;
    
  FOR_ALL(i,patch->node)
  {
    y = patch->node[i]->x;
    nrm = rms(3,x,y);
    if (LSSEQL(nrm,res))
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
  pointerEr(needle);

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

