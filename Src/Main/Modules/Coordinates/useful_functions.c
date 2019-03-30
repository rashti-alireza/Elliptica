/*
// Alireza Rashti
// June 2018
*/

#include "useful_functions.h"

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
  /* other coord sys comes here
  .
  .
  .
  */
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
  
  if (patch->coordsys == ProjectiveHemisphereUp)
    ret = x_of_X_PHUp_coord(x,X,patch);
  else if (patch->coordsys == ProjectiveHemisphereDown)
    ret = x_of_X_PHDown_coord(x,X,patch);
  else if (patch->coordsys == StereographicSphereLeft)
    ret = x_of_X_SSLeft_coord(x,X,patch);
  else if (patch->coordsys == StereographicSphereRight)
    ret = x_of_X_SSRight_coord(x,X,patch);
  else
      abortEr(NO_JOB);
 
  return ret;
}

/* find x in cartesian coord correspond to X (general coords) 
// for ProjectiveHemisphereUp 
// ->return value 1 if it is successful, otherwise 0.
*/
static int x_of_X_PHUp_coord(double *const x,const double *const X,const Patch_T *const patch)
{
  Field_T *const R1_field = patch->pool[Ind("R1_ProjectiveHemisphere")];
  Field_T *const R2_field = patch->pool[Ind("R2_ProjectiveHemisphere")];
  const double R1 = interpolation_2d_PH(R1_field,patch,X);
  const double R2 = interpolation_2d_PH(R2_field,patch,X); 
  const double r = 0.5*((R2-R1)*X[2]+(R2+R1));
  const double *const c = patch->c;/* center */
  double z2;
  int ret = 0;
  
  x[0] = r*X[0]*sqrt(1-0.5*SQR(X[1]));
  x[1] = r*X[1]*sqrt(1-0.5*SQR(X[0]));
  z2 = SQR(r)-SQR(x[0])-SQR(x[1]);
  if (EQL(z2,0))
    z2 = 0;/* avoid nan */
  x[2] = sqrt(z2);
  
  if (isnan(x[0]) != 0 || 
      isnan(x[1]) != 0 || 
      isnan(x[2]) != 0   )
  {
    ret = 0;
    abortEr("x could not been found.");
  }
  
  x[0] += c[0];
  x[1] += c[1];
  x[2] += c[2];
  
  return ret;
}

/* find x in cartesian coord correspond to X (general coords) 
// for StereographicSphereRight
// ->return value 1 if it is successful, otherwise 0.
*/
static int x_of_X_SSRight_coord(double *const x,const double *const X,const Patch_T *const patch)
{
  const double R1 = patch->CoordSysInfo->R1;
  const double R2 = patch->CoordSysInfo->R2;
  const double R0 = patch->c[1];
  double u,w;
  double r,R,A,c;
  
  r = 0.5*X[1]*(R2-R1)+0.5*(R2+R1);
  R = sqrt(SQR(r)-SQR(R0)); assert(!isnan(R));
  u = R*X[0]*sqrt(1-0.5*SQR(X[2])); assert(!isnan(u));
  w = R*X[2]*sqrt(1-0.5*SQR(X[0])); assert(!isnan(w));
  A = SQR(u/(R0-r))+SQR(w/(R0-r))+1;
  x[1] = r*(2/A-1)+R0;
  c = 2*r/(A*(r-R0));
  x[0] = c*u;
  x[2] = c*w;
  
  return 1;
}

/* find x in cartesian coord correspond to X (general coords) 
// for StereographicSphereLeft
// ->return value 1 if it is successful, otherwise 0.
*/
static int x_of_X_SSLeft_coord(double *const x,const double *const X,const Patch_T *const patch)
{
  const double R1 = patch->CoordSysInfo->R1;
  const double R2 = patch->CoordSysInfo->R2;
  const double R0 = fabs(patch->c[1]);
  double u,w;
  double r,R,A,c;
  
  r = 0.5*X[1]*(R2-R1)+0.5*(R2+R1);
  R = sqrt(SQR(r)-SQR(R0)); assert(!isnan(R));
  u = R*X[0]*sqrt(1-0.5*SQR(X[2])); assert(!isnan(u));
  w = R*X[2]*sqrt(1-0.5*SQR(X[0])); assert(!isnan(w));
  A = SQR(u/(R0-r))+SQR(w/(R0-r))+1;
  x[1] = -r*(2/A-1)-R0;
  c = 2*r/(A*(r-R0));
  x[0] = c*u;
  x[2] = c*w;
  
  return 1;
}

/* find x in cartesian coord correspond to X (general coords) 
// for ProjectiveHemisphereDown
// ->return value 1 if it is successful, otherwise 0.
*/
static int x_of_X_PHDown_coord(double *const x,const double *const X,const Patch_T *const patch)
{
  Field_T *const R1_field = patch->pool[Ind("R1_ProjectiveHemisphere")];
  Field_T *const R2_field = patch->pool[Ind("R2_ProjectiveHemisphere")];
  const double R1 = interpolation_2d_PH(R1_field,patch,X);
  const double R2 = interpolation_2d_PH(R2_field,patch,X); 
  const double r = 0.5*((R2-R1)*X[2]+(R2+R1));
  const double *const c = patch->c;/* center */
  double z2;
  int ret = 0;
  
  x[0] = r*X[0]*sqrt(1-0.5*SQR(X[1]));
  x[1] = r*X[1]*sqrt(1-0.5*SQR(X[0]));
  z2 = SQR(r)-SQR(x[0])-SQR(x[1]);
  if (EQL(z2,0))
    z2 = 0;/* avoid nan */
  x[2] = -sqrt(z2);
  
  if (isnan(x[0]) != 0 || 
      isnan(x[1]) != 0 || 
      isnan(x[2]) != 0   )
  {
    ret = 0;
    abortEr("x could not been found.");
  }
  
  x[0] += c[0];
  x[1] += c[1];
  x[2] += c[2];
  
  return ret;
}

/* find point X correspond to x for patch with Cartesian coord.
// ->return value: 1 if it is successful, otherwise 0.
*/
int X_of_x_Cartesian_coord(double *const X,const double *const x)
{
  X[0] = x[0];
  X[1] = x[1];
  X[2] = x[2];
  
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
