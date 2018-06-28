/*
// Alireza Rashti
// June 2018
*/

#include "useful_functions.h"

#define RES_EPS 1E-11

//-----------------------------------------------------------------//
/* find the point in patches which are specified needle and fill
// the needle with the answers.
// for more information about needle look at the typede_data.h in 
// Core folder. for an example to how make needle look at grid.c and
// search for Needle_T or point_finder.
// note: don't forget to free at the end of the day the answer of
/////////////////////////////////////////////////////////////////////
typedef struct NEEDLE_T
{
  double *x;
  Grid_T *grid;// the grid which is used
  int *guess;// these are guess patches searched firstly
  int *in;// force it to look only inside these patches.
          // notation: in[?] = patch number
  int *ex;// force it to not look inside these patches
           // notation: ex[?] = patch number
  int *ans;// the answers found which pointing to patch number
  int Nans;// number of answer
  int Nin;// number of included patches
  int Nex;// number of excluded patches
  int Ng;// numbef of guess patches
}Needle_T;
/////////////////////////////////////////////////////////////////////
*/
void point_finder(Needle_T *needle)
{
  int i;
  
  /* check consistency between in, ex and guess */
  IsConsistent(needle);
  
  /* it only looks inside guess patches if found it gets out */
  if (needle->Ng != 0)
  {
    find(needle,GUESS);
    if (needle->Nans > 0) return;
  }// end of if (needle->Ng != 0)
  
  /* it only looks inside include patches if found it gets out */
  if (needle->Nin > 0)
  {
    find(needle,FORCE_IN);
    if (needle->Nans > 0) return;
  }// end of if (needle->Nin > 0)
  
  /* find in all patched exluding needle->ex */
  if (needle->Nex > 0)
  {
    int j;
    Flag_T flg = NO;
    
    FOR_ALL(i,needle->grid->patch)
    {
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
  }// end of if (needle->Nex > 0)
  
  else // if non of above met
  {
    FOR_ALL(i,needle->grid->patch)
      needle_in(needle,needle->grid->patch[i]);
      
    find(needle,FORCE_IN);
  }
}

/* adding a patch to needle->ex */
void needle_ex(Needle_T *needle,Patch_T *patch)
{
  assert(needle);
  int i;
  
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
}

/* adding a patch to needle->in */
void needle_in(Needle_T *needle,Patch_T *patch)
{
  assert(needle);
  int i;
  
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
}

/* adding a patch to needle->guess */
void needle_guess(Needle_T *needle,Patch_T *patch)
{
  assert(needle);
  int i;
  
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
}

/* adding a patch to needle->ans */
void needle_ans(Needle_T *needle,Patch_T *patch)
{
  assert(needle);
  int i;
  
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
  
  needle->guess[needle->Nans] = patch->pn;
}

/* check if in != ex, ex !=guess */
static void IsConsistent(Needle_T *needle)
{
  int i,j;
  
  if (needle->in == 0)
  {
    for (i = 0; i < needle->Nex; i++)
      for (j = 0; j < needle->Ng; j++)
        if (needle->ex[i] == needle->guess[j])
          abortEr("Guess and Exclude must be mutually exclusive.\n");
  }
  else
  {
    for (i = 0; i < needle->Nex; i++)
    {
      for (j = 0; j < needle->Nin; j++)
        if (needle->ex[i] == needle->in[j])
          abortEr("Include and Exclude must be mutually exclusive.\n");
      for (j = 0; j < needle->Ng; j++)
        if (needle->ex[i] == needle->guess[j])
          abortEr("Guess and Exclude must be mutually exclusive.\n");
    }
  }
}

/* find for point in designated patches inside needle based on mode */
static void find(Needle_T *needle,Mode_T mode)
{
  int *p, np;
  int i;
  
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

/* find point X correspond to x in given patch.
// ->return value 1 if it is successful, otherwise 0.
*/
int X_of_x(double *X,double *x,Patch_T *patch)
{
  int r = 0;
  
  if (strcmp_i(patch->coordsys,"Cartesian"))
    r = X_of_x_Cartesian_coord(X,x,patch);
  /* other coord sys comes here
  .
  .
  .
  */
  else
      abortEr_s("No finder for this %s coordinate.\n",patch->coordsys);
 
  return r;
}

/* find point X correspond to x for patch with Cartesian coord.
// ->return value: 1 if it is successful, otherwise 0.
*/
int X_of_x_Cartesian_coord(double *X,double *x,Patch_T *patch)
{
  X[0] = x[0];
  X[1] = x[1];
  X[2] = x[2];
  
  return 1;
}

/* fill limits based on patch boundary */
static void fill_limits(double *lim, Patch_T *patch)
{
  lim[MIN0] = patch->min[0];
  lim[MIN1] = patch->min[1];
  lim[MIN2] = patch->min[2];
  lim[MAX0] = patch->max[0];
  lim[MAX1] = patch->max[1];
  lim[MAX2] = patch->max[2];
}

/* if x occurs inside the limits (boundary) return 1 otherwise 0 */
static int IsInside(double *x,double *lim)
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
// to that point and then return its index. otherwise return -1.
*/
int find_node(double *x, Patch_T *patch)
{
  int v = -1;
  double X[3];
  const int r = X_of_x(X,x,patch);
  
  if (r)
  {
    const double res = RES_EPS*rms(3,X,0);// resolution
    int i;
    double *y, nrm;
    
    if (strcmp_i(patch->coordsys,"Cartesian"))
    {
      FOR_ALL(i,patch->node)
      {
        y = patch->node[i]->x;
        nrm = rms(3,x,y);
        if (LSSEQL(nrm,res))
        {
          v = i;
          break;
        }
      }
    }
    else
    {
      FOR_ALL(i,patch->node)
      {
        y = patch->node[i]->X;
        nrm = rms(3,x,y);
        if (LSSEQL(nrm,res))
        {
          v = i;
          break;
        }
      }
    }
  }
  
  return v;
}
