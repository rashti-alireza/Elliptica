/*
// Alireza Rashti
// June 2018
*/

#include "useful_functions.h"

#define RES_EPS 1E-11

/* find the point in patches given specified needle and fill
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

/* find point X correspond to x in given patch.
// ->return value 1 if it is successful, otherwise 0.
*/
int X_of_x(double *const X,const double *const x,const Patch_T *const patch)
{
  int r = 0;
  
  if (strcmp_i(patch->coordsys,"Cartesian"))
    r = X_of_x_Cartesian_coord(X,x);
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
