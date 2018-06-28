/*
// Alireza Rashti
// June 2018
*/

#include "useful_functions.h"

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
    look(needle,GUESS);
    if (needle->Nans > 0) return;
  }// end of if (needle->Ng != 0)
  
  /* it only looks inside include patches if found it gets out */
  if (needle->Nin > 0)
  {
    look(needle,FORCE_IN);
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
    
    look(needle,FORCE_IN);
  }// end of if (needle->Nex > 0)
  
  else // if non of above met
  {
    FOR_ALL(i,needle->grid->patch)
      needle_in(needle,needle->grid->patch[i]);
      
    look(needle,FORCE_IN);
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

/* check if in != ex and in != guess */
static void IsConsistent(Needle_T *needle)
{
  int i,j;
  
  if (needle->in == 0) return;
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

/* look for point in designated patches inside needle based on mode */
static void look(Needle_T *needle,Mode_T mode)
{
  pr_line();
}