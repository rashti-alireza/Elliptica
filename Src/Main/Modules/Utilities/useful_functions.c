/*
// Alireza Rashti
// June 2018
*/

#include "useful_functions.h"

/* helps you to find where little tests start*/
void test_start(char *file,int line)
{
  printf("Test starts at\n"
    "File:%s\nLine:%d\n",file,line);
}

/* count the number of pointers which end to null 
// excluding the last one which is null
*/
int countf(void *p)
{
  assert(p != 0);
  
  int c = 0;
  void **pp = p;
  
  while (pp[c] != 0)
    c++;
    
  return c;
}

/* linear format to triple (i,j,k) format */
void IJK(int l, int *n, int *i, int *j, int *k)
{
  int tmp;
  
  tmp = l % (n[2]*n[1]);
  *i  = l / (n[2]*n[1]);
  *j  = tmp / n[2];
  *k  = tmp % n[2];
}

/* triple (i,j,k) format to linear format */
int L(int *n, int i, int j, int k)
{
  return (k+n[2]*(j+n[1]*i));
}

/* linear format to i component */
int I(int l, int *n)
{
  return l / (n[2]*n[1]);
}

/* linear format to j component */
int J(int l, int *n)
{
  int tmp;
  
  tmp = l % (n[2]*n[1]);
  return tmp / n[2];
}

/* linear format to k component */
int K(int l, int *n)
{
  int tmp;
  
  tmp = l % (n[2]*n[1]);
  return tmp % n[2];
}

/* changing text to enum for collocation */
Collocation_T get_collocation(char *coll)
{
  Collocation_T c;
  if (strcmp_i(coll,"EquiSpaced")) c = EquiSpaced;
  else if (strcmp_i(coll,"Chebyshev_Zero")) c = Chebyshev_Zero;
  else
    abortEr_s("There is no such %s collocation.\n",coll);
    
  return c;
}

/* find out if this point p located on an edge or not
// the algorithm is simple; if it happens at two or more interface
// it means this point is on an edge and returns 1 otherwise 0
 */
int IsThisEdge(int *n,int p)
{
  int i,j,k;
  int c, r;
  
  IJK(p,n,&i,&j,&k);
  
  c = 0;
  if (i == n[0]-1 || i == 0)  c++;
  if (j == n[1]-1 || j == 0)  c++;
  if (k == n[2]-1 || k == 0)  c++;
  
  r = 0;
  if (c > 1)  r = 1;
  
  return r;
}
