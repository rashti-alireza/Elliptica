/*
// Alireza Rashti
// June 2018
*/

#include "little_functions.h"

/* helps you to find where litte tests start*/
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
