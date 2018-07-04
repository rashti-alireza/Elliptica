/*
// Alireza Rashti
// June 2018
*/

#include "useful_functions.h"

/* helps you to find where little tests start*/
void test_start(const char *const file,const int line)
{
  printf("Test starts at\n"
    "File:%s\nLine:%d\n",file,line);
}

/* count the number of pointers which end to null 
// excluding the last one which is null
*/
unsigned countf(void *const p)
{
  assert(p != 0);
  
  unsigned c = 0;
  void **const pp = p;
  
  while (pp[c] != 0)
    c++;
    
  return c;
}

/* linear format to triple (i,j,k) format */
void IJK(const unsigned l, const unsigned *const n, unsigned *const i, unsigned *const j, unsigned *const k)
{
  unsigned tmp;
  
  tmp = l % (n[2]*n[1]);
  *i  = l / (n[2]*n[1]);
  *j  = tmp / n[2];
  *k  = tmp % n[2];
}

/* triple (i,j,k) format to linear format */
unsigned L(const unsigned *const n, const unsigned i, const unsigned j, const unsigned k)
{
  return (k+n[2]*(j+n[1]*i));
}

/* linear format to i component */
unsigned I(const unsigned l, const unsigned *const n)
{
  return l / (n[2]*n[1]);
}

/* linear format to j component */
unsigned J(const unsigned l, const unsigned *const n)
{
  unsigned tmp;
  
  tmp = l % (n[2]*n[1]);
  return tmp / n[2];
}

/* linear format to k component */
unsigned K(const unsigned l, const unsigned *const n)
{
  unsigned tmp;
  
  tmp = l % (n[2]*n[1]);
  return tmp % n[2];
}

/* changing text to enum for collocation */
Collocation_T get_collocation(const char *const coll)
{
  Collocation_T c = NDEF_COLLOCATION;
  if (strcmp_i(coll,"EquiSpaced")) c = EquiSpaced;
  else if (strcmp_i(coll,"Chebyshev_Zero")) c = Chebyshev_Zero;
  else
    abortEr_s("There is no such %s collocation.\n",coll);
    
  return c;
}

/* find out if this point p located on an edge or not
// the algorithm is simple; if it happens at two or more interfaces,
// it means this point is on an edge and returns 1 otherwise 0
// note: it's only for points which are collocated
// ->return value: 1 if found 0 if not.
*/
int IsOnEdge(const unsigned *const n,const unsigned p)
{
  unsigned i,j,k;
  int c;
  
  IJK(p,n,&i,&j,&k);
  
  c = 0;
  if (i == n[0]-1 || i == 0)  c++;
  if (j == n[1]-1 || j == 0)  c++;
  if (k == n[2]-1 || k == 0)  c++;
  
  if (c > 1)  return 1;
  
  return 0;
}

/* find out if this point p located on an face or not.
// the algorithm is simple; if it happens at two or more interfaces
// it means this point is on an face and returns number of face 
// otherwise 0.
// moreover, the found face f written like the example below:
// f[I_0] = 1, the point happens at face I_0,
// f[J_n1] = 0, the point won't happen at face J_n1 and etc.
// note: it's a general function and works for points which are not 
// collocated as well.
// ->return value: number of interface that found.
*/
int IsOnFace(const double *const x, const Patch_T *const patch,unsigned *const f)
{
  int u,c;
  double X[3];
  
  for (u = 0; u < TOT_FACE; u++)
    f[u] = 0;
  
  X_of_x(X,x,patch);
  
  c = 0;
  for (u = 0; u < TOT_FACE; u++)
  {
    f[u] = check_interface(X,patch,u);
    
    if (f[u] == 1) c++;
  }
  
  return c;
}

/* check if X is on specific interface.
// ->return value: 1 if yes, 0 otherwise.
*/
static unsigned check_interface(const double *const X, const Patch_T *const patch, const int u)
{
  
  const unsigned ind = patch->interface[u]->point[0]->ind;
  double res = EPS*rms(3,X,0);
  double *Y;
  
  res = GRT(res,EPS) ? res: EPS;
  
  if (patch->node[ind]->X != 0) Y = patch->node[ind]->X;
  else                          Y = patch->node[ind]->x;
  
  switch (u)
  {
    case I_0:
      if (LSSEQL(ABS(Y[0]-X[0]),res)) return 1;
      break;
    case I_n0:
      if (LSSEQL(ABS(Y[0]-X[0]),res)) return 1;
      break;
    case J_0:
      if (LSSEQL(ABS(Y[1]-X[1]),res)) return 1;
      break;
    case J_n1:
      if (LSSEQL(ABS(Y[1]-X[1]),res)) return 1;
      break;
    case K_0:
      if (LSSEQL(ABS(Y[2]-X[2]),res)) return 1;
      break;
    case K_n2:
      if (LSSEQL(ABS(Y[2]-X[2]),res)) return 1;
      break;
    default:
      abortEr("No such interface was defined for this function.\n");
      
  }
  
  return 0;
}

/* given Cartesian value of a point, patch and f, it finds the closest
// node to this given point and return its index.
// ->return value: node index.
*/
unsigned node_onFace(const double *const x, const unsigned f,const Patch_T *const patch)
{
  const Interface_T *const face = patch->interface[f];
  Node_T **const node = patch->node;
  double nrm, s = DBL_MAX;
  unsigned ind = UINT_MAX,i;
  
  FOR_ALL(i,face->point)
  {
    unsigned id = face->point[i]->ind;
    nrm = rms(3,x,node[id]->x);
    if (LSS(nrm,s))
    {
      s = nrm;
      ind = id;
    }
  }
  
  return ind; 
}

/* getting the subface which paired with sub
// ->return value: found paired subface.
*/
SubFace_T *get_paired_subface(const SubFace_T *const sub)
{
  const Patch_T *const patch = sub->patch->grid->patch[sub->adjPatch];
  const Interface_T *const face = patch->interface[sub->adjFace];
  
  assert(face);
  
  return face->subface[sub->adjsn];
  
}