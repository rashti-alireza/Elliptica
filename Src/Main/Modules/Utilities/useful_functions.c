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

/* changing text to enum for collocation 
// ->return value: collocation type, error if couldn't be found, 0 if null coll
*/
Collocation_T get_collocation(const char *const coll)
{
  Collocation_T c = UNDEFINED_COLLOCATION;
  
  if (!coll) return c;
  else if (strcmp_i(coll,"EquiSpaced")) c = EquiSpaced;
  else if (strcmp_i(coll,"Chebyshev_Extrema")) c = Chebyshev_Extrema;
  else if (strcmp_i(coll,"Chebyshev_Nodes")) c = Chebyshev_Nodes;
  else
    abortEr_s("There is no such %s collocation.\n",coll);
    
  return c;
}

/* changing text to enum for basis 
// ->return value: basis type, error if couldn't be found, 0 if null basis
*/
Basis_T get_basis(const char *const basis)
{
  Basis_T b = UNDEFINED_BASIS;
  
  if (!basis) return b;
  else if (strcmp_i(basis,"Chebyshev_FirstKind")) 
    b = Chebyshev_Tn_BASIS;
  else if (strcmp_i(basis,"No_Basis")) 
    b = No_BASIS;
  else
    abortEr_s("There is no such %s basis.\n",basis);
    
  return b;
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
// the algorithm is simple; if it happens at one or more interfaces
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
  Y = patch->node[ind]->X;
  
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

/* ->return value: total number of nodes on the given grid */
unsigned total_nodes_grid(const Grid_T *const grid)
{
  unsigned pa;
  unsigned sum = 0;
  
  FOR_ALL(pa,grid->patch)
    sum += total_nodes_patch(grid->patch[pa]);
  
  return sum;
}

/* ->return value: total number of nodes on the given patch */
unsigned total_nodes_patch(const Patch_T *const patch)
{
  return patch->n[0]*patch->n[1]*patch->n[2];
}

/* find coord enum based in give str.
// ->return value: found Coord_T .
*/
Coord_T find_coord(const char *const coordsys)
{
  Coord_T coord;
  
  if(strcmp_i(coordsys,"Cartesian"))
    coord = Cartesian;
  else
    abortEr_s("There is no such %s coordinates.\n",coordsys);  
    
  return coord;
}

/* generating a number number within final-initial.
// s determines weather srand be called or not.
// in fact, only when s is zero srand is called. and the reason is to avoid
// producing same random number in a loop.
// ->return value: random number in the double data type within the range.
*/
double random_double(const double initial,const double final,const unsigned s)
{
  time_t t;
  
  /* Initializes random number generator */
  if (s == 0)
    srand((unsigned) time(&t));
   
  return initial+((final-initial)/RAND_MAX)*rand();
}

/* hard copy subface s1 to subface s2 */
void copy_subface(SubFace_T *const s2,const SubFace_T *const s1)
{
  unsigned i;

  s2->patch = s1->patch;
  s2->flags_str = dup_s(s1->flags_str);
  s2->sn = s1->sn;
  s2->adjsn = s1->adjsn;
  s2->np = s1->np;
  s2->id = calloc(s2->np,sizeof(*s2->id));
  pointerEr(s2->id);
  s2->adjid = calloc(s2->np,sizeof(*s2->adjid));
  pointerEr(s2->adjid);
  for (i = 0; i < s2->np; ++i)
    s2->id[i] = s1->id[i];
    
  if (s1->copy)
    for (i = 0; i < s2->np; ++i)
      s2->adjid[i] = s1->adjid[i];
      
  s2->face     = s1->face;
  s2->adjFace  = s1->adjFace;
  s2->adjPatch = s1->adjPatch;
  s2->df_dn  = s1->df_dn;
  s2->sameX  = s1->sameX;
  s2->sameY  = s1->sameY;
  s2->sameZ  = s1->sameZ;
  s2->touch  = s1->touch;
  s2->copy   = s1->copy;
  s2->exterF = s1->exterF;
  s2->outerB = s1->outerB;
  s2->innerB = s1->innerB;
}

/* map: points on each subface has an index for example id[i] = n.
// which says point with index i correspond to node n.
// this function gets n and return i mean invers(id[i]) = i;
// ->return value : invers(id[i])
*/
unsigned subface_map_invers_id(const SubFace_T *const subface,const unsigned n)
{
  unsigned i;
  unsigned s = UINT_MAX;
  
  for (i = 0; i < subface->np; ++i)
    if (subface->id[i] == n)
      return i;
  
  return s;  
}

/* given array s and its dimension, duplicate the array.
// ->return value: duplicated unsigned array. */
unsigned *dup_UINT(const unsigned *const s,const unsigned N)
{
  if (!s)
    abortEr("The given array to be duplicated is Null!");
  
  unsigned *dup = calloc(N,sizeof(*dup));
  pointerEr(dup);
  unsigned i;
  
  for (i = 0; i < N; ++i)
    dup[i] = s[i];
    
  return dup;
}

/* calculating the magnetude of d(X,Y,Z)/d(x,y,z) */
double max_Jacobian_dX_dx(Patch_T *const patch)
{
  double max = 0;
  double abs_j;
  const unsigned nn = patch->nn;
  unsigned l; 
  
  for (l = 0; l < nn; ++l)
  {
    abs_j = dq2_dq1(patch,_a_,_x_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_a_,_y_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_a_,_z_,l);
    if (abs_j > max)
      max = abs_j;
      
    abs_j = dq2_dq1(patch,_b_,_x_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_b_,_y_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_b_,_z_,l);
    if (abs_j > max)
      max = abs_j;  
      
    abs_j = dq2_dq1(patch,_c_,_x_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_c_,_y_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_c_,_z_,l);
    if (abs_j > max)
      max = abs_j;
      
  }
  
  return max;
}

/* calculating the maximum error of spectral derivative using the fact 
// that computer is using finite number of digits.
// f : given field
// o : given order of derivative
// ->return value: error in calculation = general idea is as follow:
// 1e-14*max(func)*max(Jacobian)^(order of derivative )*n*n^(2*order of derivative)*10  */
double spectral_derivative_max_error(const Field_T *const f,const unsigned o)
{
  double e = 1e-15;
  double max_f,
         max_j;
  unsigned max_n;
  const unsigned *const n = f->patch->n;
  const char *der_par = GetParameterS("Derivative_Method");
  
  if (strstr_i(der_par,"Spectral"))
  {
    max_f = L_inf(n[0]*n[1]*n[2],f->v);
    max_j = max_Jacobian_dX_dx(f->patch);
    max_n =  n[0] > n[1] ? n[0]  : n[1];
    max_n = max_n > n[2] ? max_n : n[2];
    e = 1e-14*max_n*pow(max_n,2*o);/* 1e-14 coming from: (machine precision=1e-15)*(all other unwarranted factors=10) */
    
    if (max_j > 1.0)/* if max_j is less than 1 it won't gonna make the error lesser! */
      e *= pow(max_j,o);
    if (max_f > 1.0)/* if max_f is less than 1 it won't gonna make the error lesser! */
      e *= max_f;
    if (e == 0.0) 
      e = 1e-15;
  } 
  else
    abortEr(NO_JOB);
    
  return e;
}

/* ->return value: if the patch covers a part of the BH horizon 1, otherwise 0 */
unsigned IsItHorizonPatch(const Patch_T *const patch)
{
  unsigned ret = 0;
  
  if (strcmp_i(patch->grid->kind,"BBN_CubedSpherical_grid"))
  {
    if (regex_search("_right_BH_surrounding_.+",patch->name))
      return 1;
  }
  else if (strcmp_i(patch->grid->kind,"SNS_CubedSpherical+Box_grid"))
  {
    return 0;
  }
  else
    abortEr(NO_JOB);
  
  return ret;
}

/* ->return value: given stem name of a patch, 
// it finds the corresponding one in the given grid. */
Patch_T *GetPatch(const char *const stem,const Grid_T *const grid)
{
  Patch_T *retPatch = 0;
  char name[1000];
  unsigned p;
  
  sprintf(name,"grid%u_%s",grid->gn,stem);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (strcmp_i(name,patch->name))
    {
      retPatch = patch;
      break;
    }
  }
  
  if(!retPatch)
    abortEr_s("It could not find patch %s.\n",name);
    
  return retPatch;
}

/* ->return value: if the patch covers a part of the NS 1, otherwise 0 */
unsigned IsItNSPatch(const Patch_T *const patch)
{
  unsigned ret = 0;
  
  if (strcmp_i(patch->grid->kind,"BBN_CubedSpherical_grid"))
  {
    /* it could have used regex, but need extra attention to NS surrounding case */
    if (strstr(patch->name,"left_centeral_box") || 
        strstr(patch->name,"left_NS_up")        ||
        strstr(patch->name,"left_NS_down")      ||
        strstr(patch->name,"left_NS_back")      ||
        strstr(patch->name,"left_NS_front")     ||
        strstr(patch->name,"left_NS_left")      ||
        strstr(patch->name,"left_NS_right")     
       )
       ret = 1;
  }
  else if (strcmp_i(patch->grid->kind,"SNS_CubedSpherical+Box_grid"))
  {
    /* it could have used regex, but need extra attention to NS surrounding case */
    if (strstr(patch->name,"left_centeral_box") || 
        strstr(patch->name,"left_NS_up")        ||
        strstr(patch->name,"left_NS_down")      ||
        strstr(patch->name,"left_NS_back")      ||
        strstr(patch->name,"left_NS_front")     ||
        strstr(patch->name,"left_NS_left")      ||
        strstr(patch->name,"left_NS_right")     
       )
       ret = 1;
  }
  else
    abortEr(NO_JOB);
    
  return ret;
}

/* ->return value: if the patch covers a part of the NS surface 1, otherwise 0 */
unsigned IsItNSSurface(const Patch_T *const patch)
{
  unsigned ret = 0;
  
  if (strcmp_i(patch->grid->kind,"BBN_CubedSpherical_grid"))
  {
    /* it could have used regex, but need extra attention to NS surrounding case */
    if (strstr(patch->name,"left_NS_up")        ||
        strstr(patch->name,"left_NS_down")      ||
        strstr(patch->name,"left_NS_back")      ||
        strstr(patch->name,"left_NS_front")     ||
        strstr(patch->name,"left_NS_left")      ||
        strstr(patch->name,"left_NS_right")     
       )
       ret = 1;
  }
  else if (strcmp_i(patch->grid->kind,"SNS_CubedSpherical+Box_grid"))
  {
    /* it could have used regex, but need extra attention to NS surrounding case */
    if (strstr(patch->name,"left_NS_up")        ||
        strstr(patch->name,"left_NS_down")      ||
        strstr(patch->name,"left_NS_back")      ||
        strstr(patch->name,"left_NS_front")     ||
        strstr(patch->name,"left_NS_left")      ||
        strstr(patch->name,"left_NS_right")     
       )
       ret = 1;
  }
  else
    abortEr(NO_JOB);
    
  return ret;
}

/* ->return value: if the patch is one of the surrounding patches of the NS 1, otherwise 0 */
unsigned IsItNSSurroundingPatch(const Patch_T *const patch)
{
  unsigned ret = 0;
  
  if (strcmp_i(patch->grid->kind,"BBN_CubedSpherical_grid"))
  {
    if (strstr(patch->name,"left_NS_surrounding_up")        ||
        strstr(patch->name,"left_NS_surrounding_down")      ||
        strstr(patch->name,"left_NS_surrounding_back")      ||
        strstr(patch->name,"left_NS_surrounding_front")     ||
        strstr(patch->name,"left_NS_surrounding_left")      ||
        strstr(patch->name,"left_NS_surrounding_right")     
       )
       ret = 1;
  }
  else if (strcmp_i(patch->grid->kind,"SNS_CubedSpherical+Box_grid"))
  {
    if (strstr(patch->name,"left_NS_surrounding_up")        ||
        strstr(patch->name,"left_NS_surrounding_down")      ||
        strstr(patch->name,"left_NS_surrounding_back")      ||
        strstr(patch->name,"left_NS_surrounding_front")     ||
        strstr(patch->name,"left_NS_surrounding_left")      ||
        strstr(patch->name,"left_NS_surrounding_right")     
       )
       ret = 1;
  }
  else
    abortEr(NO_JOB);
    
  return ret;
}

/* print an array of double type with dimension n, for debuging purposes */
void dbprint(const double *v,const unsigned n,const char *const desc)
{
  unsigned i;
  
  pr_line();
  
  if (desc)/* description */
    printf("Debugging: %s\n",desc);
  else
    printf("Debugging:\n");
  
  for (i = 0; i < n; ++i)
    printf("data[%02u] = %+.15f\n",i,v[i]);
  
  pr_line();
}
