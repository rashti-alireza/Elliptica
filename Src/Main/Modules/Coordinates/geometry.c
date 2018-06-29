/*
// Alireza Rashti
// June 2018
*/

#include "geometry.h"

#define EPS 1E-5

/* 
// realizing the geometry of grid such as how patches are glued
// normal vectors at boundary, boundary of grid and etc.
*/
int realize_geometry(Grid_T *grid)
{
  int i;
  
  FOR_ALL(i,grid->patch)
  {
    Patch_T *const patch = grid->patch[i];
    
    alloc_interface(patch);// allocating interfaces
    fill_basics(patch);// filling basic elements
    fill_N(patch);// filling point[?]->N
    flush_houseK(patch);// set all of housK flag to zero
  }
  
  /* find various geometry of each point */
  fill_geometry(grid);
  
  /* printing boundary for test purposes */
  //pr_boundary(grid);
  
  return EXIT_SUCCESS;
}

/* filling the geometry of point struct */
static void fill_geometry(Grid_T *grid)
{
  sFunc_PtoV_T **func;
  int i;
  
  init_func_PtoV(&func);// initialize struct
  /* adding func to struct to be called each coord must have 
  // its own func. note, also external face and inner face must be found
  // in first place; thus, for each new coord sys one must add below 
  // the related functions for these two purposes.
  */
  add_func_PtoV(&func,FindInnerB_Cartesian_coord,"FindInnerB",Cartesian);
  add_func_PtoV(&func,FindExterF_Cartesian_coord,"FindExterF",Cartesian);
  
  /* calling function to find external and internal faces*/
  FOR_ALL(i,grid->patch)
  {
    Patch_T *const patch = grid->patch[i];
    run_func_PtoV(func,"FindExterF",patch);// find external faces
    run_func_PtoV(func,"FindInnerB",patch);// find inner boundary
  }
  free_2d(func);// freeing func struct
  
  /* realize neighbor properties */
  FOR_ALL(i,grid->patch)
    realize_neighbor(grid->patch[i]);
  
}

/* study neighbor of points to find quantities like 
// point->face, point->touch, point->adjPoint and others.
// algorithm:
// for each interface do
// 	a. seprate each interface to its inner and edge points
// 	b. realize the adjacency of inner points
// 	c. realize the adjacency of edge points
// the reason the inner and edge points are separated is they
// use different algorithm, since edge points need be treated specially
// ->return value-> EXIT_SUCCESS
*/
static int realize_neighbor(Patch_T *patch)
{
  const int nf = countf(patch->interface);
  int f;
  
  for (f = 0; f < nf; f++)
  {
    Interface_T *interface = patch->interface[f];
    PointSet_T **innerP, **edgeP;
    init_Points(interface,&innerP,&edgeP);// initializing inner and 
                                          // edge points
    realize_adj(innerP, INNER);// realizing the adjacency INNER one
    realize_adj(edgeP,EDGE);// realizing the adjacency EDGE one
    
    free_PointSet(innerP);
    free_PointSet(edgeP);
  }
  return EXIT_SUCCESS;
}

/* realizing adjacency of given points */
static void realize_adj(PointSet_T **Pnt, enum Type type)
{
  int p;
  
  FOR_ALL(p,Pnt)
  {
    if (Pnt[p]->point->houseK == 1) continue;
    
    find_adjPnt(Pnt[p],type);
    analyze_adjPnt(Pnt[p],type);
  }
}

/* analyzing the adjacent point of point pnt; 
// see which one best encompasses the adjacent boundary of the interface
// and then fill the flags of Point_T.
*/
static void analyze_adjPnt(PointSet_T *Pnt,enum Type type)
{
  Point_T *const p1 = Pnt->point;
  Point_T *p2;
  Adjacent_T *adjp1;
  
  /* this point is on outer boundary */
  if (Pnt->NadjPnt == 0)
  {
    p1->outerB = 1;
  }
  /* if none of the points are on the interface */
  else if (IsOverlap(Pnt))
  {
    p1->touch = 0;
    p1->adjPatch = Pnt->adjPnt[0].p;
    p1->copy  = 0;
  }
  /* which point best meets the normal conditon (N2dotN1 == -1)  */
  else if (IsNormalFit(Pnt))
  {
    adjp1 = &Pnt->adjPnt[Pnt->FitN];
    
    p1->touch = 1;
    p1->copy  = 1;
    p1->adjPatch = adjp1->p;
    p1->adjFace  = adjp1->FitF;
    p2 = get_p2(Pnt);
    assert(p2);
    p1->adjPoint = p2;
    
    p2->touch = 1;
    p2->copy  = 1;
    p2->adjPatch = p1->patch->pn;
    p2->adjFace  = p1->face;
    p2->adjPoint = p1;
    p2->houseK = 1;
  }
  
  p1->houseK = 1;
}

/* find the point p on the interface which has the fittest normal.
// ->return value: a pointer to point struct on interface
*/
static Point_T *get_p2(PointSet_T *Pnt)
{
  int f = Pnt->adjPnt[Pnt->FitN].FitF;
  int p = Pnt->adjPnt[Pnt->FitN].p;
  const int node = Pnt->adjPnt[Pnt->FitN].node;
  Patch_T *const patch = Pnt->point->patch->grid->patch[p];
  Interface_T *const face = patch->interface[f];
  int i;
  
  FOR_ALL(i,face->point)
    if (face->point[i]->ind == node)
      return face->point[i];
  
  return 0;
}

/* if there is an adjacent point collocated with the point,
// and their interfaces are tangent correctly, i.e. their normals
// are in expected directions (N1dotN2 == -1)
// ->retrun value: 1 if fittest normal exists, 0 otherwise
*/
static int IsNormalFit(PointSet_T *Pnt)
{
  int i;
  
  for (i = 0; i < Pnt->NadjPnt; i++)
    if (Pnt->adjPnt[i].FitF >= 0) 
    {
      Pnt->FitN = i;
      return 1;
    }
  
  return 0;
}
/* if none of the adjacent points are on an interface,
// thus, this point overlaps with other patch(es) 
// ->retrun value: 1 if it is overlap, 0 otherwise
*/
static int IsOverlap(PointSet_T *Pnt)
{
  int i;
  
  for (i = 0; i < Pnt->NadjPnt; i++)
    if (Pnt->adjPnt[i].FaceFlg == 1) return 0;
  
  return 1;
}

/* finding the adjacent points of point pnt */
static void find_adjPnt(PointSet_T *pnt,enum Type type)
{
  int ind = pnt->point->ind;
  int *found;// patches have been found
  int nfp;// number of found patches
  double *x = pnt->point->patch->node[ind]->x;
  Needle_T *needle = alloc_needle();
  needle->grid = pnt->point->patch->grid;
  needle_ex(needle,pnt->point->patch);
  needle->x = x;
  point_finder(needle);// first find the point itself
  
  if (type == INNER)
  {
    double *N = pnt->point->N;
    double eps = rms(3,x,0)*EPS;
    double q[3] = {x[0]+eps*N[0],
                   x[1]+eps*N[1],
                   x[2]+eps*N[2]};// q = pnt+eps*N
    needle->x = q;
    point_finder(needle);
  }
  else if (type == EDGE)
  {
    double eps = rms(3,x,0)*EPS;
    double *N1,N2[3];
    double q[3];
    
    tangent(pnt->point,N2);
    N1 = pnt->point->N;
    q[0] = x[0]+eps*N1[0]+EPS*N2[0];
    q[1] = x[1]+eps*N1[1]+EPS*N2[1];
    q[2] = x[2]+eps*N1[2]+EPS*N2[2];// q = pnt+eps*N1+N2
    
    needle->x = q;
    point_finder(needle);
  }
  else
    abortEr("No such type.\n");
  
  found = needle->ans;  
  nfp = needle->Nans;
  add_adjPnt(pnt,found,nfp);
  
  /* freeing */    
  free_needle(needle);
}

/* initializing and finding inner points*/
static void init_Points(Interface_T *interface,PointSet_T ***innP,PointSet_T ***edgP)
{
  PointSet_T **pnt_in,**pnt_ed;// inner and edge points
  int N_in, N_ed;// number of inner and edge points
  const int f = interface->fn;
  int *n = interface->patch->n;
  int i,j,k,im,iM,jm,jM,km,kM;// m for min and M for max
  int sum,ed,in;
  
  N_in = NumPoint(interface,INNER);
  N_ed = NumPoint(interface,EDGE);
  
  set_min_max_sum(n,f,&im,&iM,&jm,&jM,&km,&kM,&sum);
  assert(N_in+N_ed == sum);
  
  ed = in = 0;
  pnt_in = 0;
  pnt_ed = 0;
  FOR_ijk(i,j,k,im,iM,jm,jM,km,kM)
  {
    int l = L(n,i,j,k);
    int p = L2(n,f,i,j,k);
    
    /* excluding the points located on an internal interface */
    if (interface->point[p]->exterF == 0) continue;
      
    if (IsOnEdge(n,l))
    {
      alloc_PointSet(ed+1,&pnt_ed);
      pnt_ed[ed]->point = interface->point[p];
      ed++;
    }
    else
    {
      alloc_PointSet(in+1,&pnt_in);
      pnt_in[in]->point = interface->point[p];
      in++;
    }  
  }
  
  (*innP) = pnt_in;
  (*edgP) = pnt_ed;
}

/* number of inner points on an interface */
static int NumPoint(Interface_T *interface,enum Type type)
{
  int v;
  int *n = interface->patch->n;
  int f = interface->fn;
  
  if (type == INNER)
  {
    switch(f)
    {
      case I_0:
        v = (n[1]-2)*(n[2]-2);
        break;
      case I_n0:
        v = (n[1]-2)*(n[2]-2);
        break; 
      case J_0:
        v = (n[0]-2)*(n[2]-2);
        break; 
      case J_n1:
        v = (n[0]-2)*(n[2]-2);
        break; 
      case K_0:
        v = (n[0]-2)*(n[1]-2);
        break; 
      case K_n2:
        v = (n[0]-2)*(n[1]-2);
        break;
      default:
        abortEr("There is not such interface.\n");
    }
  }
  else if (type == EDGE)
  {
    switch(f)
    {
      case I_0:
        v = 2*(n[1]+n[2]-2);
        break;
      case I_n0:
        v = 2*(n[1]+n[2]-2);
        break; 
      case J_0:
        v = 2*(n[0]+n[2]-2);
        break; 
      case J_n1:
        v = 2*(n[0]+n[2]-2);
        break; 
      case K_0:
        v = 2*(n[0]+n[1]-2);
        break; 
      case K_n2:
        v = 2*(n[0]+n[1]-2);
        break;
      default:
        abortEr("There is not such interface.\n");
    }
  }
  else
    abortEr("There is no such type.\n");
  
  return v;
}

/* find inner boundary for Cartesian type */
static void FindInnerB_Cartesian_coord(Patch_T *patch)
{
  Interface_T **interface = patch->interface;
  int i,f;
  
  FOR_ALL(f,interface)
  {
    Point_T **point = interface[f]->point;
    FOR_ALL(i,point)
    {
      point[i]->innerB = 0;
    }
  }

}

/* find external faces for Cartesian type */
static void FindExterF_Cartesian_coord(Patch_T *patch)
{
  Interface_T **interface = patch->interface;
  int i,f;
  
  FOR_ALL(f,interface)
  {
    Point_T **point = interface[f]->point;
    FOR_ALL(i,point)
    {
      point[i]->exterF = 1;
    }
  }
  
}

/* filling point[?]->N */
static void fill_N(Patch_T *patch)
{
  Interface_T **interface = patch->interface;
  int i,f;
  
  FOR_ALL(f,interface)
  {
    Point_T **point = interface[f]->point;
    
    FOR_ALL(i,point)
    {
      normal_vec(point[i]);
    }
  }
}

/* 
// filling  basic elements:
// point[?]->ind, point[?]->face and point[?]->patch and
// also realize those nodes reach boundary of patch
*/
static void fill_basics(Patch_T *patch)
{
  int *n = patch->n;
  int f;
  
  for (f = I_0; f < TOT_FACE; f++)
  {
    Point_T **point;
    int i,j,k,im,iM,jm,jM,km,kM,sum;// m for min and M for max
    int t = 0;// test purposes
    
    set_min_max_sum(n,f,&im,&iM,&jm,&jM,&km,&kM,&sum);
    patch->interface[f]->point = alloc_point(sum);
    patch->interface[f]->patch = patch;
    patch->interface[f]->fn = f;
    patch->interface[f]->np = sum;
    point = patch->interface[f]->point;
    
    t = 0;
    FOR_ijk(i,j,k,im,iM,jm,jM,km,kM)
    {
      int p = L2(n,f,i,j,k);
      point[p]->ind = L(n,i,j,k);
      point[p]->face = f;
      point[p]->patch = patch;
      t++;
    }
    assert(t == sum);
  }
}

/*2d linear index format for interface points 
// L = k+n2*(j+n1*i) so depends on interface some of i,j,j,n2 and n1
// are zero or supremum as follows:
*/
static int L2(int *n,int f, int i, int j, int k)
{
  int v;
  
  switch(f)
  {
    case I_0:// i = 0
      v = k+n[2]*j;
      break;
    case I_n0:// i = n[0]-1
      v = k+n[2]*j;
      break;
    case J_0:// j = 0
      v = k+n[2]*i;
      break;
    case J_n1:// j = n[1]-1
      v = k+n[2]*i;
      break;
    case K_0:// k = 0
      v = j+n[1]*i;
      break;
    case K_n2:// k = n[2]-1
      v = j+n[1]*i;
      break;
    default:
      abortEr("Exceeding from total number of face.\n");
      break; 
  }
  
  return v;
}

/* setting min and max of total number of point on an interface.
// the purpose it serves is making loop over interface point easier.
// note: since in FOR_ijk the condition for loop is less than sign <
// then the maximum value must be the supremum; but the slice 
// should be at the exact value.
*/
static void set_min_max_sum(int *n,int f,int *im,int *iM,int *jm,int *jM,int *km,int *kM,int *sum)
{
  switch(f)
  {
    case I_0:
      *im = 0;
      *iM = 1;
      *jm = 0;
      *jM = n[1];
      *km = 0;
      *kM = n[2];
      *sum = n[1]*n[2];
      break;
    case I_n0:
      *im = n[0]-1;
      *iM = n[0];
      *jm = 0;
      *jM = n[1];
      *km = 0;
      *kM = n[2];
      *sum = n[1]*n[2];
      break;
    case J_0:
      *im = 0;
      *iM = n[0];
      *jm = 0;
      *jM = 1;
      *km = 0;
      *kM = n[2];
      *sum = n[0]*n[2];
      break;    
    case J_n1:
      *im = 0;
      *iM = n[0];
      *jm = n[1]-1;
      *jM = n[1];
      *km = 0;
      *kM = n[2];
      *sum = n[0]*n[2];
      break;    
    case K_0:
      *im = 0;
      *iM = n[0];
      *jm = 0;
      *jM = n[1];
      *km = 0;
      *kM = 1;
      *sum = n[0]*n[1]; 
      break;    
    case K_n2:
      *im = 0;
      *iM = n[0];
      *jm = 0;
      *jM = n[1];
      *km = n[2]-1;
      *kM = n[2];
      *sum = n[0]*n[1];
      break;
    default:
      abortEr("Exceeding from total number of face.\n");
      break;
  }
}

/* normal vector at interface of patch, points "OUTWARD" and "NORMALIZED";
// the normal vector is written in N at point structure;
// note: the members in point struct that must be filled before
// passing to this function are: "ind","patch","face".
// -> return value: a pointer to N.
*/
double *normal_vec(Point_T *point)
{
  if (strcmp_i(point->patch->coordsys,"Cartesian"))
  {
    normal_vec_Cartesian_coord(point);
  }
  //if (strcmp_i(point->patch->coordsys,"CubedSphere"))
  //{
    //normal_vec_CubedSphere(point);
  //}
  else
    abortEr_s("Normal for %s is not defined yet!\n",
      point->patch->coordsys);
    
  return point->N;
}

/* finding normal for Cartesian coord */
static void normal_vec_Cartesian_coord(Point_T *point)
{
  switch(point->face)
  {
    case I_0:
      point->N[0] = -1;
      point->N[1] = 0;
      point->N[2] = 0;
      break;
    case I_n0:
      point->N[0] = 1;
      point->N[1] = 0;
      point->N[2] = 0;
      break;
    case J_0:
      point->N[0] = 0;
      point->N[1] = -1;
      point->N[2] = 0;
      break;
    case J_n1:
      point->N[0] = 0;
      point->N[1] = 1;
      point->N[2] = 0;
      break;
    case K_0:
      point->N[0] = 0;
      point->N[1] = 0;
      point->N[2] = -1;
      break;
    case K_n2:
      point->N[0] = 0;
      point->N[1] = 0;
      point->N[2] = 1;
      break;
    default:
      abortEr("There is no such face.\n");
  }
}

/* reallocation mem for PointSet_T strcut with one extera block
// to put the last pointer to null.
*/
static void alloc_PointSet(int N,PointSet_T ***pnt)
{
  int i;
  
  assert(pnt);
  
  (*pnt) = realloc((*pnt),(N+1)*sizeof((*pnt)));
  pointerEr(*pnt);
  
  (*pnt)[N-1] = calloc(1,sizeof(*(*pnt)[i]));
  pointerEr((*pnt)[N-1]);
  
  (*pnt)[N] = 0;
}

/* free mem for PointSet_T strcut */
static void free_PointSet(PointSet_T **pnt)
{
  if(pnt) return;
  
  int i;
  FOR_ALL (i,pnt)
  {
    if (pnt[i]->adjPnt != 0)
    {
      free(pnt[i]->adjPnt);
    }
    
    free(pnt[i]);
  }
  
  free(pnt);
}

/* adding all of patches found in to point->adjPnt;
// and find if there is any node in adjacent patch correspond to
// this point pnt. moreover, study if this pnt located on some interfaces
// of adjacent patch, if yes put adjPnt->f[face_nun] = 1 and if not
// it is 0.
// note: it won't add any thing if the patch is already exists.
*/
static void add_adjPnt(PointSet_T *pnt,int *p, int np)
{
  int i,j;
  Flag_T flg = NONE;
  
  for (i = 0; i < np; i++)
  {
    for(j = 0; j < pnt->NadjPnt; j++)
    {
      if(p[i] == pnt->adjPnt->p)
        flg = FOUND;
    }
    
    if (flg == NONE)
    {
      /* realloc memory for new found point */
      pnt->adjPnt = 
        realloc(pnt->adjPnt,(pnt->NadjPnt+1)*sizeof(*pnt->adjPnt));
      pointerEr(pnt->adjPnt);
      
      pnt->adjPnt[pnt->NadjPnt].p = p[i];// fill p
      pnt->adjPnt[pnt->NadjPnt].FaceFlg = 0;// initialize
      find_adjNode(pnt,pnt->NadjPnt);// fill f and node and normals
      
      pnt->NadjPnt++;
     
    }
  }
}

/* study:
// 1. if this pnt collocated with some node of adjacent patch->fill node
// 2. if this pnt located on some interfaces of adjacent patch->fill f 
// and find its normal vector and its dot product with the pnt
// the last one gets used for realizing the geometry of neighbor of pnt
*/
static void find_adjNode(PointSet_T *pnt,const int N)
{
  Adjacent_T *adjPnt = &pnt->adjPnt[N];
  int ind = pnt->point->ind;
  double *x = pnt->point->patch->node[ind]->x;
  int *f = adjPnt->f;
  int *n = pnt->point->patch->n;
  Grid_T *grid = pnt->point->patch->grid;
  int ind2;
  
  ind2 = find_node(x,pnt->point->patch);
  if (ind2 >= 0)
  {
    adjPnt->node = ind2;
    if (IsOnFace(n,ind2,f))
    {
      Point_T po;
      double *N2;
      int i;
      po.ind = ind2;
      po.patch = grid->patch[adjPnt->p];
      adjPnt->FaceFlg = 1;
      
      for (i = 0; i < TOT_FACE; i++)
      {
        if (f[i])
        {
          po.face = i;
          N2 = normal_vec(&po);
          adjPnt->N2[i][0] = N2[0];
          adjPnt->N2[i][1] = N2[1];
          adjPnt->N2[i][2] = N2[2];
          adjPnt->N1dotN2[i] = dot(3,pnt->point->N,N2);
          if (EQL(adjPnt->N1dotN2[i],-1)) 
            adjPnt->FitF = i;
        }
      }
    }// if (IsOnFace(n,ind2,f))
  }// end of if (ind2 >= 0)
  
}

/* find the closest inner point to this point and then 
// put the vector made by the difference of these two points in N.
// this is the approximate tangent vector. note, the tangent
// vector won't be "NORMALIZED". I want to keep the magnetitue to be of
// the order of grid size.
*/
static void tangent(Point_T *pnt,double *N)
{
  int *const n = pnt->patch->n;
  const int ind = pnt->ind;
  const int f = pnt->face;
  double *x = pnt->patch->node[ind]->x;
  double *y;
  double s_ds = DBL_MAX;// smallest distance
  int s_in = -1;// index referring to point with smallest distance
  int i,j,k,im,iM,jm,jM,km,kM,sum;// m for min and M for max
  
  set_min_max_sum(n,f,&im,&iM,&jm,&jM,&km,&kM,&sum);
  
  FOR_ijk(i,j,k,im,iM,jm,jM,km,kM)
  {
    int l = L(n,i,j,k);
    double nrm;
    
    if (!IsOnEdge(n,l))
    {
      y = pnt->patch->node[l]->x;
      nrm = rms(3,x,y);
      
      if (LSS(nrm,s_ds)) s_in = l;
    }
  }
  
  if (s_in < 0)
    abortEr("Tangent vector could not be found.\n");
  
  /* note: it must be y-x to tilt toward interface not out of it */  
  y = pnt->patch->node[s_in]->x;
  N[0] = y[0]-x[0];
  N[1] = y[1]-x[1];
  N[2] = y[2]-x[2];
}
