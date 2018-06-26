/*
// Alireza Rashti
// June 2018
*/

#include "geometry.h"

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
    alloc_interface(patch);
    fill_basics(patch);// filling basic elements
    fill_N(patch);// filling point[?]->N
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
  /* adding func to struc to be called each coord must have 
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
    RealizeNeighbor(grid->patch[i]);
  
}

/* study neighbor of points to find quantities like 
// point->face, point->touch, point->adjPoint and others.
// algorithm:
// for each interface do
// 	a. seprate each interface to its inner and edge points
// 	b. realize the adjacency of inner points
// 	c. using 2, realize the adjacency of edge points
*/
static int RealizeNeighbor(Patch_T *patch)
{
  PointSet_T **inner_p;
  const int nf = countf(patch->interface);
  int f;
  
  for (f = 0; f < nf; f++)
  {
    Interface_T *interface = patch->interface[f];
    PointSet_T **innerP, **edgeP;
    
    init_Points(interface,&innerP,&edgeP);
    //realize_adj(innerP);
    //realize_adj(edgeP);
    
    free_PointSet(innerP);
    free_PointSet(edgeP);
  }
  return EXIT_SUCCESS;
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
  pnt_in = alloc_PointSet(N_in);
  N_ed = NumPoint(interface,EDGE);
  pnt_ed = alloc_PointSet(N_ed);
  
  set_min_max_sum(n,f,&im,&iM,&jm,&jM,&km,&kM,&sum);
  
  ed = in = 0;
  FOR_ijk(i,j,k,im,iM,jm,jM,km,kM)
  {
    int l = L(n,i,j,k);
    int p = L2(n,f,i,j,k);
    
    if (IsThisEdge(n,l))
    {
      pnt_ed[ed]->point = interface->point[p];
      pnt_ed[ed]->type = EDGE;
      ed++;
    }
    else
    {
      pnt_in[in]->point = interface->point[p];
      pnt_in[in]->type = INNER;
      in++;
    }  
  }
  assert(in+ed == sum);
  
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
      patch->node[L(n,i,j,k)]->Bpoint = point[p];
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

/* normal vector at interface of patch, it points outward and normalized;
// the normal vector is written in N at point structure;
// moreover, it returns a pointer to this N as well.
// note: the members in point struct that must be filled before
// passing to this function are: ind,patch,face.
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

/* allocation mem for PointSet_T strcut */
static void *alloc_PointSet(int N)
{
  PointSet_T **pnt;
  int i;
  
  pnt = calloc(N+1,sizeof(*pnt));
  pointerEr(pnt);
  
  for (i = 0; i < N; i++)
  {
    pnt[i] =  calloc(1,sizeof(*pnt[i]));
    pointerEr(pnt);
  }
  
  return pnt;
}

/* free mem for PointSet_T strcut */
static void free_PointSet(PointSet_T **pnt)
{
  if(pnt) return;
  
  int i;
  FOR_ALL (i,pnt)
  {
    if (pnt[i]->adj != 0)
      free(pnt[i]->adj);
    
    free(pnt[i]);
  }
  
  free(pnt);
}