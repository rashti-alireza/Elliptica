/*
// Alireza Rashti
// June 2018
*/

#include "geometry.h"

#define EPS 1E-5

/* 
// realizing the geometry of grid such as how patches are glued
// normal vectors at boundary, boundary of grid and etc.
// basically, it goes through all of points on interfaces of a patch,
// finding their normals and adjacent points and then studying
// how best this point and its neighbor matched.
// ->return value: EXIT_SUCCESS
*/
int realize_geometry(Grid_T *const grid)
{
  unsigned i;
  
  FOR_ALL(i,grid->patch)
  {
    Patch_T *const patch = grid->patch[i];
    
    alloc_interface(patch);/* allocating interfaces */
    fill_basics(patch);/* filling basic elements */
    fill_N(patch);/* filling point[?]->N */
    flush_houseK(patch);/* set all of housK flag to zero */
  }
  
  /* find various geometry of each point */
  fill_geometry(grid);
  
  /* check if all point have been found and then flush them */
  FOR_ALL(i,grid->patch)
  {
    check_houseK(grid->patch[i]);
    flush_houseK(grid->patch[i]);
  }
  
  /* printing boundary for test purposes */
  /*pr_boundary(grid);*/
  
  return EXIT_SUCCESS;
}

/* filling the geometry of point struct */
static void fill_geometry(Grid_T *const grid)
{
  sFunc_PtoV_T **func;
  unsigned i;
  
  init_func_PtoV(&func);/* initialize struct */
  /* adding func to struct to be called. each coord must have 
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
    run_func_PtoV(func,"FindExterF",patch);/* find external faces */
    run_func_PtoV(func,"FindInnerB",patch);/* find inner boundary */
  }
  free_2d(func);/* freeing func struct */
  
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
static int realize_neighbor(Patch_T *const patch)
{
  const unsigned nf = countf(patch->interface);
  unsigned f;
  
  for (f = 0; f < nf; f++)
  {
    Interface_T *interface = patch->interface[f];
    PointSet_T **innerP, **edgeP;
    init_Points(interface,&innerP,&edgeP);/* initializing inner and 
                                          // edge points */
    realize_adj(innerP);/* realizing the adjacency INNER one */
    realize_adj(edgeP);/* realizing the adjacency EDGE one */
    
    free_PointSet(innerP);
    free_PointSet(edgeP);
  }
  return EXIT_SUCCESS;
}

/* realizing adjacency of given points
// ->return value: EXIT_SUCCESS.
*/
static int realize_adj(PointSet_T **const Pnt)
{
  unsigned p;
  
  FOR_ALL(p,Pnt)
  {
    if (Pnt[p]->Pnt->houseK == 1) continue;
    
    find_adjPnt(Pnt[p]);
    analyze_adjPnt(Pnt[p]);
  }
  
  return EXIT_SUCCESS;
}

/* analyzing the adjacent point of point pnt; 
// see which one best describe the adjacent boundary of the interface
// and then fill the flags of Point_T.
*/
static void analyze_adjPnt(PointSet_T *const Pnt)
{
  Point_T *const p1 = Pnt->Pnt;
  Point_T *p2;
  AdjPoint_T *adjp1;
  
  /* this point is on outer boundary */
  if (IsOutBndry(Pnt))
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
  /* which point best meets the normal conditon (N2.N1 = -1)  */
  else if (IsNormalFit(Pnt))
  {
    adjp1 = &Pnt->adjPnt[Pnt->idFit];
    
    p1->touch = 1;
    p1->copy  = 1;
    p1->adjPatch = adjp1->p;
    p1->adjFace  = adjp1->CopyFace;
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
  /* cases in which N1.N2 = 0 and reach outerbound */
  else if (IsOrthOutBndry(Pnt))
  {
    p1->outerB = 1;
  }
  /* cases in which although N1.N2 != -1 but copy still possible */
  else if (IsMildOrth(Pnt))
  {
    adjp1 = &Pnt->adjPnt[Pnt->idFit];
    
    p1->touch = 1;
    p1->copy  = 1;
    p1->adjPatch = adjp1->p;
    p1->adjFace  = adjp1->CopyFace;
    p2 = get_p2(Pnt);
    assert(p2);
    p1->adjPoint = p2;
  }
  /* cases in which the points needs interpolation */
  else if (IsInterpolation(Pnt))
  {
    if (Pnt->overlap == 0) p1->touch = 1;
    else  		   p1->touch = 0;
    
    p1->copy = 0;
    p1->adjPatch = Pnt->adjPnt[Pnt->idInterp].p;
    
    if (Pnt->adjPnt[Pnt->idInterp].FaceFlg == 1)
    {
      p1->IntFace = 1;
      p1->adjFace = Pnt->adjPnt[Pnt->idInterp].InterpFace;
    }
  }
  else
  {
    unsigned ind = p1->ind;
    double *x = p1->patch->node[ind]->x;
    
    fprintf(stderr,"This point(%f,%f,%f) has not been found.\n",
                    x[0],x[1],x[2]);
    abortEr("Incomplete function.\n");
  }
  p1->houseK = 1;
}

/* does this need interpolation?
// algorithm:
// 1. if: adjPnt is not on an interface safly pick this patch as interpolation
// 2. else: pick all of adjPnts which are on an interface, chose the closest point
// on the same interface and get its normal, if dot product of this normal
// with Pnt normal is negative, save it and finally pick the one which
// has this dot product closest to -1.
// ->return value: 1 if yes; 0 otherwise.
*/
static int IsInterpolation(PointSet_T *const Pnt)
{
  struct Interp_S
  {
    unsigned adjid;/* adjPnt id */
    unsigned fid;/* adjPnt face id */
    double dot;/* N1 dot N2 */
  }*interp = 0;
  double tmp;
  unsigned i,j,id;
  
  /* check if Pnt reach inside any adj patch */
  for (i = 0; i < Pnt->NadjPnt; i++)
  {
    if (Pnt->adjPnt[i].FaceFlg == 0)
    {
      Pnt->idInterp = i;
      Pnt->overlap = 1;
      return 1;
    }
  }
  
  /* find the all related N1dotN2 */
  j = 0;
  for (i = 0; i < Pnt->NadjPnt; i++)
  {
    if (Pnt->adjPnt[i].FaceFlg == 1)
    {
      for (j = 0; j < TOT_FACE; j++)
      {
        struct Face_S *const f = &Pnt->adjPnt[i].fs[j];
        if (f->on_f == 1 && LSS(f->N1dotN2,0))
        {
          interp = realloc(interp,(j+1)*sizeof(*interp));
          interp[j].adjid = i;
          interp[j].fid   = j;
          interp[j].dot   = f->N1dotN2;
          j++;
        }
      }/* for (j = 0; j < TOT_FACE; j++) */
    }
  }/* end of for (i = 0; i < Pnt->NadjPnt; i++) */
  
  /* find the closest N1dotN2 to -1 */
  if (j > 0)
  {
    tmp = interp[0].dot;
    Flag_T flg = NONE;
    
    for (i = 1; i < j; i++)
      if (LSS(interp[i].dot,tmp))
      {
        id = i;
        flg = FOUND;
      }
      
    if (flg == FOUND)
    {
      Pnt->idInterp = id;
      Pnt->adjPnt[id].InterpFace = interp[id].fid;
      Pnt->overlap = 0;
    }
    else
      abortEr("It seems the algorithm is wrong!\n");
    
    free(interp);
    return 1;
  }
  
  return 0;
}

/* if pnt and adjPnt have N2 has tangentical component to N1
// less than grid size of adjPnt and N1 can goes to adjacent patch, 
// this boundary condition is still valid.
// ->return value: 1 if yes; otherwise 0.
*/
static int IsMildOrth(PointSet_T *const Pnt)
{
  const unsigned node = Pnt->Pnt->ind;
  Patch_T *const patch = Pnt->Pnt->patch;
  double *const x = patch->node[node]->x;
  const double eps = rms(3,x,0)*EPS;
  double *const N = Pnt->Pnt->N;
  unsigned i,f;
  Needle_T *needle = alloc_needle();
  double q[3] = {x[0]+eps*N[0],
                  x[1]+eps*N[1],
                   x[2]+eps*N[2]};/* q = pnt+eps*N */
  needle->grid = patch->grid;
  needle_ex(needle,patch);
  needle->x = q;
  
  for (i = 0; i < Pnt->NadjPnt; i++)
  {
    for (f = 0; f < TOT_FACE; f++)
    {
      double s = Pnt->adjPnt[i].fs[f].N1dotN2;
      /* make sure the product of normal is negative */
      if (LSS(s,0))
      {
        point_finder(needle);
        unsigned nfp = needle->Nans;
        /* make sure q is on adjPatch */
        if (nfp > 0)
        {
          double grid_size = find_grid_size(Pnt,i,f);
          
          if (LSS(sqrt(1-SQR(s)),grid_size))
          {
            Pnt->idOrth = i;
            return 1;
          }
        }
      }/* end of if (LSS(s,0)) */
    }/* end of for (f = 0; f < TOT_FACE; f++) */
  }/* end of for (i = 0; i < Pnt->NadjPnt; i++) */
  
  free_needle(needle);
  
  return 0;
}

/* finding the closest point to adjPnt in the same patch
// with index i in PointSet_T then finding its distance to adjPnt.
// ->return value: distance between adjPnt and its closest point.
*/
static double find_grid_size(const PointSet_T *const Pnt,const unsigned i,const unsigned f)
{
  Grid_T *const grid = Pnt->Pnt->patch->grid;
  const unsigned p = Pnt->adjPnt[i].p;
  Interface_T *const face = grid->patch[p]->interface[f];
  double s = DBL_MAX;
  double *x = grid->patch[p]->node[Pnt->adjPnt[i].node]->x;
  double *y;
  unsigned ind,j;
  
  FOR_ALL(j,face->point)
  {
    double r;
    ind = face->point[j]->ind;
    y = grid->patch[p]->node[ind]->x;
    r = rms(3,x,y);

    if (LSS(r,s))
    {
      ind = j;
      s = r;
    }
  }
  
  return s;
}

/* is this and outerboundary point?
// if no points are in adjPnt yes.
// ->return value = 1 if outerbound; 0 otherwise.
*/
static int IsOutBndry(PointSet_T *const Pnt)
{
  if (Pnt->NadjPnt == 0)
    return 1;
    
  return 0;
}

/* there are some cases -take place generally when interfaces of
// two patches reach outerbound - in which two normals are orthogonal 
// (or closely orthogonal) and those points can be flaged as outerbound.
// the algorithm goes like that:
// if normal orthogonal tilt each point's normal toward the other one
// if the new point couldn't be found in any other patches, flaged this
// as outerbound.
// ->return value: if their normals are orthogonal and surrounded by
// outerbound 1, otherwise 0.
*/
static int IsOrthOutBndry(PointSet_T *const Pnt)
{
  unsigned i,j;
  
  for (i = 0; i < Pnt->NadjPnt; i++)
    for (j = 0; j < TOT_FACE; j++)
      if (Pnt->adjPnt[i].fs[j].OrthFlg == 1) 
        if (ReachBnd(Pnt,i,j))
        {
          Pnt->idOrth = i;
          return 1;
        }
    
  return 0;
}

/* check if one tilts normal toward the adjPnt normal 
// there is no other patches, which means outerboundary.
// ->return value: 1 if outerboundary; 0 otherwise.
*/
static int ReachBnd(PointSet_T *const Pnt,const unsigned p,const unsigned f)
{
  unsigned node = Pnt->Pnt->ind;
  Patch_T *const patch = Pnt->Pnt->patch;
  double *x = patch->node[node]->x;
  double *N1 = Pnt->Pnt->N;
  double *N2 = Pnt->adjPnt[p].fs[f].N2;
  double eps = rms(3,x,0)*EPS;
  double q[3] = {x[0]+eps*N1[0]+EPS*N2[0],
                 x[1]+eps*N1[1]+EPS*N2[1],
                 x[2]+eps*N1[2]+EPS*N2[2]};
  Needle_T *needle = alloc_needle();
  needle->grid = patch->grid;
  needle_ex(needle,Pnt->Pnt->patch);
  needle->x = q;
  point_finder(needle);
  
  if (needle->Nans == 0)
    return 1;
  
  free_needle(needle);
  
  return 0;
}

/* find the point p on the interface which has the fittest normal.
// ->return value: a pointer to point struct on interface
*/
static Point_T *get_p2(const PointSet_T *const Pnt)
{
  unsigned f = Pnt->adjPnt[Pnt->idFit].CopyFace;
  unsigned p = Pnt->adjPnt[Pnt->idFit].p;
  const unsigned node = Pnt->adjPnt[Pnt->idFit].node;
  Patch_T *const patch = Pnt->Pnt->patch->grid->patch[p];
  Interface_T *const face = patch->interface[f];
  unsigned i;
  
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
static int IsNormalFit(PointSet_T *const Pnt)
{
  unsigned i,f;
  
  for (i = 0; i < Pnt->NadjPnt; i++)
    for (f = 0; f < TOT_FACE; f++)
    if (Pnt->adjPnt[i].fs[f].FitFlg == 1) 
    {
      Pnt->idFit = i;
      Pnt->adjPnt[i].CopyFace = f;
      return 1;
    }
  
  return 0;
}
/* if none of the adjacent points are on an interface,
// thus, this point overlaps with other patch(es) 
// ->retrun value: 1 if it is overlap, 0 otherwise
*/
static int IsOverlap(PointSet_T *const Pnt)
{
  unsigned i;
  
  for (i = 0; i < Pnt->NadjPnt; i++)
    if (Pnt->adjPnt[i].FaceFlg == 1) return 0;
  
  return 1;
}

/* finding the adjacent points of point pnt */
static void find_adjPnt(PointSet_T *const pnt)
{
  unsigned ind = pnt->Pnt->ind;
  unsigned *found;/* patches have been found */
  unsigned nfp;/* number of found patches */
  double *x = pnt->Pnt->patch->node[ind]->x;
  Needle_T *needle = alloc_needle();
  needle->grid = pnt->Pnt->patch->grid;
  needle_ex(needle,pnt->Pnt->patch);
  needle->x = x;
  point_finder(needle); 
  found = needle->ans;  
  nfp = needle->Nans;
  add_adjPnt(pnt,found,nfp);
  
  /* freeing */    
  free_needle(needle);
}

/* initializing and finding inner points*/
static void init_Points(const Interface_T *const interface,PointSet_T ***const innP,PointSet_T ***const edgP)
{
  PointSet_T **pnt_in,**pnt_ed;/* inner and edge points */
  unsigned N_in, N_ed;/* number of inner and edge points */
  const unsigned f = interface->fn;
  const unsigned *const n = interface->patch->n;
  unsigned i,j,k,im,iM,jm,jM,km,kM;/* m for min and M for max */
  unsigned sum,ed,in;
  
  N_in = NumPoint(interface,INNER);
  N_ed = NumPoint(interface,EDGE);
  
  set_min_max_sum(n,f,&im,&iM,&jm,&jM,&km,&kM,&sum);
  assert(N_in+N_ed == sum);
  
  ed = in = 0;
  pnt_in = 0;
  pnt_ed = 0;
  FOR_ijk(i,j,k,im,iM,jm,jM,km,kM)
  {
    unsigned l = L(n,i,j,k);
    unsigned p = L2(n,f,i,j,k);
    
    /* excluding the points located on an internal interface */
    if (interface->point[p]->exterF == 0) 
    {
      interface->point[p]->houseK = 1;
      continue;
    }
      
    if (IsOnEdge(n,l))
    {
      alloc_PointSet(ed+1,&pnt_ed);
      pnt_ed[ed]->Pnt = interface->point[p];
      ed++;
    }
    else
    {
      alloc_PointSet(in+1,&pnt_in);
      pnt_in[in]->Pnt = interface->point[p];
      in++;
    }  
  }
  
  (*innP) = pnt_in;
  (*edgP) = pnt_ed;
}

/* number of inner points on an interface
// ->return value: number of points on the given interface.
*/
static unsigned NumPoint(const Interface_T *const interface,const enum Type type)
{
  unsigned v = UINT_MAX;
  unsigned *n = interface->patch->n;
  unsigned f = interface->fn;
  
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
static void FindInnerB_Cartesian_coord(Patch_T *const patch)
{
  Interface_T **interface = patch->interface;
  unsigned i,f;
  
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
static void FindExterF_Cartesian_coord(Patch_T *const patch)
{
  Interface_T **interface = patch->interface;
  unsigned i,f;
  
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
static void fill_N(Patch_T *const patch)
{
  Interface_T **interface = patch->interface;
  unsigned i,f;
  
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
static void fill_basics(Patch_T *const patch)
{
  unsigned *n = patch->n;
  unsigned f;
  
  for (f = I_0; f < TOT_FACE; f++)
  {
    Point_T **point;
    unsigned i,j,k,im,iM,jm,jM,km,kM,sum;/* m for min and M for max */
    unsigned t = 0;/* test purposes */
    
    set_min_max_sum(n,f,&im,&iM,&jm,&jM,&km,&kM,&sum);
    patch->interface[f]->point = alloc_point(sum);
    patch->interface[f]->patch = patch;
    patch->interface[f]->fn = f;
    patch->interface[f]->np = sum;
    point = patch->interface[f]->point;
    
    t = 0;
    FOR_ijk(i,j,k,im,iM,jm,jM,km,kM)
    {
      unsigned p = L2(n,f,i,j,k);
      point[p]->ind = L(n,i,j,k);
      point[p]->face = f;
      point[p]->patch = patch;
      t++;
    }
    assert(t == sum);
  }
}

/* 2d linear index format for interface points.
// L = k+n2*(j+n1*i) so depends on interface some of i,j,j,n2 and n1
// are zero or supremum as follows:
// ->return value: 2d index.
*/
static unsigned L2(const unsigned *const n,const unsigned f, const unsigned i, const unsigned j, const unsigned k)
{
  unsigned v = UINT_MAX;
  
  switch(f)
  {
    case I_0:/* i = 0 */
      v = k+n[2]*j;
      break;
    case I_n0:/* i = n[0]-1 */
      v = k+n[2]*j;
      break;
    case J_0:/* j = 0 */
      v = k+n[2]*i;
      break;
    case J_n1:/* j = n[1]-1 */
      v = k+n[2]*i;
      break;
    case K_0:/* k = 0 */
      v = j+n[1]*i;
      break;
    case K_n2:/* k = n[2]-1 */
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
static void set_min_max_sum(const unsigned *const n,const unsigned f,unsigned *const im,unsigned *const iM,unsigned *const jm,unsigned *const jM,unsigned *const km,unsigned *const kM,unsigned *const sum)
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
double *normal_vec(Point_T *const point)
{
  if (strcmp_i(point->patch->coordsys,"Cartesian"))
  {
    normal_vec_Cartesian_coord(point);
  }
  /*if (strcmp_i(point->patch->coordsys,"CubedSphere"))
  //{
    //normal_vec_CubedSphere(point);
  //}*/
  else
    abortEr_s("Normal for %s is not defined yet!\n",
      point->patch->coordsys);
    
  return point->N;
}

/* finding normal for Cartesian coord */
static void normal_vec_Cartesian_coord(Point_T *const point)
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
static void alloc_PointSet(const unsigned N,PointSet_T ***const pnt)
{
  
  assert(pnt);
  
  (*pnt) = realloc((*pnt),(N+1)*sizeof((*pnt)));
  pointerEr(*pnt);
  
  (*pnt)[N-1] = calloc(1,sizeof(*(*pnt)[N-1]));
  pointerEr((*pnt)[N-1]);
  
  (*pnt)[N] = 0;
}

/* free mem for PointSet_T strcut */
static void free_PointSet(PointSet_T **const pnt)
{
  if(pnt) return;
  
  unsigned i;
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
// of adjacent patch, if yes put adjPnt->fs[face_nun].on = 1 and if not
// it is 0.
// note: it won't add any thing if the patch is already exists.
*/
static void add_adjPnt(PointSet_T *const pnt,const unsigned *const p, const unsigned np)
{
  unsigned i,j;
  
  for (i = 0; i < np; i++)
  {
    for(j = 0; j < pnt->NadjPnt; j++)
      if(p[i] == pnt->adjPnt[j].p) return;
        
    /* realloc memory for new found point */
    pnt->adjPnt = 
      realloc(pnt->adjPnt,(pnt->NadjPnt+1)*sizeof(*pnt->adjPnt));
    pointerEr(pnt->adjPnt);
    
    /* initializing */
    AdjPoint_T *const adjPnt = &pnt->adjPnt[pnt->NadjPnt];
    adjPnt->p = p[i];/* fill p */
    adjPnt->FaceFlg = 0;
    adjPnt->on_c = 0;
    
    fill_adjPnt(pnt,pnt->NadjPnt);/* filling elemetns of adjPnt struct */
    
    pnt->NadjPnt++;
  }
}

/* study:
// 1. if this pnt collocated with some node of adjacent patch
// 2. if this pnt located on some interfaces of adjacent patch
// 3. find adjPnt normal vector and its dot product with the pnt
// the last one gets used for realizing the geometry of neighbor of pnt
*/
static void fill_adjPnt(PointSet_T *const pnt,const unsigned N)
{
  AdjPoint_T *const adjPnt = &pnt->adjPnt[N];
  unsigned ind = pnt->Pnt->ind;
  const double *const x = pnt->Pnt->patch->node[ind]->x;
  unsigned f[TOT_FACE];
  const Grid_T *const grid = pnt->Pnt->patch->grid;
  Flag_T flg;
  unsigned ind2;
  
  /* find if it is on a node (collocation) */
  ind2 = find_node(x,grid->patch[adjPnt->p],&flg);
  if (flg == FOUND)
  {
    adjPnt->node = ind2;
    adjPnt->on_c = 1;
  }
  
  /* find if it is on a face or not */  
  if (IsOnFace(x,grid->patch[adjPnt->p],f))
  {
    Point_T po;
    double *N2;
    unsigned i;
    
    po.patch = grid->patch[adjPnt->p];
    adjPnt->FaceFlg = 1;
    
    for (i = 0; i < TOT_FACE; i++)
    {
      if (f[i])
      {
        po.ind = node_onFace(x,i,po.patch);
        po.face = i;
        N2 = normal_vec(&po);
        adjPnt->fs[i].on_f = 1;
        adjPnt->fs[i].N2[0] = N2[0];
        adjPnt->fs[i].N2[1] = N2[1];
        adjPnt->fs[i].N2[2] = N2[2];
        adjPnt->fs[i].N1dotN2 = dot(3,pnt->Pnt->N,N2);
        
        if (EQL(adjPnt->fs[i].N1dotN2,-1)) 
        {
          adjPnt->fs[i].FitFlg  = 1;
          adjPnt->fs[i].OrthFlg = 0;
        }
        else if (EQL(adjPnt->fs[i].N1dotN2,0)) 
        {
          adjPnt->fs[i].FitFlg  = 0;
          adjPnt->fs[i].OrthFlg = 1;
        }
      }
      else
        adjPnt->fs[i].on_f = 0;
    }
  }/* if (IsOnFace(n,ind2,f)) */
  
}

/* find the closest inner point to this point and then 
// put the vector made by the difference of these two points in N.
// this is the approximate tangent vector. note, the tangent
// vector won't be "NORMALIZED". I want to keep the magnetitue to be of
// the order of grid size.
*/
static void tangent(const Point_T *const pnt,double *const N)
{
  const unsigned *const n = pnt->patch->n;
  const unsigned ind = pnt->ind;
  const unsigned f = pnt->face;
  double *x = pnt->patch->node[ind]->x;
  double *y;
  double s_ds = DBL_MAX;/* smallest distance */
  unsigned s_in;/* index referring to point with smallest distance */
  unsigned i,j,k,im,iM,jm,jM,km,kM,sum;/* m for min and M for max */
  Flag_T flg = NONE;
  
  set_min_max_sum(n,f,&im,&iM,&jm,&jM,&km,&kM,&sum);
  
  FOR_ijk(i,j,k,im,iM,jm,jM,km,kM)
  {
    unsigned l = L(n,i,j,k);
    double nrm;
    
    if (!IsOnEdge(n,l))
    {
      y = pnt->patch->node[l]->x;
      nrm = rms(3,x,y);
      
      if (LSS(nrm,s_ds)) 
      {
        flg = FOUND;
        s_in = l;
      }
    }
  }
  
  if (flg == NONE)
    abortEr("Tangent vector could not be found.\n");
  
  /* note: it must be y-x to tilt toward interface not out of it */  
  y = pnt->patch->node[s_in]->x;
  N[0] = y[0]-x[0];
  N[1] = y[1]-x[1];
  N[2] = y[2]-x[2];
}
