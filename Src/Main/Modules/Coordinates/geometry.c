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
  
  /* 
  // allocating interface struct,
  // allocating point struct and
  // counting total number of point in each interface np
  // filling ind, N, patch and face elements of point struct 
  */
  FOR_ALL(i,grid->patch)
  {
    Patch_T *const patch = grid->patch[i];
    alloc_interface(patch);
    fill_basics(patch);// filling basic elements
    fill_N(patch);// filling point[?]->N
  }
  
  /* find various geometry of each point */
  fill_geometry(grid);
  
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
  
  /* find neighbor properties */
  //FOR_ALL(i,grid->patch)
    //FindNeighbor(grid->patch[i]);
  
}

/* study neighbor of points to find quantities like 
// point->face, point->touch, point->adjPoint and others.
*/
/*static void FindNeighbor_Cartesian_grid(Patch_T *patch)
{
  PointSet_T **inner_p;
  
  const int U = patch->interface->np;
  int i;
  
  for (i = 0; i < U; i++)
  {
    Point_T *const P = patch->interface->point[i];
    Point_T **adjP = 0;
    adjP = find_adjPoint(P);// finding adjacent points of point P
    analyze_P(P,adjP);// analyzing P based on its adjacent points
    
    if (adjP != 0) free_2d(adjP);
  }
}*/

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
// point[?]->ind, point[?]->face and point[?]->patch 
*/
static void fill_basics(Patch_T *patch)
{
  Point_T **point;
  int *n = patch->n;
  int s;
  int i,j,k,p;
  
  s = n[0]*n[1];
  patch->interface[K_0]->point = alloc_point(s);
  patch->interface[K_0]->patch = patch;
  point = patch->interface[K_0]->point;
  p = 0;
  FOR_SURFACE(i,j,k,n[0],n[1],0) // k = 0 surface
  {
    point[p]->ind = L(n,i,j,k);
    point[p]->face = K_0;
    point[p]->patch = patch;
    p++;
  }
  patch->interface[K_0]->np = p;
  
  s = n[0]*n[1];
  patch->interface[K_n2]->point = alloc_point(s);
  patch->interface[K_n2]->patch = patch;
  point = patch->interface[K_n2]->point;
  p = 0;
  FOR_SURFACE(i,j,k,n[0],n[1],n[2]-1) // k = n[2]-1 surface
  {
    point[p]->ind = L(n,i,j,k);
    point[p]->face = K_n2;
    point[p]->patch = patch;
    p++;
  }
  patch->interface[K_n2]->np = p;
  
  s = n[0]*n[2];
  patch->interface[J_0]->point = alloc_point(s);
  patch->interface[J_0]->patch = patch;
  point = patch->interface[J_0]->point;
  p = 0;  
  FOR_SURFACE(i,k,j,n[0],n[2],0) // j = 0 surface
  {
    point[p]->ind = L(n,i,j,k);
    point[p]->face = J_0;
    point[p]->patch = patch;
    p++;
  }
  patch->interface[J_0]->np = p;
  
  s = n[0]*n[2];
  patch->interface[J_n1]->point = alloc_point(s);
  patch->interface[J_n1]->patch = patch;
  point = patch->interface[J_n1]->point;
  p = 0;
  FOR_SURFACE(i,k,j,n[0],n[2],n[1]-1) // j = n[1]-1 surface
  {
    point[p]->ind = L(n,i,j,k);
    point[p]->face = J_n1;
    point[p]->patch = patch;
    p++;
  }
  patch->interface[J_n1]->np = p;  
  
  s = n[1]*n[2];
  patch->interface[I_0]->point = alloc_point(s);
  patch->interface[I_0]->patch = patch;
  point = patch->interface[I_0]->point;
  p = 0;
  FOR_SURFACE(j,k,i,n[1],n[2],0) // i = 0 surface
  {
    point[p]->ind = L(n,i,j,k);
    point[p]->face = I_0;
    point[p]->patch = patch;
    p++;
  }
  patch->interface[I_0]->np = p;
  
  s = n[1]*n[2];
  patch->interface[I_n0]->point = alloc_point(s);
  patch->interface[I_n0]->patch = patch;
  point = patch->interface[I_n0]->point;
  p = 0;
  FOR_SURFACE(j,k,i,n[1],n[2],n[0]-1) // i = n[0]-1 surface
  {
    point[p]->ind = L(n,i,j,k);
    point[p]->face = I_n0;
    point[p]->patch = patch;
    p++;
  }
  patch->interface[I_n0]->np = p;
  
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
