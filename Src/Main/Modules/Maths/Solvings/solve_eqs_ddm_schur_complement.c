/*
// Alireza Rashti
// January 2019
*/

#include "solve_eqs_ddm_schur_complement.h"

/* 
   |B E||x|   |f|
   |F C||y| = |g|
   
   => Bx+Ey=f
      Fx+Cy=g
      
   => x = B^-1*(f-Ey)
   => (C-F*B^-1*E)y = g-F*B^-1*f
   
   Algorithm:
   
   1. solve BE'= E and Bf'= f
   2. compute g' = g - Ff'
   3. compute S  = C - FE'
   4. solve Sy = g'
   5. compute x = f'-E'y
   
*/

/* using Schur Complement domain decomposition method
// to solve equation. This method is capable of using direct solver
// like UMFPACK and also it is parallelizable.
*/
int ddm_schur_complement(Grid_T *const grid)
{
  /* picking up labeling, mapping etc. */
  preparing_ingredients(grid);
  
  return EXIT_SUCCESS;
}

/* since Shur Complement needs:
// 1. specific labeling for grid 
// 2. boundary points counted only once
// 3. adjacency patch of each patch and how other 
// patches use this patch to impose their boundary conditions.
// this function provides this ingredients
*/
static void preparing_ingredients(Grid_T *const grid)
{
  unsigned p;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    DDM_Schur_Complement_T *SchurC = calloc(1,sizeof(*SchurC));
    pointerEr(SchurC);
    SchurC->map = malloc(patch->nn*sizeof(*SchurC->map));
    pointerEr(SchurC->map);
    SchurC->inv = malloc(patch->nn*sizeof(*SchurC->inv));
    pointerEr(SchurC->inv);  
    
    /* making map and inv */
    make_map_and_inv(patch);
   
   patch->solution_man->method->Schur_Complement = 1;
   patch->solution_man->method->SchurC =  SchurC;
  }
}

/* making map and inv map for Schur complement method.
// inner mesh nodes come first, then outer boundary nodes 
// and then finally inner boundar nodes.
*/
static void make_map_and_inv(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const SchurC = 
    patch->solution_man->method->SchurC;
  unsigned *const map = SchurC->map;
  unsigned *const inv = SchurC->inv;
  const unsigned nn = patch->nn;
  const unsigned *n = patch->n;
  const unsigned nintfc = countf(patch->interface);
  unsigned i,j, intfc;
  
  /* keep tracking of points, 1 means counted, 0 means not */
  unsigned *flag_point = calloc(nn,sizeof(*flag_point));
  pointerEr(flag_point);
  
  /* filling inner points */
  j = 0;
  for (i = 0; i < nn; ++i)
  {
    if (!IsOnEdge(n,i))
    {
      map[i] = j;
      inv[j] = i;
      flag_point[i] = 1;
      j++;
    }
  }
  
  /* filling outer boundary points */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    unsigned nsfc = interface->ns;
    unsigned sfc;
    
    /* loop over all subfaces to filling outerbound */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      const unsigned *const B = subface->id;
      
      if (subface->outerB)/* if it reaches outer boundary */
      {
        for (i = 0; i < subface->np; ++i)
          if (!flag_point[B[i]])
          {
            map[B[i]] = j;
            inv[j] = B[i];
            flag_point[B[i]] = 1;
            j++;
          }
      }
    }
  }
  
  /* filling the other remaining points */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    unsigned nsfc = interface->ns;
    unsigned sfc;
    
    /* loop over all subfaces to filling outerbound */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      const unsigned *const T = subface->id;
      
      if (!subface->innerB && subface->exterF)/* not boundary and internal faces */
      {
        for (i = 0; i < subface->np; ++i)
          if (!flag_point[T[i]])
          {
            map[T[i]] = j;
            inv[j] = T[i];
            flag_point[T[i]] = 1;
            j++;
          }
      }
    }
  }
  
  if (j != nn)
  {
    abortEr("Not all points are mapped.\n");
  }
  
  free(flag_point);
}
