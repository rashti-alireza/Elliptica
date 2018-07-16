/*
// Alireza Rashti
// June 2018
*/

#include "coordinate_systems.h"

/* making coordinates of nodes */
int fill_nodes(Grid_T *const grid)
{
  int i;
  FOR_ALL(i,grid->patch)
  {
    Patch_T *patch = grid->patch[i];
    
    /* if coord is Cartesian */
    if (strcmp_i(patch->coordsys,"Cartesian"))
    {
        make_coords_Cartesian_coord(patch);
    }
    
    else
      abortEr_s("There is no such %s coordinate.\n",patch->coordsys);
  }
  
  return EXIT_SUCCESS;
}

/* making value of coords. it is a general function for Cartesian type */
static void make_coords_Cartesian_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3];
  const unsigned U = countf(patch->node);
  unsigned i,j,k,l,*n;
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  initialize_collocation_struct(patch,&coll_s[2],2);
  n = patch->n;
  
  for (l = 0; l < U; l++)
  {
    double *x = patch->node[l]->x;
    
    IJK(l,n,&i,&j,&k);
    x[0] = point(i,&coll_s[0]);
    x[1] = point(j,&coll_s[1]);
    x[2] = point(k,&coll_s[2]);
    
    patch->node[l]->X = 0;
  }
}

/* initializing collocation struct 
// for making coords based on type of collocation.
//
// transformations:
// ========================
//
// note: a and b refer to the minimum and maximum of a line respectively,
// and n refers to number of points on that line.
//
// o. EquiSpaced:
//      dividing a line in equal sizes
//	-> grid space = (b-a)/(n-1)
//
// o. Chebyshev_extrema:
//	mapping a line [a,b] to [-1,1] and then using Chebyshev extrema
//	in [0,Pi] to find the nodes on the line. the following map is used:
//	x = 0.5*(a-b)*y + 0.5*(a+b) where x in [a,b] and y in [-1,1].
// 	y = cos(t).
// 	note that x and y are in reverse order, and the reason is that
// 	I wanted to have same increasing behavior between t and x.
//
//	
*/
static void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const coll_s,const unsigned dir)
{
  /* some assertions */
  assert(dir < 3);
  assert(patch->min[dir] < patch->max[dir]);
  
  coll_s->c = patch->collocation[dir];
  coll_s->n = patch->n[dir];
  coll_s->min = patch->min[dir];
  coll_s->max = patch->max[dir];
  assert(coll_s->n-1 != 0);

  if (coll_s->c == EquiSpaced)
  {
    coll_s->stp = (coll_s->max-coll_s->min)/(coll_s->n-1);
  }
  else if (coll_s->c == Chebyshev_Extrema)
  {
    coll_s->a = 0.5*(coll_s->min-coll_s->max);
    coll_s->b = 0.5*(coll_s->min+coll_s->max);
  }
  else
    abortEr("There is no such COLLOCATION.\n");
}

/* point value based on collocation */
static double point(const unsigned i, const struct Collocation_s *const coll_s)
{
  double v = DBL_MAX;
   
  if (coll_s->c == EquiSpaced)
  {
    v = coll_s->min + coll_s->stp*i;
  }
  else if (coll_s->c == Chebyshev_Extrema)
  {
    double phi = i*M_PI/(coll_s->n-1);
    
    v = coll_s->a*cos(phi)+coll_s->b;
  }
  else
    abortEr("There is no such collocation.\n");
  
  return v;
}

/* making collocation point based on patch properties and direction
// when x is normalized i.e x in [-1,1].
// Note: it allocates memory.
//-> return value: collocations point for given patch and direction */
double *make_normalized_collocation_1d(const Patch_T *const patch,const unsigned dir)
{
  const unsigned N = patch->n[dir];
  double *const x = alloc_double(N);
  unsigned i;
  
  if (patch->collocation[dir] == EquiSpaced)
    for (i = 0; i < N; i++)
      x[i] = -1+i*2.0/(N-1);
  else if (patch->collocation[dir] == Chebyshev_Extrema)
  {
    double phi0 = M_PI/(N-1);
    for (i = 0; i < N; i++)
      x[i] = cos(i*phi0);
  }
  else
    abortEr("No such collocation exists.\n");

  return x;
}
