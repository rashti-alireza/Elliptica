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
// note: min and max refer to the minimum and maximum of a line respectively,
// and n refers to number of points on that line.
//
// o. EquiSpaced:
//      dividing a line in equal sizes
//	-> grid space = (max-min)/(n-1)
//
// o. Chebyshev_Extrema:
//	mapping a line [min,max] to [-1,1] and then using Chebyshev extrema
//	in [0,Pi] to find the nodes on the line. the following map is used:
// 	x = a*N+b = a*cos(t)+b i.e.
//	x = 0.5*(max-min)*N + 0.5*(max+min) where x in [min,max] and N in [-1,1].
// 	N = cos(t).
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
    coll_s->a = 0.5*(coll_s->max-coll_s->min);
    coll_s->b = 0.5*(coll_s->max+coll_s->min);
  }
  else
    abortEr("There is no such COLLOCATION.\n");
}

/* point value based on collocation 
// ->return value: collocation point on i-th location.
*/
static double point(const unsigned i, const struct Collocation_s *const coll_s)
{
  double x = DBL_MAX;
  
  /* x = min + i*grid_size */ 
  if (coll_s->c == EquiSpaced)
  {
    x = coll_s->min + coll_s->stp*i;
  }
  /* x =  a*N+b => x = a*cos(t)+b */
  else if (coll_s->c == Chebyshev_Extrema)
  {
    double t = i*M_PI/(coll_s->n-1);
    
    x = coll_s->a*cos(t)+coll_s->b;
  }
  else
    abortEr("There is no such collocation.\n");
  
  return x;
}

/* making collocation point based on patch properties and direction
// when x in [min,max].
// Note: it allocates memory.
//-> return value: collocations point for given patch and direction */
double *make_collocation_1d(const Patch_T *const patch,const unsigned dir,const double min,const double max)
{
  const unsigned N = patch->n[dir];
  double *const x = alloc_double(N);
  unsigned i;
  
  if (patch->collocation[dir] == EquiSpaced)
    for (i = 0; i < N; i++)
      x[i] = min+i*(max-min)/(N-1);
  else if (patch->collocation[dir] == Chebyshev_Extrema)
  {
    double t0 = M_PI/(N-1);
    for (i = 0; i < N; i++)
      x[i] = 0.5*(max-min)*cos(i*t0) +0.5*(max+min);
  }
  else
    abortEr("No such collocation exists.\n");

  return x;
}
