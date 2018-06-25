/*
// Alireza Rashti
// June 2018
*/

#include "coordinate_systems.h"

/* making coordinates of nodes */
int fill_nodes(Grid_T *grid)
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
static void make_coords_Cartesian_coord(Patch_T *patch)
{
  struct Collocation_s coll_s[3];
  const int U = countf(patch->node);
  int i,j,k,l,*n;
  
  coll_s[0].f = coll_s[1].f = coll_s[2].f = patch->collocation;
  initialize_collocation_struct(patch,coll_s);
  n = patch->n;
  
  for (l = 0; l < U; l++)
  {
    double *x = patch->node[l]->cart;
    
    IJK(l,n,&i,&j,&k);
    x[0] = point(i,&coll_s[0]);
    x[1] = point(j,&coll_s[1]);
    x[2] = point(k,&coll_s[2]);
    
    patch->node[l]->curv = 0;
  }
}

/* initializing collocation struct 
// for making coords based on type of collocation.
*/
static void initialize_collocation_struct(Patch_T *patch,struct Collocation_s *coll_s)
{
  /* some assertions */
  assert(patch->min[0] < patch->max[0]);
  assert(patch->min[1] < patch->max[1]);
  assert(patch->min[2] < patch->max[2]);
  
  coll_s[0].n = patch->n[0];
  coll_s[0].min = patch->min[0];
  coll_s[0].max = patch->max[0];
  coll_s[1].n = patch->n[1];
  coll_s[1].min = patch->min[1];
  coll_s[1].max = patch->max[1];
  coll_s[2].n = patch->n[2];
  coll_s[2].min = patch->min[2];
  coll_s[2].max = patch->max[2];
  
  if (coll_s[0].f == EquiSpaced)
  {
    assert(coll_s[0].n-1 != 0);
    assert(coll_s[1].n-1 != 0);
    assert(coll_s[2].n-1 != 0);
    coll_s[0].stp = (coll_s[0].max-coll_s[0].min)/(coll_s[0].n-1);
    coll_s[1].stp = (coll_s[1].max-coll_s[1].min)/(coll_s[1].n-1);
    coll_s[2].stp = (coll_s[2].max-coll_s[2].min)/(coll_s[2].n-1);
  }
  else if (coll_s[0].f == Chebyshev_Zero)
  {
    coll_s[0].phi_in = 0.5*M_PI/coll_s[0].n;
    coll_s[0].phi_fi = M_PI-coll_s[0].phi_in;
    coll_s[0].a = 
      (coll_s[0].min-coll_s[0].max)/
        (cos(coll_s[0].phi_in)-cos(coll_s[0].phi_fi));
    coll_s[0].b = 
      (coll_s[0].max*cos(coll_s[0].phi_in)-
        coll_s[0].min*cos(coll_s[0].phi_fi))/
        (cos(coll_s[0].phi_in)-cos(coll_s[0].phi_fi));
    
    coll_s[1].phi_in = 0.5*M_PI/coll_s[1].n;
    coll_s[1].phi_fi = M_PI-coll_s[1].phi_in;
    coll_s[1].a = 
      (coll_s[1].min-coll_s[1].max)/
        (cos(coll_s[1].phi_in)-cos(coll_s[1].phi_fi));
        coll_s[1].b = 
    (coll_s[1].max*cos(coll_s[1].phi_in)-
        coll_s[1].min*cos(coll_s[1].phi_fi))/
        (cos(coll_s[1].phi_in)-cos(coll_s[1].phi_fi));

    coll_s[2].phi_in = 0.5*M_PI/coll_s[2].n;
    coll_s[2].phi_fi = M_PI-coll_s[2].phi_in;
    coll_s[2].a = 
      (coll_s[2].min-coll_s[2].max)/
        (cos(coll_s[2].phi_in)-cos(coll_s[2].phi_fi));
    coll_s[2].b = 
      (coll_s[2].max*cos(coll_s[2].phi_in)-
        coll_s[2].min*cos(coll_s[2].phi_fi))/
        (cos(coll_s[2].phi_in)-cos(coll_s[2].phi_fi));
  }
  else
    abortEr("There is no such COLLOCATION.\n");
}

/* point value based on collocation */
static double point(int i, struct Collocation_s *coll_s)
{
  double v;
   
  if (coll_s->f == EquiSpaced)
  {
    v = coll_s->min + coll_s->stp*i;
  }
  else if (coll_s->f == Chebyshev_Zero)
  {
    double phi = (i+0.5)*M_PI/coll_s->n;
    
    v = coll_s->a*cos(phi)+coll_s->b;
  }
  else
    abortEr("There is no such collocation.\n");
  
  return v;
}

