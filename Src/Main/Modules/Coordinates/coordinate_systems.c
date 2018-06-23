/*
// Alireza Rashti
// June 2018
*/

#include "coordinate_systems.h"

/* making coordinates of nodes */
int make_coordinates(Grid_T *grid)
{
  int i;
  for_all_patches_macro(i,grid)
  {
    Patch_T *patch = grid->patch[i];
    int U = countf(patch->node);
    
    /* if coord is Cartesian */
    if (!strcmp(patch->coordsys,"Cartesian"))
    {
    
      if(!strcmp(patch->collocation,"EquiSpaced"))
        make_coords_Cartesian(patch,U,EQUISPACED);
    
      //else if(!strcmp(patch->collocation,"Chebyshev_Zero"))
        //make_coords_Cartesian(patch,U,CHEBYSHEV_ZERO);
      else
        abortEr("There is no such collocation.\n");
    }
    
    /* if coord is CubedSphere */
    /*
    else if (!strcmp(patch->coord,"CubedSphere"))
    {
      if(!strcmp(patch->collocation,"Chebyshev_Zero"))
        make_coords_CubedSphere(patch,U,CHEBYSHEV_ZERO);
      else
        abortEr("There is no such collocation.\n");
    }*/
    
  }
  
  return EXIT_SUCCESS;
}

/* making value of coords. it is a general function for Cartesian type */
static void make_coords_Cartesian(Patch_T *patch,const int U, enum Flow coll_e)
{
  struct Collocation coll_s[3];
  int i,j,k,l,*n;
  
  coll_s[0].f = coll_s[1].f = coll_s[2].f = coll_e;
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
static void initialize_collocation_struct(Patch_T *patch,struct Collocation *coll_s)
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
  
  if (coll_s[0].f == EQUISPACED)
  {
    coll_s[0].stp = (coll_s[0].max-coll_s[0].min)/(coll_s[0].n-1);
    coll_s[1].stp = (coll_s[1].max-coll_s[1].min)/(coll_s[1].n-1);
    coll_s[2].stp = (coll_s[2].max-coll_s[2].min)/(coll_s[2].n-1);
  }
  else if (coll_s[0].f == CHEBYSHEV_ZERO)
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
static double point(int i, struct Collocation *coll_s)
{
  double v;
   
  if (coll_s->f == EQUISPACED)
  {
    v = coll_s->min + coll_s->stp*i;
  }
  else if (coll_s->f == CHEBYSHEV_ZERO)
  {
    double phi = (i+0.5)*M_PI/coll_s->n;
    
    v = coll_s->a*cos(phi)+coll_s->b;
  }
  else
    abortEr("There is no such collocation.\n");
  
  return v;
}

/* if it needs to fill up double *curv 
  if (!strcmp(patch->coord,"Cartesian")
  {
    allocating mem for double *curv
    for (i = 0; i < U; i++)
    {
      patch->node[i]->curv = malloc(3*sizeof(*patch->node[i]->curv));
      pointerEr(patch->node[i]->curv);
    }
  }
*/