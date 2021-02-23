/*
// Alireza Rashti
// June 2018
*/

#include "cartesian_coordinate.h"

/* making value of coords. it is a general function for Cartesian type */
void make_nodes_Cartesian_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const Uint U = patch->nn;
  const Uint *const n = patch->n;
  Uint i,j,k,l;
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  initialize_collocation_struct(patch,&coll_s[2],2);
  
  for (l = 0; l < U; l++)
  {
    double *x = patch->node[l]->x;
    
    ijk_to_i_j_k(l,n,&i,&j,&k);
    x[0] = point_value(i,&coll_s[0]);
    x[1] = point_value(j,&coll_s[1]);
    x[2] = point_value(k,&coll_s[2]);
    
    /* since X and x are the same we have: */
    patch->node[l]->X = x;
  }
}

/* making Jacobian transformation for Cartesian coord. */
void make_JacobianT_Cartesian_coord(Patch_T *const patch)
{
  patch->JacobianT->j      = JT_Cartesian_patch;
}

/* Jacobian transformation for Cartesian patch.
// ->return value: dq2/dq1
*/
double JT_Cartesian_patch(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  double j = 0;
  
  if (q2_e%3 == q1_e%3)/* e.g. _y_%3 = 1 and _b_%3 = 1 */
    j = 1;
    
  return j;
  
  UNUSED(patch);
  UNUSED(p);
}

/*filling patch struct for Cartesian*/
void fill_patches_Cartesian_grid(Grid_T *const grid)
{
  char name[20] = {'\0'};
  Collocation_T c;
  Basis_T b;
  Uint i;
  
  FOR_ALL(i,grid->patch)
  {
    struct Ret_S ret;
    Patch_T *const patch = grid->patch[i];
    c = UNDEFINED_COLLOCATION;
    b = UNDEFINED_BASIS;
    Uint n;
    
    /* filling grid */
    patch->grid = grid;
    
    /* filling patch number */
    patch->pn = i;
    
    /* filling inner boundary */
    patch->innerB = 0;
    
    /* filling name */
    sprintf(name,"box%d",i);
    patch->name = dup_s(name);
    
    /* filling n */
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    
    /* check for override */
    make_keyword_parameter(&ret,name,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;

    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
      
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    /* filling center */
    make_keyword_parameter(&ret,name,"center");
    patch->c[0] = Pgetd(ret.s0);
    patch->c[1] = Pgetd(ret.s1);
    patch->c[2] = Pgetd(ret.s2);
    
    /* filling size */
    make_keyword_parameter(&ret,name,"size");
    patch->s[0] = Pgetd(ret.s0);
    patch->s[1] = Pgetd(ret.s1);
    patch->s[2] = Pgetd(ret.s2);
    
    /* filling min: min = center-l/2 */
    patch->min[0] = patch->c[0]-patch->s[0]/2;
    patch->min[1] = patch->c[1]-patch->s[1]/2;
    patch->min[2] = patch->c[2]-patch->s[2]/2;
    
    /* filling max: max = center+l/2 */
    patch->max[0] = patch->c[0]+patch->s[0]/2;
    patch->max[1] = patch->c[1]+patch->s[1]/2;
    patch->max[2] = patch->c[2]+patch->s[2]/2;
    
    /* filling flags */
    patch->coordsys = Cartesian;
    
    /* collocation */
    patch->collocation[0] = get_collocation(PgetsEZ("collocation_a"));
    patch->collocation[1] = get_collocation(PgetsEZ("collocation_b"));
    patch->collocation[2] = get_collocation(PgetsEZ("collocation_c"));
  
    /* check for override */
    make_keyword_parameter(&ret,name,"collocation");
    c = get_collocation(PgetsEZ(ret.s0));
    if (c != UNDEFINED_COLLOCATION)
      patch->collocation[0] = c;
    c = get_collocation(PgetsEZ(ret.s1));
    if (c != UNDEFINED_COLLOCATION)
      patch->collocation[1] = c;
    c = get_collocation(PgetsEZ(ret.s2));
    if (c != UNDEFINED_COLLOCATION)
      patch->collocation[2] = c;
    
    assert(patch->collocation[0] != UNDEFINED_COLLOCATION);
    assert(patch->collocation[1] != UNDEFINED_COLLOCATION);
    assert(patch->collocation[2] != UNDEFINED_COLLOCATION);
    
    /* basis */
    patch->basis[0] = get_basis(PgetsEZ("basis_a"));
    patch->basis[1] = get_basis(PgetsEZ("basis_b"));
    patch->basis[2] = get_basis(PgetsEZ("basis_c"));
  
    /* check for override */
    make_keyword_parameter(&ret,name,"basis");
    b = get_basis(PgetsEZ(ret.s0));
    if ( b != UNDEFINED_BASIS)
      patch->basis[0] = b;
    b = get_basis(PgetsEZ(ret.s1));
    if ( b != UNDEFINED_BASIS)
      patch->basis[1] = b;
    b = get_basis(PgetsEZ(ret.s2));
    if ( b != UNDEFINED_BASIS)
      patch->basis[2] = b;
    
    assert(patch->basis[0] != UNDEFINED_BASIS);
    assert(patch->basis[1] != UNDEFINED_BASIS);
    assert(patch->basis[2] != UNDEFINED_BASIS);
    
  }
  
}

/* populating properties of a box patch for split cubed spherical.
// like: filling_box, central_box. */
void 
populate_box_patch_SplitCS
  (
  Grid_T *const grid,
  const char *const obj0,/* filling_box,central_box. */
  const Flag_T dir0,/* direction */
  const char *const region/* covering region */
  )

{
  const Uint Nsd[3] = {(Uint)Pgeti(P_"SplitCS_Nsplit_a"),
                       (Uint)Pgeti(P_"SplitCS_Nsplit_b"),
                       (Uint)Pgeti(P_"SplitCS_Nsplit_c")};
  char par[STR_SIZE3]  = {'\0'};
  char obj[STR_SIZE1]  = {'\0'};
  char name[STR_SIZE3] = {'\0'};
  const char *const dir = StrSide[dir0];
  Uint d0,d1,d2;
  
  /* object name */
  set_object_name_split_CS(obj,obj0);
  
  for (d0 = 0; d0 < Nsd[0]; d0++)
  {
    for (d1 = 0; d1 <  Nsd[1]; d1++)
    {
      for (d2 = 0; d2 <  Nsd[2]; d2++)
      {
        Patch_T *const patch = calloc(1,sizeof(*patch));
        IsNull(patch);
        grid->patch = 
          realloc(grid->patch,(grid->np+2)*sizeof(*grid->patch));
        IsNull(grid->patch);
        grid->patch[grid->np]   = patch;
        grid->patch[grid->np+1] = 0;
        
        /* filling grid */
        patch->grid = grid;
        
        /* filling patch number */
        patch->pn = grid->np;
        
        /* increase number patch */
        grid->np += 1;
          
        Flag_T side = dir0;
        
        assert(StrSide[side]);
        
        /* cover region */
        if (strcmp_i(region,"NS") || strcmp_i(region,"BH"))
        {
          sprintf(patch->CoordSysInfo->region,
            "(%s)(%s_%s)(%s_%s)",region,dir,region,dir,obj);
        }
        else if (strcmp_i(region,"NS1") || strcmp_i(region,"NS2"))
        {
          sprintf(patch->CoordSysInfo->region,
            "(%s)(%s_%s)(%s_%s)",region,dir,region,dir,obj);
        }
        else if (strcmp_i(region,"BH1") || strcmp_i(region,"BH2"))
        {
          sprintf(patch->CoordSysInfo->region,
            "(%s)(%s_%s)(%s_%s)",region,dir,region,dir,obj);
        }
        else if (strcmp_i(region,"filling_box"))
        {
          sprintf(patch->CoordSysInfo->region,
            "(%s)",region);
        }
        else
          Error0(NO_OPTION);
        
        
        /* filling n */
        patch->n[0] = (Uint)Pgeti(P_"SplitCS_n_a");
        patch->n[1] = (Uint)Pgeti(P_"SplitCS_n_b");
        patch->n[2] = (Uint)Pgeti(P_"SplitCS_n_c");
        
        /* filling nn */
        patch->nn = total_nodes_patch(patch);
        
        /* filling number of split */
        patch->nsplit[0] = Nsd[0];
        patch->nsplit[1] = Nsd[1];
        patch->nsplit[2] = Nsd[2];

        /* filling name */
        SCS_par_name(name);
        patch->name = dup_s(name);
        
        /* filling center */
        SCS_par_box_center(par,"a");
        patch->c[0] = Pgetd(par);
        SCS_par_box_center(par,"b");
        patch->c[1] = Pgetd(par);
        SCS_par_box_center(par,"c");
        patch->c[2] = Pgetd(par);
        
        /* filling size */
        SCS_par_box_length(par,"l");
        patch->s[0] = Pgetd(par);
        SCS_par_box_length(par,"w");
        patch->s[1] = Pgetd(par);
        SCS_par_box_length(par,"h");
        patch->s[2] = Pgetd(par);
        
        /* filling min: min = center-l/2 */
        patch->min[0] = patch->c[0]-patch->s[0]/2;
        patch->min[1] = patch->c[1]-patch->s[1]/2;
        patch->min[2] = patch->c[2]-patch->s[2]/2;

        /* filling max: max = center+l/2 */
        patch->max[0] = patch->c[0]+patch->s[0]/2;
        patch->max[1] = patch->c[1]+patch->s[1]/2;
        patch->max[2] = patch->c[2]+patch->s[2]/2;

        /* filling flags */
        patch->coordsys = Cartesian;
        
        /* collocation */
        patch->collocation[0] = Chebyshev_Extrema;
        patch->collocation[1] = Chebyshev_Extrema;
        patch->collocation[2] = Chebyshev_Extrema;
        
        /* basis */
        patch->basis[0] = Chebyshev_Tn_BASIS;
        patch->basis[1] = Chebyshev_Tn_BASIS;
        patch->basis[2] = Chebyshev_Tn_BASIS;
      }
    }
  }
}

/* populating properties of the box at the middle of left NS */
void populate_left_NS_central_box(Grid_T *const grid,const Uint pn)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,PATCH_NAME_PRT_P_"left_central_box",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  sprintf(var,PATCH_NAME_PRT_P_"left_central_box_n_a",grid->gn);
  patch->n[0] = (Uint)Pgeti(var);
  
  sprintf(var,PATCH_NAME_PRT_P_"left_central_box_n_b",grid->gn);
  patch->n[1] = (Uint)Pgeti(var);
  
  sprintf(var,PATCH_NAME_PRT_P_"left_central_box_n_c",grid->gn);
  patch->n[2] = (Uint)Pgeti(var);
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_a",grid->gn);
  patch->c[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_b",grid->gn);
  patch->c[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_c",grid->gn);
  patch->c[2] = Pgetd(var);
  
  /* filling size */
  sprintf(var,PATCH_NAME_PRT_P_"left_central_box_size_a",grid->gn);
  patch->s[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"left_central_box_size_b",grid->gn);
  patch->s[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"left_central_box_size_c",grid->gn);
  patch->s[2] = Pgetd(var);
  
  /* filling min: min = center-l/2 */
  patch->min[0] = patch->c[0]-patch->s[0]/2;
  patch->min[1] = patch->c[1]-patch->s[1]/2;
  patch->min[2] = patch->c[2]-patch->s[2]/2;
  
  /* filling max: max = center+l/2 */
  patch->max[0] = patch->c[0]+patch->s[0]/2;
  patch->max[1] = patch->c[1]+patch->s[1]/2;
  patch->max[2] = patch->c[2]+patch->s[2]/2;
  
  /* filling flags */
  patch->coordsys = Cartesian;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of the box at the middle of right BH */
void populate_right_BH_central_box(Grid_T *const grid,const Uint pn)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,PATCH_NAME_PRT_P_"right_central_box",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_n_a",grid->gn);
  patch->n[0] = (Uint)Pgeti(var);
  
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_n_b",grid->gn);
  patch->n[1] = (Uint)Pgeti(var);
  
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_n_c",grid->gn);
  patch->n[2] = (Uint)Pgeti(var);
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,PATCH_NAME_PRT_P_"right_BH_center_a",grid->gn);
  patch->c[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_BH_center_b",grid->gn);
  patch->c[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_BH_center_c",grid->gn);
  patch->c[2] = Pgetd(var);
  
  /* filling size */
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_a",grid->gn);
  patch->s[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_b",grid->gn);
  patch->s[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_c",grid->gn);
  patch->s[2] = Pgetd(var);
  
  /* filling min: min = center-l/2 */
  patch->min[0] = patch->c[0]-patch->s[0]/2;
  patch->min[1] = patch->c[1]-patch->s[1]/2;
  patch->min[2] = patch->c[2]-patch->s[2]/2;
  
  /* filling max: max = center+l/2 */
  patch->max[0] = patch->c[0]+patch->s[0]/2;
  patch->max[1] = patch->c[1]+patch->s[1]/2;
  patch->max[2] = patch->c[2]+patch->s[2]/2;
  
  /* filling flags */
  patch->coordsys = Cartesian;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}


/* populating properties of the box at the middle of central NS */
void populate_central_NS_central_box(Grid_T *const grid,const Uint pn)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,PATCH_NAME_PRT_P_"central_box",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  sprintf(var,PATCH_NAME_PRT_P_"central_box_n_a",grid->gn);
  patch->n[0] = (Uint)Pgeti(var);
  
  sprintf(var,PATCH_NAME_PRT_P_"central_box_n_b",grid->gn);
  patch->n[1] = (Uint)Pgeti(var);
  
  sprintf(var,PATCH_NAME_PRT_P_"central_box_n_c",grid->gn);
  patch->n[2] = (Uint)Pgeti(var);
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,PATCH_NAME_PRT_P_"NS_center_a",grid->gn);
  patch->c[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"NS_center_b",grid->gn);
  patch->c[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"NS_center_c",grid->gn);
  patch->c[2] = Pgetd(var);
  
  /* filling size */
  sprintf(var,PATCH_NAME_PRT_P_"central_box_size_a",grid->gn);
  patch->s[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"central_box_size_b",grid->gn);
  patch->s[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"central_box_size_c",grid->gn);
  patch->s[2] = Pgetd(var);
  
  /* filling min: min = center-l/2 */
  patch->min[0] = patch->c[0]-patch->s[0]/2;
  patch->min[1] = patch->c[1]-patch->s[1]/2;
  patch->min[2] = patch->c[2]-patch->s[2]/2;
  
  /* filling max: max = center+l/2 */
  patch->max[0] = patch->c[0]+patch->s[0]/2;
  patch->max[1] = patch->c[1]+patch->s[1]/2;
  patch->max[2] = patch->c[2]+patch->s[2]/2;
  
  /* filling flags */
  patch->coordsys = Cartesian;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of right box in for single neutron star */
void populate_right_box_sns(Grid_T *const grid,const Uint pn)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,PATCH_NAME_PRT_P_"right_box",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  sprintf(var,PATCH_NAME_PRT_P_"right_box_n_a",grid->gn);
  patch->n[0] = (Uint)Pgeti(var);
  
  sprintf(var,PATCH_NAME_PRT_P_"right_box_n_b",grid->gn);
  patch->n[1] = (Uint)Pgeti(var);
  
  sprintf(var,PATCH_NAME_PRT_P_"right_box_n_c",grid->gn);
  patch->n[2] = (Uint)Pgeti(var);
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,PATCH_NAME_PRT_P_"right_box_center_a",grid->gn);
  patch->c[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_box_center_b",grid->gn);
  patch->c[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_box_center_c",grid->gn);
  patch->c[2] = Pgetd(var);
  
  /* filling size */
  sprintf(var,PATCH_NAME_PRT_P_"right_box_size_a",grid->gn);
  patch->s[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_box_size_b",grid->gn);
  patch->s[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_box_size_c",grid->gn);
  patch->s[2] = Pgetd(var);
  
  /* filling min: min = center-l/2 */
  patch->min[0] = patch->c[0]-patch->s[0]/2;
  patch->min[1] = patch->c[1]-patch->s[1]/2;
  patch->min[2] = patch->c[2]-patch->s[2]/2;
  
  /* filling max: max = center+l/2 */
  patch->max[0] = patch->c[0]+patch->s[0]/2;
  patch->max[1] = patch->c[1]+patch->s[1]/2;
  patch->max[2] = patch->c[2]+patch->s[2]/2;
  
  /* filling flags */
  patch->coordsys = Cartesian;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}


/* populating properties of the filling box in cubed spherical grid */
void populate_filling_box_CubedSpherical(Grid_T *const grid,const Uint pn,const Flag_T side)
{
  Patch_T *const patch = grid->patch[pn];
  Uint n;
  double l;/* length */
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  struct Ret_S ret;
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling n */
  patch->n[0] = (Uint)PgetiEZ("n_a");
  patch->n[1] = (Uint)PgetiEZ("n_b");
  patch->n[2] = (Uint)PgetiEZ("n_c");
  
  switch(side)
  {
    case UP:
    /* filling name */
    sprintf(name,PATCH_NAME_PRT_P_"filling_box_up",grid->gn);
    patch->name = dup_s(name);
    
    sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
    l = Pgetd(var);
    
    /* filling center */
    patch->c[0] = 0;
    patch->c[1] = 0;
    patch->c[2] = 3./4.*l;
    
    /* filling size */
    patch->s[0] = l;
    patch->s[1] = 2*l;
    patch->s[2] = 0.5*l;
  
    /* check for override n*/
    sprintf(var,"Outermost0");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)   patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)   patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)   patch->n[2] = n/2;
    
    break;
    case DOWN:
    /* filling name */
    sprintf(name,PATCH_NAME_PRT_P_"filling_box_down",grid->gn);
    patch->name = dup_s(name);
    
    sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
    l = Pgetd(var);
    
    /* filling center */
    patch->c[0] = 0;
    patch->c[1] = 0;
    patch->c[2] = -3./4.*l;
    
    /* filling size */
    patch->s[0] = l;
    patch->s[1] = 2*l;
    patch->s[2] = 0.5*l;
  
    /* check for override n*/
    sprintf(var,"Outermost0");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)   patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)   patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)   patch->n[2] = n/2;
    
    break;
    case BACK:
    /* filling name */
    sprintf(name,PATCH_NAME_PRT_P_"filling_box_back",grid->gn);
    patch->name = dup_s(name);
    
    sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
    l = Pgetd(var);
    
    /* filling center */
    patch->c[0] = -3./4.*l;
    patch->c[1] = 0;
    patch->c[2] = 0;
    
    /* filling size */
    patch->s[0] = 0.5*l;
    patch->s[1] = 2*l;
    patch->s[2] = 2*l;
  
    /* check for override n*/
    sprintf(var,"Outermost0");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)   patch->n[0] = n/2;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)   patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)   patch->n[2] = n;
    
    break;
    case FRONT:
    /* filling name */
    sprintf(name,PATCH_NAME_PRT_P_"filling_box_front",grid->gn);
    patch->name = dup_s(name);
    
    sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
    l = Pgetd(var);
    
    /* filling center */
    patch->c[0] = 3./4.*l;
    patch->c[1] = 0;
    patch->c[2] = 0;
    
    /* filling size */
    patch->s[0] = 0.5*l;
    patch->s[1] = 2*l;
    patch->s[2] = 2*l;
    
    /* check for override n*/
    sprintf(var,"Outermost0");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)   patch->n[0] = n/2;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)   patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)   patch->n[2] = n;
    
    break;
    default:
      Error0(NO_OPTION);
  }
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling min: min = center-l/2 */
  patch->min[0] = patch->c[0]-patch->s[0]/2;
  patch->min[1] = patch->c[1]-patch->s[1]/2;
  patch->min[2] = patch->c[2]-patch->s[2]/2;
  
  /* filling max: max = center+l/2 */
  patch->max[0] = patch->c[0]+patch->s[0]/2;
  patch->max[1] = patch->c[1]+patch->s[1]/2;
  patch->max[2] = patch->c[2]+patch->s[2]/2;
  
  /* filling flags */
  patch->coordsys = Cartesian;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* populating properties of the box at the middle of right NS  */
void populate_right_NS_central_box(Grid_T *const grid,const Uint pn)
{
  Patch_T *const patch = grid->patch[pn];
  char name[100] = {'\0'};
  char var[100] = {'\0'};
  
  /* filling grid */
  patch->grid = grid;
  
  /* filling patch number */
  patch->pn = pn;
  
  /* filling inner boundary */
  patch->innerB = 0;
  
  /* filling name */
  sprintf(name,PATCH_NAME_PRT_P_"right_central_box",grid->gn);
  patch->name = dup_s(name);
  
  /* filling n */
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_n_a",grid->gn);
  patch->n[0] = (Uint)Pgeti(var);
  
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_n_b",grid->gn);
  patch->n[1] = (Uint)Pgeti(var);
  
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_n_c",grid->gn);
  patch->n[2] = (Uint)Pgeti(var);
  
  /* filling nn */
  patch->nn = total_nodes_patch(patch);
  
  /* filling center */
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_a",grid->gn);
  patch->c[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_b",grid->gn);
  patch->c[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_c",grid->gn);
  patch->c[2] = Pgetd(var);
  
  /* filling size */
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_a",grid->gn);
  patch->s[0] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_b",grid->gn);
  patch->s[1] = Pgetd(var);
  sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_c",grid->gn);
  patch->s[2] = Pgetd(var);
  
  /* filling min: min = center-l/2 */
  patch->min[0] = patch->c[0]-patch->s[0]/2;
  patch->min[1] = patch->c[1]-patch->s[1]/2;
  patch->min[2] = patch->c[2]-patch->s[2]/2;
  
  /* filling max: max = center+l/2 */
  patch->max[0] = patch->c[0]+patch->s[0]/2;
  patch->max[1] = patch->c[1]+patch->s[1]/2;
  patch->max[2] = patch->c[2]+patch->s[2]/2;
  
  /* filling flags */
  patch->coordsys = Cartesian;
  
 /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;
  
  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
}

/* memory alloc patches for Cartesian grid type */
void alloc_patches_Cartesian_grid(Grid_T *const grid)
{
  Uint Nboxes;/* number of boxes */
  Uint i;
  
  if (get_parameter("number_of_boxes") == 0)
    Error0("\"number_of_boxes\" parameter is not defined!\n");
    
  Nboxes = (Uint) Pgeti("number_of_boxes");
  
  grid->patch = calloc((Nboxes+1),sizeof(*grid->patch));
  IsNull(grid->patch);
  
  for (i = 0; i < Nboxes; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    IsNull(grid->patch[i]);
  }
  
}
