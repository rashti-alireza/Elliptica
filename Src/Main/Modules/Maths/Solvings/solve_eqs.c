/*
// Alireza Rashti
// August 2018
*/

#include "solve_eqs.h"

/* solving equations
// ->return value: EXIT_SUCCESS
*/
int solve_eqs(Grid_T *const grid)
{
  fSolve_T *solve = 0;
  
  /* choosing solving method */
  if (strcmp_i(GetParameterS_E("Solving_Method"),"Parallel_Patch"))
    solve = parallel_patch_method;
  else
    abortEr_s("No such method \"%s\" defined for this function.\n",
      GetParameterS("Solving_Method"));
  
  /* call the specific solving method */
  solve(grid);
  
  return EXIT_SUCCESS;
}

/* threads are sprawned over patches and 
// each solves equations in their region (patch)
// ->return value: EXIT_SUCCESS
*/
static int parallel_patch_method(Grid_T *const grid)
{
  unsigned p;
  const double tol = GetParameterD_E("Solving_Tolerance");
  double current_tol_grid = 0;//calculate_tolerance_on_grid(grid);
  
  while (!GRT(current_tol_grid,tol))
  {
    unsigned cn = 0;
    OMP_PARALLEL_PATCH(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *patch = grid->patch[p];
      /* making b in ax = b */
      b_in_ax_b_whole_ppm(patch,cn);
      b_in_ax_b_bndry_ppm(patch,cn);
      /* making a in ax = b */
      a_in_ax_b_whole_ppm(patch,cn);
      //a_in_ax_b_bndry_ppm(patch,cn);
      /* solve ax = b */
      //solve_ax_b_ppm(patch,cn);
    }
    break;
    //current_tol_grid = calculate_tolerance_on_grid(grid);
  }
  
  return EXIT_SUCCESS;
}

/* filling a in a.x = b for whole points of the given patch
// for the given colloction of fields cn according to the equations for a.
// THREAD SAFE.
// ->return value: EXIT_SUCCESS
*/
static int a_in_ax_b_whole_ppm(Patch_T *const patch,const unsigned cn)
{
  Solve_T *const slv = patch->solution_man->solve[cn];
  unsigned i;
  
  for (i = 0; i < slv->nf; ++i)
  {
    Field_T *field        = slv->field[i];
    fEquation_T *field_eq = slv->field_eq[i];
    double *b             = &slv->b[slv->f_occupy[i]];
    
    /* fill b in ax = b for the specified field */
    field_eq(field,b);
  }
  
  return EXIT_SUCCESS;
}

/* filling b in a.x = b for whole points of the given patch
// for the given colloction of fields cn according to field equations.
// THREAD SAFE.
// ->return value: EXIT_SUCCESS
*/
static int b_in_ax_b_whole_ppm(Patch_T *const patch,const unsigned cn)
{
  Solve_T *const slv = patch->solution_man->solve[cn];
  unsigned i;
  
  for (i = 0; i < slv->nf; ++i)
  {
    Field_T *field        = slv->field[i];
    fEquation_T *field_eq = slv->field_eq[i];
    double *b             = &slv->b[slv->f_occupy[i]];
    
    /* fill b in ax = b for the specified field */
    field_eq(field,b);
  }
  
  return EXIT_SUCCESS;
}

/* filling b in a.x = b for boundary points of the given patch
// for the given collocation of fields cn according to the boundary 
// conditions between adjacent patches and boundary conditions equations.
// THREAD SAFE.
// ->return value: EXIT_SUCCESS
*/
static int b_in_ax_b_bndry_ppm(Patch_T *const patch,const unsigned cn)
{
  const unsigned nintfc = countf(patch->interface);
  Solve_T *const solve = patch->solution_man->solve[cn];
  unsigned intfc;
  
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    unsigned nsfc = interface->ns;
    unsigned sfc;
    
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      Boundary_Condition_T BC;
      
      BC.subface = subface;
      BC.solve   = solve;
      BC.cn      = cn;
      
      if (!subface->exterF)/* if subface is internal */
      {
        continue;
      }
      else if (subface->innerB)/* if there is inner boundary */
      {
        continue;/* this part can be develope for general cases */
      }
      else if (subface->outerB)/* if it reaches outer boundary */
      {
        bndry_outerB_ppm(&BC);
      }
      else if (subface->touch)/* if two patches are in touch */
      {
        if (subface->copy)/* if the collocated point */
        {
          /* copy values */
          bndry_copy_ppm(&BC);
        }
        else
        {
          /* interpolate values */
          bndry_interpolate_ppm(&BC);
        }
      }
      else /* if there is an overlap case */
      {
        /* interpolate values */
        bndry_interpolate_ppm(&BC);
      }
      
    }/* end of for (sfc = 0; sfc < nsfc; ++sfc) */
    
  }/* end of for (intfc = 0; intfc < nintfc; ++intfc) */
  
  return EXIT_SUCCESS;
}

/* copy values of fields or their derivative
// along the normaln vector between collocated points 
// at the boundary between two patches for parallel patch method.
// ->return value: EXIT_SUCCESS
*/
static int bndry_copy_ppm(Boundary_Condition_T *const bc)
{
  SubFace_T *const subface = bc->subface;
  Solve_T *const solve     = bc->solve;
  Grid_T *const grid       = subface->patch->grid;
  Patch_T *const patch_adj = grid->patch[subface->adjPatch];
  const unsigned np        = subface->np;
  unsigned *const id       = subface->id;
  unsigned *const adjid    = subface->adjid;
  Field_T *field;
  Field_T *field_adj;
  double *b;
  unsigned n;
  unsigned i;
  
  /* df/dn = df/dn|adjacent */
  if (subface->df_dn)
  {
    Patch_T patch_tmp = make_temp_patch(patch_adj);/* for thread safety purposes */
    Field_T *field_tmp;/* for thread safety purposes */
    double *Nvec;/* normal vector */
    const char *der0 = "x",*der1 = "y",*der2 = "z";
    
    /* if one wants to solve whole equations on curvilinear patches */
    if (!strcmp_i(GetParameterS("Solving_Interpolation_Normal"),"Cartesian_Normal"))
    {
      der0 = "a";
      der1 = "b";
      der2 = "c";
      normal_vec_curvilinear(0);/* NOTE: it should be defined a pointer to function 
                                // since it is called for each points */
    }
    
    for (i = 0; i < solve->nf; ++i)
    {
      double *f_a = 0,*f_b = 0,*f_c = 0;
      double *f_a_adj = 0,*f_b_adj = 0,*f_c_adj = 0;
      field     = solve->field[i];
      field_adj = patch_adj->solution_man->solve[bc->cn]->field[i];
      field_tmp = add_field("tmp_field","(3dim)",&patch_tmp,NO);
      field_tmp->v = field_adj->v;/* soft copying values of fields */
      b = &solve->b[solve->f_occupy[i]];
      
      f_a     = Partial_Derivative(field,der0);
      f_b     = Partial_Derivative(field,der1);
      f_c     = Partial_Derivative(field,der2);
      f_a_adj = Partial_Derivative(field_tmp,der0);
      f_b_adj = Partial_Derivative(field_tmp,der1);
      f_c_adj = Partial_Derivative(field_tmp,der2);
      
      for (n = 0; n < np; ++n)
      {
        Point_T point;
        point.ind   = id[n];
        point.patch = subface->patch;
        point.face  = subface->face;
        Nvec = normal_vec(&point);
        
        b[id[n]]   = Nvec[0]*(f_a[id[n]] - f_a_adj[adjid[n]]) +
                     Nvec[1]*(f_b[id[n]] - f_a_adj[adjid[n]]) +
                     Nvec[2]*(f_c[id[n]] - f_a_adj[adjid[n]]) ;
      }
      /* freeing memories */
      field_tmp ->v = 0;/* since pointing to field_adj->v */
      remove_field(field_tmp);
      free(field_tmp);
      free(f_a);
      free(f_b);
      free(f_c);
      free(f_a_adj);
      free(f_b_adj);
      free(f_c_adj);
    }/* end of for (i = 0; i < solve->nf; ++i) */
  }/* if (subface->df_dn) */
  
  /* f = f|adjacent */
  else
  {
    for (i = 0; i < solve->nf; ++i)
    {
      field     = solve->field[i];
      field_adj = patch_adj->solution_man->solve[bc->cn]->field[i];
      b         = &solve->b[solve->f_occupy[i]];
      
      for (n = 0; n < np; ++n)
      {
        b[id[n]]   = field->v[id[n]] - field_adj->v[adjid[n]];
      }
      
    }/* end of for (i = 0; i < solve->nf; ++i) */
  }
  
  return EXIT_SUCCESS;
}


/* interpolating values of fields at the interface of 
// two patches for parallel patch method.
// ->return value: EXIT_SUCCESS
*/
static int bndry_interpolate_ppm(Boundary_Condition_T *const bc)
{
  SubFace_T *const subface   = bc->subface;
  Solve_T   *const solve     = bc->solve;
  Grid_T    *const grid      = subface->patch->grid;
  Patch_T   *const patch_adj = grid->patch[subface->adjPatch];
  unsigned   const np        = subface->np;
  unsigned  *const id        = subface->id;
  Patch_T patch_tmp = make_temp_patch(patch_adj);/* for thread safety purposes */
  Field_T *field_tmp;/* for thread safety purposes */
  Interpolation_T *interp_s = init_interpolation();
  unsigned i;
  
  interp_s->point->patch    = patch_adj;
  interp_s->point->touch    = subface->touch;
  interp_s->point->face     = subface->adjFace;
  interp_s->point->sameX    = subface->sameX;
  interp_s->point->sameY    = subface->sameY;
  interp_s->point->sameZ    = subface->sameZ;
  
  for (i = 0; i < solve->nf; ++i)
  {
    Field_T *field     = solve->field[i];
    Field_T *field_adj = patch_adj->solution_man->solve[bc->cn]->field[i];
    double *b = &solve->b[solve->f_occupy[i]];
    unsigned n;
    
    field_tmp = add_field("tmp_field","(3dim)",&patch_tmp,NO);
    field_tmp->v = field_adj->v;/* soft copying values of fields */
    interp_s->field = field_tmp;
     
    for (n = 0; n < np; ++n)
    {
      interp_s->point->x = 0;
      b[id[n]] = field->v[id[n]] - interpolation(interp_s);
    }
    
  }/* end of for (i = 0; i < solve->nf; ++i) */
  
  free_interpolation(interp_s);
  
  return EXIT_SUCCESS;
}

/* filling values of b in a.x = b for those points reach outer boundary
// according to boundary condition imposed by the equation for
// parallel patch method.
// ->return value: EXIT_SUCCESS
*/
static int bndry_outerB_ppm(Boundary_Condition_T *const bc)
{
  Solve_T *const solve = bc->solve;
  unsigned i;
  
  for (i = 0; i < solve->nf; ++i)
  {
    bc->field        	 = solve->field[i];
    bc->node		 = bc->subface->id;
    bc->nn 		 = bc->subface->np;
    fEquation_T *bc_eq   = solve->bc_eq[i];
    double *b            = &solve->b[solve->f_occupy[i]];
    
    /* fill b in ax = b for the specified field */
    bc_eq(bc,b);
  }
  
  return EXIT_SUCCESS;
}

/* finding normal vector at the interface of each patch 
// in a general coords which is not necessarily Cartesian coords.
// the idea is if one can cover the whole manifold with a coord sys
// such that each patch can be written in this coord sys,
// then the equations alos would be preferred to be written in this 
// coord sys; thus, this type of normal vector is needed.
// ->return value: the normal vector reported in a general coord.
*/
static double *normal_vec_curvilinear(Point_T *const point)
{
  UNUSED(point);
  abortEr("This part is not developed yet!\n");
  return 0;
}
