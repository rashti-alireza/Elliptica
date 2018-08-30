/*
// Alireza Rashti
// August 2018
*/

#include "solve_eqs_ppm.h"

/* threads are sprawned over patches and 
// each solves equations in their region (patch)
// ->return value: EXIT_SUCCESS
*/
int parallel_patch_method(Grid_T *const grid)
{
  unsigned p;
  /* residual determined in the input file */
  const double res_input = GetParameterD_E("Solving_Residual");
  
  PARALLEL_PATCH_METHOD_OpenMP(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    unsigned cn;/* collection number, which refers to 
                // collection of fields to be solved.
                */
    
    for (cn = 0; cn < patch->solution_man->ns; ++cn)
    {
      double res_patch;/* residual of each collection at this patch */
      double *b = patch->solution_man->solve[cn]->b;
      /* number of components of b */
      unsigned nb = patch->nn*patch->solution_man->solve[cn]->nf;
      
      /* making b in ax = b */
      b_in_ax_b_whole_ppm(patch,cn);
      b_in_ax_b_bndry_ppm(patch,cn);
      
      /* if root means square of b is less than res_input skip otherwise keep going */
      res_patch = rms(nb,b,0);
      
      if (LSS(res_patch,res_input))
        continue;
      
      /* making a in ax = b */
      a_in_ax_b_whole_ppm(patch,cn);
      a_in_ax_b_bndry_ppm(patch,cn);
      
      /* solve ax = b */
      solve_ax_b_ppm(patch,cn);
    }
  }/* end of for (p = 0; p < grid->np; ++p) */
  
  return EXIT_SUCCESS;
}

/* solving a.x = b.
// THREAD SAFE.
// ->return value: EXIT_SUCCESS
*/
static int solve_ax_b_ppm(Patch_T *const patch,const unsigned cn)
{
  Solve_T *const solve = patch->solution_man->solve[cn];
  
  if (strcmp_i(GetParameterS_E("Linear_Solver"),"UMFPACK"))
  {
    UmfPack_T umfpack[1] = {0};
    umfpack->a = solve->a;
    umfpack->b = solve->b;
    umfpack->x = solve->x;
    direct_solver_umfpack_di(umfpack);
    
  }
  else if (strcmp_i(GetParameterS_E("Linear_Solver"),"UMFPACK_long"))
  {
    UmfPack_T umfpack[1] = {0};
    umfpack->a = solve->a;
    umfpack->b = solve->b;
    umfpack->x = solve->x;
    direct_solver_umfpack_dl(umfpack);
  }
  else
    abortEr_s("There is no such \"%s\" solver defined.\n",
        GetParameterS_E("Linear_Solver"));
        
  return EXIT_SUCCESS;
}

/* filling a in a.x = b for the bndry points of the given patch
// for the given colloction of fields cn according to the equations for a.
// THREAD SAFE.
// ->return value: EXIT_SUCCESS
*/
static int a_in_ax_b_bndry_ppm(Patch_T *const patch,const unsigned cn)
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
        abortEr(INCOMPLETE_FUNC);
      }
      else if (subface->outerB)/* if it reaches outer boundary */
      {
        a_bndry_outerB_ppm(&BC);
      }
      else if (subface->touch)/* if two patches are in touch */
      {
        if (subface->copy)/* if the collocated point */
        {
          /* copy values */
          a_bndry_copy_ppm(&BC);
        }
        else
        {
          /* interpolate values */
          a_bndry_interpolate_ppm(&BC);
        }
      }
      else /* if there is an overlap case */
      {
        /* interpolate values */
        a_bndry_interpolate_ppm(&BC);
      }
      
    }/* end of for (sfc = 0; sfc < nsfc; ++sfc) */
    
  }/* end of for (intfc = 0; intfc < nintfc; ++intfc) */
  
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
  Jacobian_Eq_T *jac = calloc(1,sizeof(*jac));
  pointerEr(jac);
  
  /* filling up the jacobian struct */
  jac->field = slv->field;
  jac->nf    = slv->nf;
  
  /* getting the pertinent equation for jacobian */
  fEquation_T *const jacobian  = slv->jacobian_eq;
  
  /* fill a in a.x = b for the specified field */
  slv->a = jacobian(jac,patch);
  
  /* freeing memory */
  free(jac);
  
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
        abortEr(INCOMPLETE_FUNC);
      }
      else if (subface->outerB)/* if it reaches outer boundary */
      {
        b_bndry_outerB_ppm(&BC);
      }
      else if (subface->touch)/* if two patches are in touch */
      {
        if (subface->copy)/* if the collocated point */
        {
          /* copy values */
          b_bndry_copy_ppm(&BC);
        }
        else
        {
          /* interpolate values */
          b_bndry_interpolate_ppm(&BC);
        }
      }
      else /* if there is an overlap case */
      {
        /* interpolate values */
        b_bndry_interpolate_ppm(&BC);
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
static int b_bndry_copy_ppm(Boundary_Condition_T *const bc)
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


/* making those parts of jacobian related to the section in b in which 
// it copys values of fields or their derivative
// along the normaln vector between collocated points 
// at the boundary between two patches for parallel patch method.
// ->return value: EXIT_SUCCESS
*/
static int a_bndry_copy_ppm(Boundary_Condition_T *const bc)
{
  SubFace_T *const subface = bc->subface;
  Solve_T *const solve     = bc->solve;
  Patch_T *const patch = subface->patch;
  const unsigned nb        = subface->np;
  const unsigned nn = bc->subface->patch->nn;
  const unsigned *const bndry = subface->id;
  double **const a = solve->a->reg->A;
  double *Nvec = 0;
  unsigned i;
  
  /* df/dn = df/dn|adjacent */
  if (subface->df_dn)
  {
    fJs_T *const jf = patch->solution_man->j_func;
    Matrix_T *j_x = 0,*j_y = 0,*j_z = 0;
    Point_T point;
    point.patch = patch;
    point.face  = subface->face;
    
    /* if one wants to solve whole equations on curvilinear patches */
    if (!strcmp_i(GetParameterS("Solving_Interpolation_Normal"),"Cartesian_Normal"))
    {
      abortEr(INCOMPLETE_FUNC);
      
      const char *types[] = {"j_a","j_b","j_c",0};
      prepare_Js_jacobian_eq(patch,types);
      j_x = get_j_matrix(patch,"j_a");
      j_y = get_j_matrix(patch,"j_b");
      j_z = get_j_matrix(patch,"j_c");
  
    }
    else
    {
      const char *types[] = {"j_x","j_y","j_z",0};
      prepare_Js_jacobian_eq(patch,types);
      j_x = get_j_matrix(patch,"j_x");
      j_y = get_j_matrix(patch,"j_y");
      j_z = get_j_matrix(patch,"j_z");
    }
    
    for (i = 0; i < solve->nf; ++i)
    {
      unsigned initial = solve->f_occupy[i];
      unsigned final   = initial+nn;
      unsigned ijk,lmn;
      
      /* fill a in a.x = b for the specified field */
      for (ijk = 0; ijk < nb; ++ijk)
      {
        point.ind   = bndry[ijk];
        Nvec = normal_vec(&point);
        for (lmn = initial; lmn < final;++lmn)
          a[initial+bndry[ijk]][lmn] = 
              Nvec[0]*jf(j_x,bndry[ijk],lmn-initial) +
              Nvec[1]*jf(j_y,bndry[ijk],lmn-initial) +
              Nvec[2]*jf(j_z,bndry[ijk],lmn-initial) ;
      }
    }/* end of for (i = 0; i < solve->nf; ++i) */
  }/* if (subface->df_dn) */
  
  /* f = f|adjacent */
  else
  {
    for (i = 0; i < solve->nf; ++i)
    {
      unsigned initial = solve->f_occupy[i];
      unsigned final   = initial+nn;
      unsigned ijk,lmn;
      
      /* fill a in a.x = b for the specified field */
      for (ijk = 0; ijk < nb; ++ijk)
      {
        for (lmn = initial; lmn < final;++lmn)
          a[initial+bndry[ijk]][lmn] = 0;
          
        a[initial+bndry[ijk]][initial+bndry[ijk]] = 1;

      }
      
    }/* end of for (i = 0; i < solve->nf; ++i) */
  }
  
  return EXIT_SUCCESS;
}

/* making a in a.x = b for sections which needs 
// interpolating values of fields at the interface of 
// two patches for parallel patch method.
// ->return value: EXIT_SUCCESS
*/
static int a_bndry_interpolate_ppm(Boundary_Condition_T *const bc)
{
  Solve_T *const solve = bc->solve;
  double **a = solve->a->reg->A;
  const unsigned *const bndry = bc->subface->id;
  const unsigned Nb = bc->subface->np;
  const unsigned nn = bc->subface->patch->nn;
  unsigned i;
  
  for (i = 0; i < solve->nf; ++i)
  {
    unsigned initial = solve->f_occupy[i];
    unsigned final   = initial+nn;
    unsigned ijk,lmn;
    
    /* fill a in a.x = b for the specified field */
    for (ijk = 0; ijk < Nb; ++ijk)
    {
      for (lmn = initial; lmn < final;++lmn)
        a[initial+bndry[ijk]][lmn] = 0;
      
      a[initial+bndry[ijk]][initial+bndry[ijk]] = 1;
    }
  }

  return EXIT_SUCCESS;
}

/* interpolating values of fields at the interface of 
// two patches for parallel patch method.
// ->return value: EXIT_SUCCESS
*/
static int b_bndry_interpolate_ppm(Boundary_Condition_T *const bc)
{
  SubFace_T *const subface   = bc->subface;
  Solve_T   *const solve     = bc->solve;
  Grid_T    *const grid      = subface->patch->grid;
  Patch_T   *const adj_patch = grid->patch[subface->adjPatch];
  unsigned   const np        = subface->np;
  unsigned  *const id        = subface->id;
  Patch_T adj_patch_tmp = make_temp_patch(adj_patch);/* for thread safety purposes */
  Field_T *field_tmp;/* for thread safety purposes */
  Interpolation_T *interp_s = init_interpolation();
  double X[3],*x;
  unsigned i;
  
  for (i = 0; i < solve->nf; ++i)
  {
    Field_T *field     = solve->field[i];
    Field_T *field_adj = adj_patch->solution_man->solve[bc->cn]->field[i];
    double *b = &solve->b[solve->f_occupy[i]];
    unsigned n;
    
    field_tmp = add_field("tmp_field","(3dim)",&adj_patch_tmp,NO);
    field_tmp->v = field_adj->v;/* soft copying values of fields */
    interp_s->field = field_tmp;
    
    fill_interpolation_flags(interp_s,subface);
    plan_interpolation(interp_s);
    
    for (n = 0; n < np; ++n)
    {
      x = subface->patch->node[id[n]]->x;
      X_of_x(X,x,&adj_patch_tmp);
      interp_s->X = X[0];
      interp_s->Y = X[1];
      interp_s->Z = X[2];
      b[id[n]] = field->v[id[n]] - execute_interpolation(interp_s);
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
static int b_bndry_outerB_ppm(Boundary_Condition_T *const bc)
{
  Solve_T *const solve = bc->solve;
  bc->node	       = bc->subface->id;
  bc->nn               = bc->subface->np;
  unsigned i;
  
  for (i = 0; i < solve->nf; ++i)
  {
    bc->field        	 = solve->field[i];
    fEquation_T *bc_eq   = solve->bc_eq[i];
    double *b            = &solve->b[solve->f_occupy[i]];
    
    /* fill b in a.x = b for the specified field */
    bc_eq(bc,b);
  }
  
  return EXIT_SUCCESS;
}

/* filling values of a in a.x = b for those points reach outer boundary
// for parallel patch method. note that since b equals to eq|i = f[i]-bc
// then d(eq|i)/df(j) = Dirac_Delta(i,j).
// ->return value: EXIT_SUCCESS
*/
static int a_bndry_outerB_ppm(Boundary_Condition_T *const bc)
{
  Solve_T *const solve = bc->solve;
  double **a = solve->a->reg->A;
  const unsigned *const bndry = bc->subface->id;
  const unsigned nb = bc->subface->np;
  const unsigned nn = bc->subface->patch->nn;
  unsigned i;
  
  for (i = 0; i < solve->nf; ++i)
  {
    unsigned initial = solve->f_occupy[i];
    unsigned final   = initial+nn;
    unsigned ijk,lmn;
    
    /* fill a in a.x = b for the specified field */
    for (ijk = 0; ijk < nb; ++ijk)
    {
      for (lmn = initial; lmn < final;++lmn)
        a[initial+bndry[ijk]][lmn] = 0;
      
      a[initial+bndry[ijk]][initial+bndry[ijk]] = 1;
    }
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

/* filling flags of iterpolation and based on subface. */
static void fill_interpolation_flags(Interpolation_T *const it,const SubFace_T *const sf)
{
  if (sf->sameX)
  {
    it->YZ_dir_flag = 1;
    it->I = const_index_on_AdjFace(sf);
  }
  else if (sf->sameY)
  {
    it->XZ_dir_flag = 1;
    it->J = const_index_on_AdjFace(sf);
  }
  else if (sf->sameZ)
  {
    it->XY_dir_flag = 1;
    it->K = const_index_on_AdjFace(sf);
  }
  else if (!sf->sameX && !sf->sameY && !sf->sameZ)
  {
    it->XYZ_dir_flag = 1;
  }
  else if (sf->sameX && sf->sameY && sf->sameZ)
  {
    abortEr("How come to have all the sameX,Y,Z flags of subface be the same.\n"
    "It means there is something wrong at finding of adjacent patches.\n");
  }
  
}

/* according to the adjFace where the subface touches it,
// it returns the index of those points on the "adjacent face" which
// not varying on the adjFace.
// ->return value: constant index(coords) of a given face.
*/
static unsigned const_index_on_AdjFace(const SubFace_T *const sf)
{
  const unsigned f = sf->adjFace;
  const Patch_T *const patch = sf->patch->grid->patch[sf->adjPatch];
  const unsigned *const n = patch->n;
  unsigned C = 0;/* constant value */
  
  switch(f)
  {
    case I_0:
      C = 0;
      break;
    case I_n0:
      C = n[0]-1;
      break; 
    case J_0:
      C = 0;
      break; 
    case J_n1:
      C = n[1]-1;
      break; 
    case K_0:
      C = 0;
      break; 
    case K_n2:
      C = n[2]-1;
      break;
    default:
      abortEr("There is not such interface.\n");
  }
  
  return C;
}