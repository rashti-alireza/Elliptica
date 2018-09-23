/*
// Alireza Rashti
// August 2018
*/

#include "solve_eqs_ppm.h"

/* Test: only uses Dirichlet BC */

/* threads are sprawned over patches and 
// each solves equations in their region (patch)
// ->return value: EXIT_SUCCESS
*/
int parallel_patch_method(Grid_T *const grid)
{
  /* residual determined in the input file */
  const double res_input = fabs(GetParameterD_E("Solving_Residual"));
  Flag_T IsItSolved = NO;
  const int NumIter = GetParameterI_E("Linear_Solver_Number_of_Iteration");
  int iter = 0;
  
  while (IsItSolved == NO && iter < NumIter)
  {
    unsigned p;
    
    /* initializing fields to be solved. */
    initialize_ppm(grid);
    
    IsItSolved = check_residual(grid,res_input);
    if (IsItSolved == YES)
      break;
      
    PARALLEL_PATCH_METHOD_OpenMP(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *patch = grid->patch[p];
      unsigned cn;/* collection number, which refers to 
                  // collection of fields to be solved.
                  */
      printf("Solving Equation(s) at patch = %s\n",patch->name);
      
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
        printf("L2-norm of F at J.x = -F in %s = %0.15f\n",patch->name,res_patch);
        if (LSS(res_patch,res_input))
          continue;
        
        /* making a in ax = b */
        if (strstr_i(GetParameterS_E("Making_Jacobian_For_Newton_Method"),"Finite_Difference")) 
        {
          a_in_ax_b_finite_diff_ppm(patch,cn);
        }
        else if (strstr_i(GetParameterS_E("Making_Jacobian_For_Newton_Method"),"Spectral"))
        {
          a_in_ax_b_whole_ppm(patch,cn);
          a_in_ax_b_bndry_ppm(patch,cn);
        }
        else
          abortEr("There is no such option for \" Making_Jacobian_For_Newton_Method\" \n");
        
        /* solve ax = b */
        solve_ax_b_ppm(patch,cn);
      }
    }/* end of for (p = 0; p < grid->np; ++p) */
    
    update_fields_ppm(grid);/*  update fields with new solution */
    iter++;
  }
  
  return EXIT_SUCCESS;
}

/* using finite difference to find a in a.x = b.
// a[i][j] = { (b2(u1,u2,...,uj+eps,..,un)-b1(u1,u2,...,un))/eps } |i.
*/
static int a_in_ax_b_finite_diff_ppm(Patch_T *const patch,const unsigned cn)
{
  abortEr(INCOMPLETE_FUNC"\n NOT THREAD SAFE.\n");
  
  Solve_T *const slv = patch->solution_man->solve[cn];
  unsigned nn = patch->nn;
  const double EPS = 1.0/nn;
  unsigned i,j,f;
  double *b1 = alloc_double(nn);
  
  if (slv->nf > 1)
    abortEr(INCOMPLETE_FUNC);
    
  for (f = 0; f < slv->nf; ++f)
  {
    Field_T *field  = slv->field[f];
    Matrix_T *jac   = alloc_matrix(REG_SF,nn,nn);
    double **J      = jac->reg->A;
    double *b2      = &slv->b[slv->f_occupy[f]];
    
    //b_in_ax_b_whole_ppm(patch,cn);
    //b_in_ax_b_bndry_ppm(patch,cn);
    
    for(i = 0; i < nn; ++i)
      b1[i] = b2[i];
     
    /* varying b by fields and filling up jacobian */
    for (j = 0; j < nn; ++j)
    {
      /* ERRORRRRR: NOT THREAD SAFE
      // not thread safe since when one changes 
      // the following field value, since the other thread needs this value
      // for interpolation and other boundary conditions, it messes up
      // thread safety.
      */
      field->v[j] += EPS;
      free_coeffs(field);
      
      b_in_ax_b_whole_ppm(patch,cn);
      b_in_ax_b_bndry_ppm(patch,cn);
      
      for (i = 0; i < nn; ++i)
      {
        //test
        //if (i == 13)
          //printf("b2=%f,b1=%f,diff=%f\n",b2[i],b1[i],b2[i]-b1[i]);
        //end
        J[i][j] = (b2[i]-b1[i])/EPS;
      }
        
      field->v[j] -= EPS;
      free_coeffs(field);
    }/* end of for (j = 0; j < nn; ++j) */

    slv->a = jac;
    
    /* since b is modified we need to retrieve it again */
    b_in_ax_b_whole_ppm(patch,cn);
    b_in_ax_b_bndry_ppm(patch,cn);
  }/* end of for (f = 0; f < slv->nf; ++f) */
  
  free(b1);
  return EXIT_SUCCESS;
}

/* find out the residual of each patch and decide weather the equations
// are already solved or not.
// ->return value: YES if EQs are solved, NO otherwise.
*/
static Flag_T check_residual(const Grid_T *const grid,const double res_input)
{
  Flag_T flg = YES;
  unsigned p;
  
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
      //test
      static int si = 0;
      
      fprintf(stderr,"%d %s,%g\n",si,patch->name,res_patch);
      si++;
      //end
      if (GRT(res_patch,res_input))
      {
        flg = NO;
        break;
      }
      
    }/* end of for (cn = 0; cn < patch->solution_man->ns; ++cn) */
    if (flg == NO)
      break;
  }/* end of for (p = 0; p < grid->np; ++p) */
  
  return flg;
}

/* having solved for fields, now one needs to 
// update them with found value.
// Note: in ppm method solving algorith is Newton so for update
// we have: J.(x_2-x_1) = -F -> J.(x_1-x_2) = F -> J.x = F
// => u_new = u_old-x.
*/
static void update_fields_ppm(Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Solution_Man_T *solution_man = patch->solution_man;
    unsigned s;
    
    for (s = 0; s < solution_man->ns; ++s)
    {
      Solve_T *solve = solution_man->solve[s];
      unsigned f;
      
      for (f = 0; f < solve->nf; ++f)
      {
        Field_T *field = patch->pool[Ind(solve->field[f]->name)];
        double *x = &solve->x[solve->f_occupy[f]];
        double *u_old = field->v;
        double *u_new = field->v;
        unsigned i;
        
        free_coeffs(field);
        
        for (i = 0; i < patch->nn; ++i)
          u_new[i] = u_old[i]-x[i];
      }
      
    }/* end of for (cn = 0; cn < solution_man->ns; ++cn) */
    
  }/* end of FOR_ALL_PATCHES(p,grid) */
}

/* initializing some fields in each Solve_T 
// and duplicating some memories to prohibit race condition/
*/
static void initialize_ppm(Grid_T *const grid)
{
  copy_initial_values_ppm(grid);
}

/* copying initial value of fields in patch->pool
// to solve->fields in each patch.
*/
static void copy_initial_values_ppm(Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Solution_Man_T *solution_man = patch->solution_man;
    unsigned s;
    
    for (s = 0; s < solution_man->ns; ++s)
    {
      Solve_T *solve = solution_man->solve[s];
      unsigned f;
      
      for (f = 0; f < solve->nf; ++f)
      {
        Field_T *field_rec = solve->field[f];
        Field_T *field_src = patch->pool[Ind(field_rec->name)];
        unsigned i;
        
        free_coeffs(field_rec);
        for (i = 0; i < patch->nn; ++i)
          field_rec->v[i] = field_src->v[i];
      }
      
    }/* end of for (cn = 0; cn < solution_man->ns; ++cn) */
    
  }/* end of FOR_ALL_PATCHES(p,grid) */
}

/* solving a.x = b.
// THREAD SAFE.
// ->return value: EXIT_SUCCESS
*/
static int solve_ax_b_ppm(Patch_T *const patch,const unsigned cn)
{
  Solve_T *const solve = patch->solution_man->solve[cn];
  Matrix_T *m = 0;
  
  if (strstr_i(GetParameterS("Solving_Preconditioning"),"yes"))
    precondition(solve->a,solve->b);
  if (strcmp_i(GetParameterS_E("Linear_Solver"),"UMFPACK"))
  {
    UmfPack_T umfpack[1] = {0};
    if (!solve->a->ccs_f)
    {
      //test
      //pr_matrix(solve->a,"J_ij");
      //end
      m = cast_matrix_ccs(solve->a);
      free_matrix(solve->a);
      solve->a = m;
    }
    //test
    /*Node_T **node = patch->node;
    for(unsigned nn=0;nn<patch->nn;++nn)
      fprintf(stderr,"F(%f,%f,%f|%u)=%f\n",
      node[nn]->x[0],node[nn]->x[1],node[nn]->x[2],nn,
      solve->b[nn]);
      */
    //end
    umfpack->a = solve->a;
    umfpack->b = solve->b;
    umfpack->x = solve->x;
    direct_solver_umfpack_di(umfpack);
    free_matrix(m);
  }
  else if (strcmp_i(GetParameterS_E("Linear_Solver"),"UMFPACK_long"))
  {
    UmfPack_T umfpack[1] = {0};
    if (!solve->a->ccs_l_f)
    {
      abortEr(INCOMPLETE_FUNC);
      //m = cast_matrix_ccs(solve->a);
      //free_matrix(solve->a);
      //solve->a = m;
    }
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
// for the given collection of fields cn according to the equations for a.
// THREAD SAFE.
// ->return value: EXIT_SUCCESS
*/
static int a_in_ax_b_bndry_ppm(Patch_T *const patch,const unsigned cn)
{
  const unsigned nintfc = countf(patch->interface);
  Solve_T *const solve = patch->solution_man->solve[cn];
  unsigned intfc;
  
  /* loop over all interfaces */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    unsigned nsfc = interface->ns;
    unsigned sfc;
    
    /* loop over all subfaces */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      Boundary_Condition_T BC;
      
      BC.subface = subface;
      BC.solve   = solve;
      BC.cn      = cn;
      
      if (!subface->exterF)/* if subface is internal */
      {
        abortEr(INCOMPLETE_FUNC);
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
// for the given collection of fields cn according to the equations for a.
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
// for the given collection of fields cn according to field equations.
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
    
    /* fill b in a.x = b for the specified field */
    field_eq(field,b);
  }
  
  return EXIT_SUCCESS;
}

/* filling b in a.x = b for boundary points of the given patch
// for the given collection of fields cn according to the boundary 
// conditions between adjacent patches and boundary conditions equations.
// THREAD SAFE.
// ->return value: EXIT_SUCCESS
*/
static int b_in_ax_b_bndry_ppm(Patch_T *const patch,const unsigned cn)
{
  const unsigned nintfc = countf(patch->interface);
  Solve_T *const solve = patch->solution_man->solve[cn];
  unsigned intfc;
  
  /* loop over all interfaces */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    unsigned nsfc = interface->ns;
    unsigned sfc;
    
    /* loop over all subfaces */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      Boundary_Condition_T BC;
      
      BC.subface = subface;
      BC.solve   = solve;
      BC.cn      = cn;
      
      if (!subface->exterF)/* if subface is internal */
      {
        abortEr(INCOMPLETE_FUNC);
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

/* copy values of fields or their derivatives
// along the normal vector at collocated points occured at
// the boundary between two patches for parallel patch method.
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
  unsigned xyz1,xyz2,boundary;
  unsigned n;
  unsigned i;
  
  /* df/dn = df/dn|adjacent */
  if (subface->df_dn&&0)
  {
    Patch_T patch_tmp = make_temp_patch(patch_adj);/* for thread safety purposes */
    double *Nvec;/* normal vector */
    const char *der0 = "x",*der1 = "y",*der2 = "z";
    
    /* if one wants to solve the whole equations on curvilinear patches */
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
      Field_T *field_tmp = add_field("tmp_field","(3dim)",&patch_tmp,NO);
      double *f_a = 0,*f_b = 0,*f_c = 0;
      double *f_a_adj = 0,*f_b_adj = 0,*f_c_adj = 0;
      field     = solve->field[i];
      field_adj = patch_adj->solution_man->solve[bc->cn]->field[i];
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
        xyz1 = id[n];
        xyz2 = adjid[n];
        boundary = id[n];
        
        b[boundary] = Nvec[0]*(f_a[xyz1] - f_a_adj[xyz2]) +
                      Nvec[1]*(f_b[xyz1] - f_b_adj[xyz2]) +
                      Nvec[2]*(f_c[xyz1] - f_c_adj[xyz2]) ;
          //test
          
          //fprintf(stderr,"b[%u]=%g\n",boundary,b[boundary]);
          //if (EQL(b[boundary],4))
          {
            //double *pp = subface->patch->node[xyz1]->x;
            //fprintf(stderr,"at(%f,%f,%f) f_b=%g,f_b_adj=%g\n",
            //pp[0],pp[1],pp[2],f_b[xyz1],f_b_adj[xyz2]);
          }//*/
          //end
      }
      /* freeing memories */
      field_tmp ->v = 0;/* since pointing to field_adj->v */
      remove_field(field_tmp);
      free_temp_patch(&patch_tmp);
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
        xyz1 = id[n];
        xyz2 = adjid[n];
        boundary = id[n];
        
        b[boundary] = field->v[xyz1] - field_adj->v[xyz2];
      }
      
    }/* end of for (i = 0; i < solve->nf; ++i) */
  }
  
  return EXIT_SUCCESS;
}


/* making those parts of jacobian related to the section in b in which 
// it copys values of fields or their derivatives
// along the normal vector between collocated points 
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
  if (subface->df_dn&&0)
  {
    fJs_T *j_x = 0,*j_y = 0,*j_z = 0;
    Matrix_T *j0 = 0,*j1 = 0,*j2 = 0;
    Point_T point;
    point.patch = patch;
    point.face  = subface->face;
    
    /* if one wants to solve whole equations on curvilinear patches */
    if (!strcmp_i(GetParameterS("Solving_Interpolation_Normal"),"Cartesian_Normal"))
    {
      abortEr(INCOMPLETE_FUNC);
      
      const char *types[] = {"j_a","j_b","j_c",0};
      prepare_Js_jacobian_eq(patch,types);
      j0  = get_j_matrix(patch,"j_a");
      j1  = get_j_matrix(patch,"j_b");
      j2  = get_j_matrix(patch,"j_c");
      j_x = get_j_reader(j0);
      j_y = get_j_reader(j1);
      j_z = get_j_reader(j2);
  
    }
    else
    {
      const char *types[] = {"j_x","j_y","j_z",0};
      prepare_Js_jacobian_eq(patch,types);
      j0  = get_j_matrix(patch,"j_x");
      j1  = get_j_matrix(patch,"j_y");
      j2  = get_j_matrix(patch,"j_z");
      j_x = get_j_reader(j0);
      j_y = get_j_reader(j1);
      j_z = get_j_reader(j2);
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
              Nvec[0]*j_x(j0,bndry[ijk],lmn-initial) +
              Nvec[1]*j_y(j1,bndry[ijk],lmn-initial) +
              Nvec[2]*j_z(j2,bndry[ijk],lmn-initial) ;
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
  if (subface->df_dn&&0)
  {
    fJs_T *j_x = 0,*j_y = 0,*j_z = 0;
    Matrix_T *j0 = 0,*j1 = 0,*j2 = 0;
    Point_T point;
    point.patch = patch;
    point.face  = subface->face;
    
    /* if one wants to solve whole equations on curvilinear patches */
    if (!strcmp_i(GetParameterS("Solving_Interpolation_Normal"),"Cartesian_Normal"))
    {
      abortEr(INCOMPLETE_FUNC);
      
      const char *types[] = {"j_a","j_b","j_c",0};
      prepare_Js_jacobian_eq(patch,types);
      j0  = get_j_matrix(patch,"j_a");
      j1  = get_j_matrix(patch,"j_b");
      j2  = get_j_matrix(patch,"j_c");
      j_x = get_j_reader(j0);
      j_y = get_j_reader(j1);
      j_z = get_j_reader(j2);
  
    }
    else
    {
      const char *types[] = {"j_x","j_y","j_z",0};
      prepare_Js_jacobian_eq(patch,types);
      j0  = get_j_matrix(patch,"j_x");
      j1  = get_j_matrix(patch,"j_y");
      j2  = get_j_matrix(patch,"j_z");
      j_x = get_j_reader(j0);
      j_y = get_j_reader(j1);
      j_z = get_j_reader(j2);
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
              Nvec[0]*j_x(j0,bndry[ijk],lmn-initial) +
              Nvec[1]*j_y(j1,bndry[ijk],lmn-initial) +
              Nvec[2]*j_z(j2,bndry[ijk],lmn-initial) ;
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
  unsigned boundary;
  double X[3],*x;
  unsigned i;
  
  /* df/dn = df/dn|adjacent */
  if (subface->df_dn&&0)
  {
    double *Nvec;/* normal vector */
    const char *der0 = "x",*der1 = "y",*der2 = "z";
    Patch_T adj_patch_tmp = make_temp_patch(adj_patch);/* for thread safety purposes */
    
    /* if one wants to solve the whole equations on curvilinear patches */
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
      double *b = &solve->b[solve->f_occupy[i]];
      Interpolation_T *interp_s_a_adj = init_interpolation();
      Interpolation_T *interp_s_b_adj = init_interpolation();
      Interpolation_T *interp_s_c_adj = init_interpolation();
      Field_T *field     = solve->field[i];
      Field_T *field_adj = adj_patch->solution_man->solve[bc->cn]->field[i];
      Field_T *field_tmp = add_field("tmp_field","(3dim)",&adj_patch_tmp,NO);
      Field_T *f_a_adj   = add_field("f_a_adj","(3dim)",&adj_patch_tmp,NO),
              *f_b_adj   = add_field("f_b_adj","(3dim)",&adj_patch_tmp,NO), 
              *f_c_adj   = add_field("f_c_adj","(3dim)",&adj_patch_tmp,NO);
      double *f_a = 0, *f_b = 0, *f_c = 0;
      unsigned n;
      
      field_tmp->v = field_adj->v;/* soft copying values of fields */
      f_a_adj->v = Partial_Derivative(field_tmp,der0);
      f_b_adj->v = Partial_Derivative(field_tmp,der1);
      f_c_adj->v = Partial_Derivative(field_tmp,der2);
      
      f_a = Partial_Derivative(field,der0);
      f_b = Partial_Derivative(field,der1);
      f_c = Partial_Derivative(field,der2);
      
      interp_s_a_adj->field = f_a_adj;
      interp_s_b_adj->field = f_b_adj;
      interp_s_c_adj->field = f_c_adj;
      
      fill_interpolation_flags(interp_s_a_adj,subface);
      plan_interpolation(interp_s_a_adj);
      fill_interpolation_flags(interp_s_b_adj,subface);
      plan_interpolation(interp_s_b_adj);
      fill_interpolation_flags(interp_s_c_adj,subface);
      plan_interpolation(interp_s_c_adj);
      
      for (n = 0; n < np; ++n)
      {
        Point_T point;
        point.ind   = id[n];
        point.patch = subface->patch;
        point.face  = subface->face;
        Nvec = normal_vec(&point);
        x = subface->patch->node[id[n]]->x;
        X_of_x(X,x,&adj_patch_tmp);
        interp_s_a_adj->X = X[0];
        interp_s_a_adj->Y = X[1];
        interp_s_a_adj->Z = X[2];
        interp_s_b_adj->X = X[0];
        interp_s_b_adj->Y = X[1];
        interp_s_b_adj->Z = X[2];
        interp_s_c_adj->X = X[0];
        interp_s_c_adj->Y = X[1];
        interp_s_c_adj->Z = X[2];
        
        boundary = id[n];
        
        b[boundary] = Nvec[0]*(f_a[boundary] - execute_interpolation(interp_s_a_adj)) +
                      Nvec[1]*(f_b[boundary] - execute_interpolation(interp_s_b_adj)) +
                      Nvec[2]*(f_c[boundary] - execute_interpolation(interp_s_c_adj));
      }
      /* freeing */
      field_tmp->v = 0;
      free_interpolation(interp_s_a_adj);
      free_interpolation(interp_s_b_adj);
      free_interpolation(interp_s_c_adj);
      remove_field(f_a_adj);
      remove_field(f_b_adj);
      remove_field(f_c_adj);
      remove_field(field_tmp);
      free_temp_patch(&adj_patch_tmp);
      free(f_a);
      free(f_b);
      free(f_c);
    }/* end of for (i = 0; i < solve->nf; ++i) */
  }/* end of if (subface->df_dn) */
  /* f = f|adjacent */
  else
  {
    Patch_T adj_patch_tmp = make_temp_patch(adj_patch);/* for thread safety purposes */
    
    for (i = 0; i < solve->nf; ++i)
    {
      double *b = &solve->b[solve->f_occupy[i]];
      Interpolation_T *interp_s = init_interpolation();
      Field_T *field     = solve->field[i];
      Field_T *field_adj = adj_patch->solution_man->solve[bc->cn]->field[i];
      Field_T *field_tmp = add_field("tmp_field","(3dim)",&adj_patch_tmp,NO);
      unsigned n;
      
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
        boundary = id[n];
        
        b[boundary] = field->v[boundary] - execute_interpolation(interp_s);
      }
      /* freeing */
      field_tmp->v = 0;
      remove_field(field_tmp);
      free_temp_patch(&adj_patch_tmp);
      free_interpolation(interp_s);
    }/* end of for (i = 0; i < solve->nf; ++i) */
  }/* end of else */
  
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