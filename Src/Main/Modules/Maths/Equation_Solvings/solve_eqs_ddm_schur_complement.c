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
   
   1. solve BE'  = E and Bf'= f where E' = B^(-1)*E and f' = B^(-1)*f
   2. compute g' = g - Ff'
   3. compute S  = C - FE'
   4. solve Sy = g'
   5. compute x = f'-E'y
  
// some definitions used here:
// outer boundary = boundary at the outer most part of the patch, like infinity
// inner boundary = boundary at some internal part of the patch, like black hole horizon.
// boundary points = points (nodes) on outer boundary + inner boundary
// subdomain = inner points + boundary points
// interface = all nodes(points) - subdomain.
// note: interface struct includes all boundary points.
*/

/* using Schur Complement domain decomposition method
// to solve equation. This method is using direct solver
// ,like UMFPACK solver, in parallel fashion to solve 
// the whole systme of equations Ax=B; it's fast and robust.
*/
int ddm_schur_complement(Solve_Equations_T *const SolveEqs)
{
  FUNC_TIC
  
  Grid_T *grid;
  char **field_name = 0;/* name of all fields (equations) to be solved */
  Uint nf = 0;/* number of all fields */
  Uint f;/* dummy index */
  
  pr_intro_ddm_schur_complement();
  
  /* read order of fields to be solved from input */
  field_name = get_solving_field_name(SolveEqs->solving_order,&nf);
  
  /* synchronize pools with the default grid */
  sync_patch_pools(SolveEqs->grid,SolveEqs);
  
  /* solving fields in order */
  for (f = 0; f < nf; ++f)
  {
    pr_half_line_custom('-');
    printf("{ Solving equation for field: \"%s\" ...\n",field_name[f]);
    pr_half_line_custom('-');
    
    /* set the name of the field we are solving it */
    SolveEqs->field_name = field_name[f];
    
    /* get the computational grid */
    grid = get_grid_solve_equations(SolveEqs);
    
    /* set solving_man->cf */
    set_solving_man_cf(SolveEqs);
    
    /* set solving_man->settings */
    set_solving_man_settings(SolveEqs);
    
    /* picking up labeling, mapping etc. */
    preparing_ingredients(SolveEqs);
    
    /* solve field[f] on grid_target */
    solve_field(SolveEqs);
    
    /* updating source if any has been set */
    if (SolveEqs->SourceUpdate)
      SolveEqs->SourceUpdate(grid,field_name[f]);
    
    sync_patch_pools(grid,SolveEqs);
      
    printf("\n");
    pr_half_line_custom('-');
    printf("} Solving equation for field: \"%s\" ==> Done.\n",field_name[f]);
    pr_half_line_custom('-');
  }
  
  /* free names */
  free_2d_mem(field_name,nf);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}  

/* solving the field.
// ->return value: EXIT_SUCCESS. */
static int solve_field(Solve_Equations_T *const SolveEqs)
{
  Grid_T *const grid = get_grid_solve_equations(SolveEqs);
  int CONTINUE = 1;
  int step = 0;
  
  /* if number grid only has one patch */
  if (grid->np == 1)
  {
    while (CONTINUE)
    {
      Patch_T *patch = grid->patch[0];
      DDM_Schur_Complement_T *Schur = patch->solving_man->method->SchurC;
      double time1 = get_time_sec();
      
      /* set current step */
      patch->solving_man->settings->solver_step = step;
      
      make_f(patch);/* making f */
      
      if (!step)/* only at first step */
        set_solving_man_settings_Frms_i_single_patch(patch);
        
      /* calculate the current residual and set it in patch->solving_man->Frms */  
      calculate_residual_single_patch(patch);
      
      /* check the stop criteria */
      CONTINUE = SolveEqs->StopCriteria(grid,SolveEqs->field_name);
      if (!CONTINUE)
      {
        free(Schur->f);/* free{f} */
        free_solving_man_settings(grid);/* free (solving_man->settings) */
        break;
      }

      pr_line_custom('~');
      printf("{---> %s equation:\n",SolveEqs->field_name);
      printf("      |---> Newton step '%d':\n",step+1);
      pr_line_custom('~');
      
      making_B_single_patch(patch);/* making B */
      solve_Bx_f(patch);/* solve Bx=f, free{B,f} */
      update_field_single_patch(patch);
      free(Schur->x);/* free{x} */
      /* updating fields and their derivative and related */
      if (SolveEqs->FieldUpdate)/* if any FieldUpdate set */
        SolveEqs->FieldUpdate(patch,SolveEqs->field_name);
      
      pr_line_custom('~');
      printf("}---> %s equation:\n",SolveEqs->field_name);
      printf("      |---> Newton step '%d' is done.\n",step+1);
      printf("            |---> Elapsed seconds = %.0f .\n",get_time_sec()-time1);
      pr_clock();
      pr_line_custom('~');
      step++;
    }
  }
  else/* multi-domains grid */
  {
    while (CONTINUE)
    {
      const Uint npatch = grid->np;
      double *g_prime = 0;
      Matrix_T *S;
      double time1 = get_time_sec();
      double time_temp;
      char *pr_msg = 0;
      Uint p;
      
      /* set current step */
      set_solving_man_settings_solver_step(grid,step);
      
      printf("{ Compute f and g ...\n");
      time_temp = get_time_sec();
      DDM_SCHUR_OpenMP(omp parallel for)
      for (p = 0; p < npatch; ++p)
      {
        Patch_T *patch = grid->patch[p];
        make_f(patch);
        make_partial_g(patch);
      }
      make_g(grid);/* free pg */
      printf("} Compute f and g --> Done. ( Wall-Clock = %.0fs )\n",
              get_time_sec()-time_temp);
      fflush(stdout);
      if (!step)/* only at first step */
        set_solving_man_settings_Frms_i(grid);
      
      /* calculate the current residual and set it in patch->solving_man->Frms */  
      calculate_residual(grid);
      
      /* check the stop criteria */
      CONTINUE = SolveEqs->StopCriteria(grid,SolveEqs->field_name);
      if (!CONTINUE)
      {
        free_schur_f_g(grid);/* free {f,g} */
        free_solving_man_settings(grid);/* free (solving_man->settings) */
        break;
      }
      
      pr_line_custom('~');
      printf("{---> %s equation:\n",SolveEqs->field_name);
      printf("      |---> Newton step '%d':\n",step+1);
      pr_line_custom('~');
      
      DDM_SCHUR_OpenMP_SET_NUM_THREADS(omp parallel for)
      for (p = 0; p < npatch; ++p)
      {
        Patch_T *patch = grid->patch[p];
        double tic = get_time_sec();
        char *msg  = calloc(MSG_SIZE2,1); 
        IsNull(msg);
        
        char *msg1 = making_B_and_E(patch);
        char *msg2 = making_E_prime_and_f_prime(patch);/* free{B,E} */
        char *msg3 = making_F_and_C(patch);
        char *msg4 = making_F_by_f_prime(patch);
        char *msg5 = making_F_by_E_prime(patch);/* free {F} */
        char *msg6 = making_sub_S_matrix(patch);/* free C and 
                                                // F_by_E_prime
                                                // and save subS */
        
        sprintf(msg,"{ %s ...\n\n"
                    "%s\n"
                    "%s\n"
                    "%s\n"
                    "%s\n"
                    "%s\n"
                    "%s\n"
                    "} %s --> Done. ( Wall-Clock = %.0fs )\n",
                patch->name,
                msg1,msg2,msg3,msg4,msg5,msg6,
                patch->name,
                get_time_sec()-tic);
                
        printf("%s\n",msg);
        fflush(stdout);
        
        free(msg1);
        free(msg2);
        free(msg3);
        free(msg4);
        free(msg5);
        free(msg6);
        free(msg);
      }
      DDM_SCHUR_OpenMP_UNSET_NUM_THREADS
      
      /* compute g' */
      printf("{ Compute g' ...\n");
      time_temp = get_time_sec();
      g_prime = compute_g_prime(grid);/* free {g,Ff'} */
      printf("} Compute g' --> Done. ( Wall-Clock = %.0fs )\n\n",
              get_time_sec()-time_temp);
      fflush(stdout);
      
      /* compute S */
      printf("{ Compute S ...\n");
      time_temp = get_time_sec();
      S = compute_S(grid);/* free {subS} */
      printf("} Compute S --> Done. ( Wall-Clock = %.0fs )\n\n",get_time_sec()-time_temp);
      fflush(stdout);
      
      /* solve Sy = g' */
      printf("{ Solve S y = g' ...\n");
      time_temp = get_time_sec();
      pr_msg = solve_Sy_g_prime(S,g_prime,grid);/* free{S,g'} */
      printf("%s",pr_msg);
      printf("} Solve S y = g' --> Done. ( Wall-Clock = %.0fs )\n\n",get_time_sec()-time_temp);
      fflush(stdout);
      free(pr_msg);
      
      DDM_SCHUR_OpenMP(omp parallel for)
      for (p = 0; p < npatch; ++p)
      {
        Patch_T *patch    = grid->patch[p];
        Solving_Man_T *SM = patch->solving_man;
        char *field_name  = SM->field_name[SM->cf];
        char msg1[MSG_SIZE1],msg2[MSG_SIZE1],msg3[MSG_SIZE1],
            msg4[MSG_SIZE1],msg5[MSG_SIZE1];
        double tic;
        double w_lambda = patch->solving_man->settings->relaxation_factor;
        
        sprintf(msg1,"{ %s ...\n",patch->name);
        sprintf(msg4,"} %s --> Done.\n",patch->name);
        
        /* x = f'-E'y */
        tic = get_time_sec();
        compute_x(patch);
        
        sprintf(msg2,"  --> Compute x. ( Wall-Clock = %.0fs )\n",
                get_time_sec()-tic);
        
        free_E_Trans_prime(patch);/* free{E^T'} */
        
        /* update field */
        tic = get_time_sec();
        update_field(patch);
        sprintf(msg3,"  --> Update '%s', weight = %g. ( Wall-Clock = %.0fs )\n",
                field_name,w_lambda,get_time_sec()-tic);
        
        free_x(patch);/* free{x} */
        
        msg5[0] = '\0';
        /* updating fields and their derivatives and related */
        if (SolveEqs->FieldUpdate)/* if defined any FieldUpdate function */
        {
          tic = get_time_sec();
          SolveEqs->FieldUpdate(patch,field_name);
          sprintf(msg5,"  --> Update '%s' derivatives. ( Wall-Clock = %.0fs )\n",
                  field_name,get_time_sec()-tic);
        }
        
        printf("%s%s%s%s%s\n",msg1,msg2,msg3,msg5,msg4);
        fflush(stdout);
      }
      free_y(grid);/* free{y} */
      
      pr_line_custom('~');
      printf("}---> %s equation:\n",SolveEqs->field_name);
      printf("      |---> Newton step '%d' is done .\n",step+1);
      printf("            |---> Wall-Clock = %.0fs .\n",get_time_sec()-time1);
      pr_clock();
      pr_line_custom('~');
      
      step++;
    }/* end of while (CONTINUE) */
  }/* end of else */
  
  return EXIT_SUCCESS;
}

/* finding x and y now time to update field to reduce residual */
static void update_field(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  const Uint NS = Schur->NS;
  const Uint NI = Schur->NI;
  const Uint *const inv = Schur->inv;
  const Uint *const Iinv = Schur->Iinv;
  const double *const y = Schur->y;
  const double *const x = Schur->x;
  const Uint cf = patch->solving_man->cf;
  const char *const field_aliased = patch->solving_man->field_aliased[cf];  
  Field_T *const f = patch->fields[Ind(field_aliased)];
  const double lambda = patch->solving_man->settings->relaxation_factor;
  const double *const u_old = f->v;
  double *const u_new = f->v;
  double *const u_bckup = alloc_double(patch->nn);
  Uint s,s_node,i,i_node;
  
  free_coeffs(f);
  for (s = 0; s < NS; ++s)
  {
    s_node = inv[s];
    u_bckup[s_node] = u_old[s_node];
    u_new[s_node]   = u_old[s_node]-lambda*x[s];
  }
  for (i = 0; i < NI; ++i)
  {
    i_node = Iinv[i];
    u_bckup[i_node] = u_old[i_node];
    u_new[i_node]   = u_old[i_node]-lambda*y[i];
  }
  Free(patch->solving_man->settings->last_sol);
  patch->solving_man->settings->last_sol = 0;
  patch->solving_man->settings->last_sol = u_bckup;
  
}

/* updating field when the grid only has single patch */
static void update_field_single_patch(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  const Uint NS = Schur->NS;
  const Uint *const inv = Schur->inv;
  const double *const x = Schur->x;
  const Uint cf = patch->solving_man->cf;
  const char *const field_aliased = patch->solving_man->field_aliased[cf];
  Field_T *const f = patch->fields[Ind(field_aliased)];
  const double lambda = patch->solving_man->settings->relaxation_factor;
  const double *const u_old = f->v;
  double *const u_new = f->v;
  double *const u_bckup = alloc_double(patch->nn);
  Uint s,s_node;
  
  free_coeffs(f);
  for (s = 0; s < NS; ++s)
  {
    s_node = inv[s];
    u_bckup[s_node] = u_old[s_node];
    u_new[s_node]   = u_old[s_node]-lambda*x[s];
  }
  Free(patch->solving_man->settings->last_sol);
  patch->solving_man->settings->last_sol = 0;
  patch->solving_man->settings->last_sol = u_bckup;
}

/* free x in Schur */
static void free_x(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const Schur = 
                          patch->solving_man->method->SchurC;
  free(Schur->x);
}

/* free y in Schur.
// NOTE: since y memory is shared between the patches 
// ONLY free one of them. */
static void free_y(Grid_T *const grid)
{
  DDM_Schur_Complement_T *const Schur = 
                          grid->patch[0]->solving_man->method->SchurC;
  
  free(Schur->y); 
}

/* free E_Trans_prime matrix in each patch. */
static void free_E_Trans_prime(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  if (Schur->E_Trans_prime)
    free_matrix(Schur->E_Trans_prime);
}

/* x = f'-E'y */
static void compute_x(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  const Uint NS = Schur->NS;
  double *f_prime  = Schur->f_prime;
  const double *const y = Schur->y;
  double *const x = alloc_double(NS);
  double Ey;/* E_Trans_prime by y */

# if E_Trans_prime_IS_RMO_FORMAT
  const Uint NI = Schur->NI;
  double *const E_Trans_prime = Schur->E_Trans_prime->rmo->A;
  const Uint ETp_ncol         = (Uint)Schur->E_Trans_prime->col;
  Uint i,s;
  
  for(s = 0; s < NS; ++s)
  {
    Ey = 0;/* E_Trans_prime by y */
    for (i = 0; i < NI; ++i)
      Ey += E_Trans_prime[i_j_to_ij(ETp_ncol,i,s)]*y[i];
    x[s] = f_prime[s]-Ey;
  }

# elif E_Trans_prime_IS_CCS_FORMAT
  double *const E_Ax = Schur->E_Trans_prime->ccs->Ax;
  int *const E_Ap = Schur->E_Trans_prime->ccs->Ap;
  int *const E_Ai = Schur->E_Trans_prime->ccs->Ai;
  const int Nc    = (int)Schur->E_Trans_prime->col;/* == NS */
  
  for(int c = 0; c < Nc; ++c)
  {
    Ey = 0;/* E_Trans_prime by y */
    for (int i = E_Ap[c]; i < E_Ap[c+1]; ++i)
      Ey += E_Ax[i]*y[E_Ai[i]];
    x[c] = f_prime[c]-Ey;
  }

# endif

  free(Schur->f_prime);
  Schur->x = x;
}

/* allocate y and solve Sy = g' and then populate SchurC->y;
// note: it frees S and g_prime too. */
static char *solve_Sy_g_prime(Matrix_T *const S,double *const g_prime,Grid_T *const grid)
{
  /* since NI_total and NI_p are unique we pick one of them from patch[0] */
  const Uint NI_total    = 
                  grid->patch[0]->solving_man->method->SchurC->NI_total;
  const Uint *const NI_p = 
                  grid->patch[0]->solving_man->method->SchurC->NI_p;
  double *y = alloc_double(NI_total);
  Umfpack_T *umfpack = init_umfpack();
  DDM_Schur_Complement_T *Schur;
  Uint R = 0;
  Uint p;
  char *msg = calloc(MSG_SIZE2,1);
  IsNull(msg);
  
  umfpack->a = S;
  umfpack->b = g_prime;
  umfpack->x = y;
  umfpack->refinement = 
    grid->patch[0]->solving_man->settings->umfpack_refine;
  
  if (S->ccs_f)
    direct_solver_umfpack_di(umfpack);
  
  else if (S->ccs_l_f)
    direct_solver_umfpack_dl(umfpack);
    
  else
    Error0(NO_OPTION);

  /* populate SchurC->y */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Schur = patch->solving_man->method->SchurC;
    Schur->y = &y[R];
    R += NI_p[p];
  }
  
  sprintf(msg,"%s",umfpack->description);
  
  free_matrix(S);
  free(g_prime);
  free_umfpack(umfpack);
  
  return msg;
}

/* allocate and compute S in either ccs or ccs long format.
// ->return value: S matrix. */
static Matrix_T *compute_S(Grid_T *const grid)
{
  if (0)
    return compute_S_CCS(grid);
  
  if (1)/* if the S matrix gets very hug then use this */
    return compute_S_CCS_long(grid);
  
}

/* ->: some info.
// having FE' and C compute S = C - FE' which is Schur->subS. */
static char *making_sub_S_matrix(Patch_T *const patch)
{
  const double time1 = get_time_sec();
  DDM_Schur_Complement_T *const Schur = 
              patch->solving_man->method->SchurC;
  const Uint Np          = Schur->np;
  const Uint *const NI_p = Schur->NI_p;
  const Uint NI          = Schur->NI;
  const Uint NI_total    = Schur->NI_total;
  Matrix_T ** FxE  = Schur->F_by_E_prime;
  Matrix_T ** C    = Schur->C;
  Matrix_T **stack = calloc(Np,sizeof(*stack)); IsNull(stack);
  double *a  = 0;
  double **b = 0;
  double **c = 0;
  long a_row,a_col,b_row,b_col;
  long row,col;
  Uint p;  
  char *msg = calloc(MSG_SIZE1,1);
  IsNull(msg);
  
  for (p = 0; p < Np; ++p)
  {
    if (FxE[p]&&C[p])/* if both FE' and C are none empty */
    {
      a_row = FxE[p]->row;
      a_col = FxE[p]->col;
      b_row = C[p]->row;
      b_col = C[p]->col;
      /* some checks */
      assert(a_row == b_row);
      assert(a_col == b_col);
      assert(FxE[p]->rmo_f && C[p]->reg_f);
      
      stack[p] = alloc_matrix(REG_SF,a_row,a_col);
      a = FxE[p]->rmo->A;
      b = C[p]->reg->A;
      c = stack[p]->reg->A;
      
      for (row = 0; row < a_row; ++row)
        for (col = 0; col < a_col; ++col)
          c[row][col] = b[row][col]-a[i_j_to_ij(a_col,row,col)];
      
      free_matrix(C[p]);
      free_matrix(FxE[p]);
    }
    else if (FxE[p])
    {
      a_row = FxE[p]->row;
      a_col = FxE[p]->col;
      assert(FxE[p]->rmo_f);
      stack[p] = FxE[p];
      a = FxE[p]->rmo->A;
      c = stack[p]->reg->A;
      
      for (row = 0; row < a_row; ++row)
        for (col = 0; col < a_col; ++col)
          c[row][col] = -a[i_j_to_ij(a_col,row,col)];
      
      FxE[p] = 0;
    }
    else if (C[p])
    {
      b_row = C[p]->row;
      b_col = C[p]->col;
      assert(C[p]->reg_f);
      stack[p] = C[p];
      C[p] = 0;
    }
    
  }
  
  Schur->subS = compress_stack2ccs(stack,Np,NI_p,NI_total,NI,YES);
  
  /* free */
  Free(C);
  Free(FxE);
  Schur->F_by_E_prime = 0;
  Schur->C = 0;
  Free(stack);
  
  sprintf(msg,"{ Make subS ...\n"
              "} Make subS --> Done. ( Wall-Clock = %.0fs )\n",
              get_time_sec()-time1);
  
  return msg;
}

/* allocate and compute S ccs long format storage
// ->return value: S matrix. */
static Matrix_T *compute_S_CCS_long(Grid_T *const grid)
{
  /* since NI_total and NI_p are unique we pick one of them from patch[0] */
  const Uint NI_total    = 
                  grid->patch[0]->solving_man->method->SchurC->NI_total;
  const Uint npatch = grid->np;
  Matrix_T *S = alloc_matrix(CCS_L_SF,NI_total,NI_total);
  long *Ap = 0;
  long *Ai = 0;
  double *Ax = 0;
  long i,j,nnz,R;
  Uint p;
  
  /* to be safe we used long format data type */
  Ap = calloc(NI_total+1,sizeof(*Ap));
  IsNull(Ap);
  R = nnz = 0;
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    DDM_Schur_Complement_T *Schur = patch->solving_man->method->SchurC;
    Matrix_T *subS = Schur->subS;
    int *Ap1    = subS->ccs->Ap;
    int *Ai1    = subS->ccs->Ai;
    double *Ax1 = subS->ccs->Ax;
    long Nc1    = subS->col;
    
    Ai = realloc(Ai,(long Uint)(Ap[R]+Ap1[Nc1])*sizeof(*Ai));
    IsNull(Ai);
    Ax = realloc(Ax,(long Uint)(Ap[R]+Ap1[Nc1])*sizeof(*Ax));
    IsNull(Ax);
    
    for (i = 0; i < Nc1; ++i)
    {
      Ap[i+R] = Ap1[i]+nnz;
      for (j = Ap1[i]; j < Ap1[i+1]; ++j)
      {
        Ai[nnz+j] = Ai1[j];
        Ax[nnz+j] = Ax1[j];
      }
    }
    R    += Nc1;
    nnz  += Ap1[Nc1];
    Ap[R] = nnz;
    
    free_matrix(subS);
  }
  
  S->ccs_long->Ap = Ap;
  S->ccs_long->Ai = Ai;
  S->ccs_long->Ax = Ax;
  
  return S;
}

/* allocate and compute S ccs format storage
// ->return value: S matrix. */
static Matrix_T *compute_S_CCS(Grid_T *const grid)
{
  /* since NI_total and NI_p are unique we pick one of them from patch[0] */
  const Uint NI_total    = 
                  grid->patch[0]->solving_man->method->SchurC->NI_total;
  const Uint npatch = grid->np;
  Matrix_T *S = alloc_matrix(CCS_SF,NI_total,NI_total);
  int *Ap = 0;
  int *Ai = 0;
  double *Ax = 0;
  int nnz;
  long i,j,R;
  Uint p;
  
  Ap = calloc(NI_total+1,sizeof(*Ap));
  IsNull(Ap);
  R = nnz = 0;
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    DDM_Schur_Complement_T *Schur = patch->solving_man->method->SchurC;
    Matrix_T *subS = Schur->subS;
    int *Ap1    = subS->ccs->Ap;
    int *Ai1    = subS->ccs->Ai;
    double *Ax1 = subS->ccs->Ax;
    long Nc1    = subS->col;
    
    Ai = realloc(Ai,(long Uint)(Ap[R]+Ap1[Nc1])*sizeof(*Ai));
    IsNull(Ai);
    Ax = realloc(Ax,(long Uint)(Ap[R]+Ap1[Nc1])*sizeof(*Ax));
    IsNull(Ax);
    
    for (i = 0; i < Nc1; ++i)
    {
      Ap[i+R] = Ap1[i]+nnz;
      for (j = Ap1[i]; j < Ap1[i+1]; ++j)
      {
        Ai[nnz+j] = Ai1[j];
        Ax[nnz+j] = Ax1[j];
      }
    }
    R    += Nc1;
    nnz  += Ap1[Nc1];
    Ap[R] = nnz;
    
    free_matrix(subS);
  }
  
  S->ccs->Ap = Ap;
  S->ccs->Ai = Ai;
  S->ccs->Ax = Ax;
  
  return S;
}

/* compute g'.
// ->return value : g'. */
static double *compute_g_prime(Grid_T *const grid)
{
  /* since all of SchurC->F_by_f_prime's have same dimensions we pick one
  // of them thus we not only don't need to alloc memory but also we skip one of
  // the calculations and then start from patch[1]. same for NI_total. */
  double *const Ff = grid->patch[0]->solving_man->method->SchurC->F_by_f_prime;
  const Uint NI_total = grid->patch[0]->solving_man->method->SchurC->NI_total;
  DDM_Schur_Complement_T *Schur;
  double *const g_prime = Ff;/* to save some memory */
  double *g;
  double *Ffp;
  Uint NI;
  const Uint np = grid->np;
  Uint R = 0;/* reference */
  Uint p,i;
  
  for (p = 1; p < np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Schur = patch->solving_man->method->SchurC;
    
    NI = Schur->NI;
    Ffp = Schur->F_by_f_prime;
    
    for (i = 0; i < NI_total; ++i)
      Ff[i] += Ffp[i];
    
    free(Ffp);
  }
  /* since the memory handed over g_prime we don't need this pointer anymore */
  grid->patch[0]->solving_man->method->SchurC->F_by_f_prime = 0;
  R = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Schur = patch->solving_man->method->SchurC;
    
    NI = Schur->NI;
    g  = Schur->g;
    
    for (i = 0; i < NI; ++i)
      g_prime[R+i] = g[i]-g_prime[R+i];
    R += NI;
    
    free(g);
  }
  
  return g_prime;
}

/* computing matrix multiplication: F[p][j]xE'[p].
// it also frees Fs. 
// ->return value: some info */
static char *making_F_by_E_prime(Patch_T *const patch)
{
  const double time1 = get_time_sec();
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  const Uint np = Schur->np;
  Matrix_T *const E_Trans_prime = Schur->E_Trans_prime;
  
# if E_Trans_prime_IS_RMO_FORMAT  
  Matrix_T *ETp_ccs = cast_matrix_ccs(E_Trans_prime);
  
# elif E_Trans_prime_IS_CCS_FORMAT
  Matrix_T *const ETp_ccs = E_Trans_prime;
  
# endif

  Uint p;
  char *msg = calloc(MSG_SIZE1,1);
  IsNull(msg);
  
  Schur->F_by_E_prime = calloc(np,sizeof(*Schur->F_by_E_prime));
  IsNull(Schur->F_by_E_prime);
  
  for (p = 0; p < np; ++p)
  {
    Matrix_T *F     = Schur->F[p];
    Matrix_T **FxEprime = Schur->F_by_E_prime;
    
    if (F)
    {
      #if defined(MxM_GSL_BLAS)
      
        FxEprime[p] = alloc_matrix(RMO_SF,F->row,E_Trans_prime->row);
        assert(F->rmo_f);
        assert(E_Trans_prime->rmo_f);
        assert(FxEprime[p]->rmo_f);
        
        gsl_matrix_view gslF    = gsl_matrix_view_array
          (F->rmo->A, (Uint)F->row, (Uint)F->col);
        
        gsl_matrix_view gslEtp  = gsl_matrix_view_array
          (E_Trans_prime->rmo->A, (Uint)E_Trans_prime->row, 
                                  (Uint)E_Trans_prime->col);
           
        gsl_matrix_view gslFxEp = gsl_matrix_view_array
          (FxEprime[p]->rmo->A, (Uint)FxEprime[p]->row, 
                                (Uint)FxEprime[p]->col);
        
        /* Compute F[p][j]xE'[p] */
        gsl_blas_dgemm (CblasNoTrans, CblasTrans,
                  1.0, &gslF.matrix, &gslEtp.matrix,
                  0.0, &gslFxEp.matrix);
      
      #elif defined(MxM_MKL_BLAS) || defined(MxM_C_BLAS)

        FxEprime[p] = alloc_matrix(RMO_SF,F->row,E_Trans_prime->row);
        assert(F->rmo_f);
        assert(E_Trans_prime->rmo_f);
        assert(FxEprime[p]->rmo_f);
        
        /* Compute F[p][j]xE'[p]:
        // for more info about syntax and documentations 
        // see Math Kernel Library: 
        // software.intel.com/content/www/us/en/develop/documentation/
        //                    onemkl-developer-reference-c/top.html.
        // notation: for row-major order we have a_{ij} = [j+ i*lda], 
        // where lda is the leading dimension for the array. */
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                F->row, E_Trans_prime->row, F->col,
                1.0/* alpha */,
                F->rmo->A, F->col/* lda */, 
                E_Trans_prime->rmo->A, E_Trans_prime->col/* lda */, 
                0.0/* beta */, 
                FxEprime[p]->rmo->A, FxEprime[p]->col/* lda */);

      #else
      
        /* cast to ccs format */
        Matrix_T* F_ccs  = cast_matrix_ccs(F);
        
        /* only if has any non zero entries */
        if (F_ccs->ccs->Ap[F_ccs->col] && ETp_ccs->ccs->Ap[ETp_ccs->col])
        {
          FxEprime[p] = alloc_matrix(RMO_SF,F_ccs->row,ETp_ccs->row);
          matrix_by_matrix(F_ccs,ETp_ccs,FxEprime[p],"a*transpose(b)");
        }
        free_matrix(F_ccs);
        
      #endif
      
      free_matrix(F);
    }
  }
  free(Schur->F);
  
# if E_Trans_prime_IS_RMO_FORMAT
  free_matrix(ETp_ccs);
# endif

  sprintf(msg,"{ Make F*E' ...\n"
              "} Make F*E' --> Done. ( Wall-Clock = %.0fs )\n",
              get_time_sec()-time1);
  return msg;
}

/* computing:
  |F[0][p]f'[p]|
  |F[1][p]f'[p]|
  |F[2][p]f'[p]|
  |...         |

// ->return value: some info */
static char *making_F_by_f_prime(Patch_T *const patch)
{
  const double time1 = get_time_sec();
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  const Uint np = Schur->np;
  const Uint *const NI_p = Schur->NI_p;
  Uint R = 0;/* reference */
  const double *const f_prime = Schur->f_prime;
  double *const F_by_f_prime = alloc_double(Schur->NI_total);
  Uint p;
  char *msg = calloc(MSG_SIZE1,1);
  IsNull(msg);
  
  for (p = 0; p < np; ++p)
  {
    Matrix_T *F = Schur->F[p];
    
    if (F)
    {
      matrix_by_vector(F,f_prime,&F_by_f_prime[R],NOT_INITIALIZE);
    }
    R += NI_p[p];
  }
  
  Schur->F_by_f_prime = F_by_f_prime;
  
  sprintf(msg,"{ Make F*f' ...\n"
              "} Make F*f' --> Done. ( Wall-Clock = %.0fs )\n",
              get_time_sec()-time1);
  return msg;
}

/* making F and C parts. refer to the note on the very top.
// ->return value: some info */
static char *making_F_and_C(Patch_T *const patch)
{
  const double time1 = get_time_sec();
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  Sewing_T **const sewing = Schur->sewing;
  const Uint np = Schur->np;
  Uint p;
  char *msg = calloc(MSG_SIZE1,1);
  IsNull(msg);
  
  /* allocation matrices */
  Schur->F = calloc(np,sizeof(*Schur->F));
  IsNull(Schur->F);
  Schur->C = calloc(np,sizeof(*Schur->C));
  IsNull(Schur->C);
  
  /* go thru all of sewings  */
  for (p = 0; p < np; ++p)
  {
    Uint pr;
    Uint rowF,colF,rowC,colC;
    
    /* if there is no connection bewteen the patches, skip */
    if (!sewing[p])
      continue;
    
    /* allocate F and C matrices */
    rowF = sewing[p]->NI;/* must be number of interface point for patch[p] */
    rowC = sewing[p]->NI;/* must be number of interface point for patch[p] */
    colF = Schur->NS;/* must be number of subdomain point for patch */
    colC = Schur->NI;/* must be number of interface point for patch */
    
    Schur->F[p] = alloc_matrix(RMO_SF,rowF,colF);
    Schur->C[p] = alloc_matrix(REG_SF,rowC,colC);
    
    /* go thru all of pairs in each sewings */
    for (pr = 0; pr < sewing[p]->npair; ++pr)
    {
      Pair_T *const pair = sewing[p]->pair[pr];
      
      populate_F_and_C(patch,pair);/* filling F and C pertinent to this pair */
    }
  }
  
  sprintf(msg,"{ Populate F and C ...\n"
              "} Populate F and C --> Done. ( Wall-Clock = %.0fs )\n",
              get_time_sec()-time1);
  return msg;
}

/* making F and C parts. refer to the note on the very top. 
// C kept in regular format. */
static void making_F_and_C_Regular(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  Sewing_T **const sewing = Schur->sewing;
  const Uint np = Schur->np;
  Uint p;
  
  /* allocation matrices */
  Schur->F = calloc(np,sizeof(*Schur->F));
  IsNull(Schur->F);
  Schur->C = calloc(np,sizeof(*Schur->C));
  IsNull(Schur->C);
  
  /* go thru all of sewings  */
  for (p = 0; p < np; ++p)
  {
    Uint pr;
    Uint rowF,colF,rowC,colC;
    
    /* if there is no connection bewteen the patches, skip */
    if (!sewing[p])
      continue;
    
    /* allocate F and C matrices */
    rowF = sewing[p]->NI;/* must be number of interface point for patch[p] */
    rowC = sewing[p]->NI;/* must be number of interface point for patch[p] */
    colF = Schur->NS;/* must be number of subdomain point for patch */
    colC = Schur->NI;/* must be number of interface point for patch */
    
    Schur->F[p] = alloc_matrix(REG_SF,rowF,colF);
    Schur->C[p] = alloc_matrix(REG_SF,rowC,colC);
    
    /* go thru all of pairs in each sewings */
    for (pr = 0; pr < sewing[p]->npair; ++pr)
    {
      Pair_T *const pair = sewing[p]->pair[pr];
      
      populate_F_and_C(patch,pair);/* filling F and C pertinent to this pair */
    }
  }
  
}

/* filling F and C pertinent to this pair */
static void populate_F_and_C(Patch_T *const patch, Pair_T *const pair)
{
  SubFace_T *const subface = pair->subface;
  
  if (!subface->exterF)/* if subface is internal */
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (subface->innerB)/* if there is inner boundary */
  {
    Error0("Wrong subface:\n"
        "It isn't suppoed to have this subface here!\n");
  }
  else if (subface->outerB)/* if it reaches outer boundary */
  {
    Error0("Wrong subface:\n"
        "It isn't suppoed to have this subface here!\n");
  }
  else if (subface->touch)/* if two patches are in touch */
  {
    if (subface->copy)/* if it is collocated point */
    {
      fill_C_F_collocation(patch,pair);
    }
    else
    {
      fill_C_F_interpolation(patch,pair);
    }
  }
  else /* if there is an overlap case */
  {
    fill_C_F_interpolation(patch,pair);
  }
}

/* filling C and F matrices whose entries coming from interpolation points.
// algorithm:
// find variation the following with respect to field(Jacobian):
// ' means derivative with respect to x, y or z
// and N means outward normal on the interface at that point
// y1 = interpolation(y2)
// interpolation(y1) = y2
// Nx.y1_x+Ny.y1_y+Nz.y1_z = Nx.interpolation(y2_x)+Ny.interpolation(y2_y)+Nz.interpolation(y2_z)
// Nx.y2_x+Ny.y2_y+Nz.y2_z = Nx.interpolation(y1_x)+Ny.interpolation(y1_y)+Nz.interpolation(y1_z)
*/
static void fill_C_F_interpolation(Patch_T *const patch, Pair_T *const pair)
{
  /* variable notation convention:
  // let's say we have two patches patch1 and patch2, all of the subscripts 1
  // refer to patch1 and all of the subscripts 2 refer to patch2.
  // some time to improve the readability of the code I defined
  // same variables but with two names. furthermore, we always assume
  // that we are in patch1 and y1 and y2 refer to a generic field
  // that we want to find their F and C; again 1 and 2 refer
  // to their patch numbers. */
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  const SubFace_T *const subface = pair->subface;
  const Uint ppn = pair->patchN;
  const Uint NsubFP1    = subface->np;
  const Uint NsubFP2    = subface->np;
  const Uint NsubM1     = Schur->NS;
  const Uint NinterFP1  = Schur->NI;
  const Uint Fncol = (Uint)Schur->F[ppn]->col;
  double *const F  = Schur->F[ppn]->rmo->A;
  double **const C = Schur->C[ppn]->reg->A;
  double *i2_point;
  double sign;
  Uint subfp1,subfp2,i1,i1_node,i2,s1,s1_node,ip1,ip1_node;
  
  /* "if normal derivatives of the interpolated fields must be continuous" */
  if (subface->df_dn)
  {
    double *N;
    
    /* -(Nx.y1_x+Ny.y1_y+Nz.y1_z)+Nx.interpolation(y2_x)+Ny.interpolation(y2_y)+Nz.interpolation(y2_z) = 0*/
    if (patch->pn == ppn)
    {
      sign = -1;
      const Uint *const inv1 = Schur->inv;
      const Uint *const Imap1 = Schur->Imap;
      const Uint *const Iinv1 = Schur->Iinv;
      const Uint *const node1 = subface->id;
      
      Header_Jacobian
      Init_Jacobian(Jf_D0)
      Init_Jacobian(Jf_D1)
      Init_Jacobian(Jf_D2)

      for (subfp1 = 0; subfp1 < NsubFP1; ++subfp1)
      {
        N = pair->nv[subfp1].N;
        i1_node = node1[subfp1];
        i1 = Imap1[i1_node];
        
        if (i1 == UINT_MAX)
          continue;
        
        Workspace_ijk_Jacobian(i1_node);
        /* F part */
        for (s1 = 0; s1 < NsubM1; ++s1)
        {
          s1_node = inv1[s1];
          Workspace_lmn_Jacobian(s1_node);
          double dfx_df = d2f_dxdu_Jacobian(patch,0,i1_node,s1_node,Jf_D0);
          double dfy_df = d2f_dxdu_Jacobian(patch,1,i1_node,s1_node,Jf_D1);
          double dfz_df = d2f_dxdu_Jacobian(patch,2,i1_node,s1_node,Jf_D2);
          
          F[i_j_to_ij(Fncol,i1,s1)] += sign*(N[0]*dfx_df+N[1]*dfy_df+N[2]*dfz_df);
        }
        /* C part */
        for (ip1 = 0; ip1 < NinterFP1; ++ip1)
        {
          ip1_node = Iinv1[ip1];
          Workspace_lmn_Jacobian(ip1_node);
          
          double dfx_df = d2f_dxdu_Jacobian(patch,0,i1_node,ip1_node,Jf_D0);
          double dfy_df = d2f_dxdu_Jacobian(patch,1,i1_node,ip1_node,Jf_D1);
          double dfz_df = d2f_dxdu_Jacobian(patch,2,i1_node,ip1_node,Jf_D2);
          
          C[i1][ip1] += sign*(N[0]*dfx_df+N[1]*dfy_df+N[2]*dfz_df);
        }
      }
      
      Free_Jacobian(Jf_D0)
      Free_Jacobian(Jf_D1)
      Free_Jacobian(Jf_D2)
      Footer_Jacobian
    }/* end of if (patch->pn == ppn) */
    /* -(Nx.y2_x+Ny.y2_y+Nz.y2_z)+Nx.interpolation(y1_x)+Ny.interpolation(y1_y)+Nz.interpolation(y1_z) = 0 */
    else
    {
      sign = 1;
      const Uint *const Imap2 = pair->sewing->Imap;
      const Uint *const inv1  = Schur->inv;
      const Uint *const Iinv1 = Schur->Iinv;
      const Uint *const node2 = subface->id;
      const Uint plane        = const_index_of_face(patch,subface);
      fdInterp_dfs_T *const dInterp_df_x = get_dInterp_df(patch,subface,"x derivative");
      fdInterp_dfs_T *const dInterp_df_y = get_dInterp_df(patch,subface,"y derivative");
      fdInterp_dfs_T *const dInterp_df_z = get_dInterp_df(patch,subface,"z derivative");
      
      for (subfp2 = 0; subfp2 < NsubFP2; ++subfp2)
      {
        N = pair->nv[subfp2].N;
        i2_point = pair->ip[subfp2].X;
        i2 = Imap2[node2[subfp2]];
        
        if (i2 == UINT_MAX)
          continue;
        
        /* F part */
        for (s1 = 0; s1 < NsubM1; s1++)
        {
          s1_node = inv1[s1];
          F[i_j_to_ij(Fncol,i2,s1)] += sign*(
                       N[0]*dInterp_df_x(patch,i2_point,s1_node,plane)
                       +
                       N[1]*dInterp_df_y(patch,i2_point,s1_node,plane)
                       +
                       N[2]*dInterp_df_z(patch,i2_point,s1_node,plane));
        }
        /* C part */
        for (i1 = 0; i1 < NinterFP1; ++i1)
        {
          i1_node = Iinv1[i1];
          C[i2][i1] += sign*(
                       N[0]*dInterp_df_x(patch,i2_point,i1_node,plane)
                       +
                       N[1]*dInterp_df_y(patch,i2_point,i1_node,plane)
                       +
                       N[2]*dInterp_df_z(patch,i2_point,i1_node,plane));
        }
          
      }
    }/* end of else */
  }/* end of if (subface->df_dn) */
  else/* "if field should be continuous" */
  {
    /* y1-interpolation(y2) = 0 */
    if (patch->pn == ppn)
    {
      sign = 1;
      const Uint *const Imap1 = Schur->Imap;
      const Uint *const node1 = subface->id;
      
      /* F should be zero */
      /* C part */
      for (subfp1 = 0; subfp1 < NsubFP1; ++subfp1)
      {
        i1 = Imap1[node1[subfp1]];
        
        if (i1 == UINT_MAX)
          continue;
          
        C[i1][i1] += 1;
      }
    }
    /* y2-interpolation(y1) = 0 */
    else
    {
      sign = -1;
      const Uint *const Imap2 = pair->sewing->Imap;
      const Uint *const inv1 = Schur->inv;
      const Uint *const Iinv1 = Schur->Iinv;
      const Uint *const node2 = subface->id;
      fdInterp_dfs_T *const dInterp_df = get_dInterp_df(patch,subface,"none");
      
      for (subfp2 = 0; subfp2 < NsubFP2; ++subfp2)
      {
        i2_point = pair->ip[subfp2].X;
        i2 = Imap2[node2[subfp2]];
        
        if (i2 == UINT_MAX)
          continue;
        
        /* F part */
        for (s1 = 0; s1 < NsubM1; s1++)
        {
          s1_node = inv1[s1];
          F[i_j_to_ij(Fncol,i2,s1)] += sign*dInterp_df(patch,i2_point,s1_node,0);
        }
        /* C part */
        for (i1 = 0; i1 < NinterFP1; ++i1)
        {
          i1_node = Iinv1[i1];
          C[i2][i1] += sign*dInterp_df(patch,i2_point,i1_node,0);
        }
      }
    }/* end of else */
  }/* end of else "if field should be continuous" */
}

/* filling C and F matrices whose entries coming from collocation points */
static void fill_C_F_collocation(Patch_T *const patch, Pair_T *const pair)
{
  /* variable notation convention:
  // let's say we have two patches patch1 and patch2, all of the subscripts 1
  // refer to patch1 and all of the subscripts 2 refer to patch2.
  // some time to improve the readability of the code I defined
  // same variables but with two names. furthermore, we always assume
  // that we are in patch 1 and y1 and y2 refer to a generic field
  // that we want to find their F and C; again 1 and 2 refer
  // to their patch numbers. */
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  const SubFace_T *const subface = pair->subface;
  const Uint *const Imap1 = Schur->Imap;
  const Uint *const Imap2 = pair->sewing->Imap;
  const Uint *const inv1 = Schur->inv;
  const Uint *const Iinv1 = Schur->Iinv;
  const Uint ppn = pair->patchN;
  const Uint NsubFP2 = subface->np;
  const Uint NsubM1 = Schur->NS;
  const Uint NinterFP1 = Schur->NI;
  const Uint *node1 = 0,*node2 = 0;
  const Uint Fncol = (Uint)Schur->F[ppn]->col;
  double *const F  = Schur->F[ppn]->rmo->A;
  double **const C = Schur->C[ppn]->reg->A;
  double sign;
  Uint subfp2,i1,i2,s1,s1_node,i1_node,i2_node;
  
  /* if this pair is for the same patch */
  if (patch->pn == ppn)
  {
    node1 = subface->id;
    node2 = subface->id;
  }
  else
  {
    node1 = subface->adjid;
    node2 = subface->id;
  }
    
  if (subface->df_dn)/* n.y1' - n.y2' = 0 => F and C != 0*/
  {
    if (patch->pn == ppn) sign = -1;/* -n.y1'+n.y2'=0 */
    else		  sign =  1;/* -n.y2'+n.y1'=0 */
    
    Header_Jacobian
    Init_Jacobian(Jf_D0)
    Init_Jacobian(Jf_D1)
    Init_Jacobian(Jf_D2)
    
    for (subfp2 = 0; subfp2 < NsubFP2; ++subfp2)
    {
      double *N = pair->nv[subfp2].N;
      i2 = Imap2[node2[subfp2]];
      
      if (i2 == UINT_MAX)
          continue;
      
      i2_node = node1[subfp2];
      Workspace_ijk_Jacobian(i2_node);
      for (s1 = 0; s1 < NsubM1; ++s1)
      {
        s1_node = inv1[s1];
        Workspace_lmn_Jacobian(s1_node);
        double dfx_df = d2f_dxdu_Jacobian(patch,0,i2_node,s1_node,Jf_D0);
        double dfy_df = d2f_dxdu_Jacobian(patch,1,i2_node,s1_node,Jf_D1);
        double dfz_df = d2f_dxdu_Jacobian(patch,2,i2_node,s1_node,Jf_D2);
        
        F[i_j_to_ij(Fncol,i2,s1)] += sign*(N[0]*dfx_df+N[1]*dfy_df+N[2]*dfz_df);
      }
      for (i1 = 0; i1 < NinterFP1; ++i1)
      {
        i1_node = Iinv1[i1];
        Workspace_lmn_Jacobian(i1_node);
        double dfx_df = d2f_dxdu_Jacobian(patch,0,i2_node,i1_node,Jf_D0);
        double dfy_df = d2f_dxdu_Jacobian(patch,1,i2_node,i1_node,Jf_D1);
        double dfz_df = d2f_dxdu_Jacobian(patch,2,i2_node,i1_node,Jf_D2);

        C[i2][i1] += sign*(N[0]*dfx_df+N[1]*dfy_df+N[2]*dfz_df);
      }
    }/* end of for (subfp2 = 0; subfp2 < NsubFP2; ++subfp2) */
    
    Free_Jacobian(Jf_D0)
    Free_Jacobian(Jf_D1)
    Free_Jacobian(Jf_D2)
    Footer_Jacobian
  }/* end of if (subface->df_dn) */
  else/* y1 - y2 = 0 => F is zero */
  {
    if (patch->pn == ppn) sign =  1;/* y1-y2=0 */
    else		  sign = -1;/* y2-y1=0 */
    
    for (subfp2 = 0; subfp2 < NsubFP2; ++subfp2)
    {
      i2 = Imap2[node2[subfp2]];
      
      if (i2 == UINT_MAX)
          continue;
      
      i1 = Imap1[node1[subfp2]];
      /* note: there is a small possibility that i1 becomes UINT_MAX, it means that
      // in finding of the subfaces some points could have been considered as outerboundary
      // but they didn't. it is not really a bug for subfaces, but it causes segfault for this
      // algorithm. if you ever get segfault at this part, be aware of this note! */
      C[i2][i1] += sign*1;
    }
  }
}

/* making E prime and f prime. refer to the note on the very top
// NOTE: if you wanna switch to long version put LONG_VERSION to 1.
// ->return value: some info */
static char *making_E_prime_and_f_prime(Patch_T *const patch)
{
  const double time1 = get_time_sec();
  const int LONG_VERSION = patch->solving_man->settings->umfpack_size;
  DDM_Schur_Complement_T *const S = patch->solving_man->method->SchurC;
  double **E_Trans = 0;
  Matrix_T *a  = 0;
  double *const f = S->f;
  double **xs = 0, **bs = 0;
  Umfpack_T *umfpack = init_umfpack();
  Uint ns = 1;
  Uint i;
  char *msg = calloc(MSG_SIZE2,1);
  IsNull(msg);
  
  /* normal version (= 0), otherwise long */
  if (LONG_VERSION)
    a = cast_matrix_ccs_long(S->B);
  else
    a = cast_matrix_ccs(S->B);
  
  /* free unwanted memories */
  free_matrix(S->B);
  S->B = 0;/* making sure B refers to null */
  
  /* if there is any interface points */
  ns += (Uint)S->E_Trans->row;
    
  xs = calloc(ns,sizeof(*xs));
  IsNull(xs);
  bs = calloc(ns,sizeof(*bs));
  IsNull(bs);
  
  E_Trans = S->E_Trans->reg->A;
  for (i = 0; i < ns-1; ++i)
  {
    bs[i] = E_Trans[i];
    xs[i] = calloc(S->NS,sizeof(*xs[i]));
    IsNull(xs[i]);
  }
  
  bs[ns-1] = f;
  xs[ns-1] = calloc(S->NS,sizeof(*xs[ns-1]));
  IsNull(xs[ns-1]);
  
  umfpack->a = a;
  umfpack->bs = bs;
  umfpack->xs = xs;
  umfpack->ns = ns;
  umfpack->refinement = patch->solving_man->settings->umfpack_refine;
  
  if (LONG_VERSION)
    direct_solver_series_umfpack_dl(umfpack);
  else
    direct_solver_series_umfpack_di(umfpack);
  
  S->f_prime = xs[ns-1];
  /* cast 2d array to row major order */
  
  /* this is an ad-hoc solution, later one can improve it */
  Matrix_T *m_xs = calloc(1,sizeof(*m_xs));
  m_xs->row      = (long)S->E_Trans->row;
  m_xs->col      = (long)S->E_Trans->col;
  m_xs->reg_f    = 1;
  m_xs->reg->A   = xs;

# if E_Trans_prime_IS_RMO_FORMAT
  S->E_Trans_prime = cast_matrix_rmo(m_xs);
  
# elif E_Trans_prime_IS_CCS_FORMAT
  S->E_Trans_prime = cast_matrix_ccs(m_xs);
  
# endif

  free_matrix(m_xs);
  
  sprintf(msg,"{ Solve BE' = E and Bf' = f ...\n"
              "%s"
              "} Solve BE' = E and Bf' = f --> Done. ( Wall-Clock = %.0fs )\n",
              umfpack->description,
              get_time_sec()-time1);
              
  free_matrix(S->E_Trans);
  free_matrix(a);
  free(S->f);
  free(bs);
  free_umfpack(umfpack);
  return msg;
}

/* making B and E matrices. refer to the note on the very top.
// ->return value: some info */
static char *making_B_and_E(Patch_T *const patch)
{
  const double time1 = get_time_sec();
  Solving_Man_T *const S      = patch->solving_man;
  const Uint cf           = S->cf;
  fEquation_T *const jacobian_field_eq = S->jacobian_field_eq[cf];
  fEquation_T *const jacobian_bc_eq = S->jacobian_bc_eq[cf];
  const long Brow = (long)S->method->SchurC->NS;
  const long Bcol = (long)S->method->SchurC->NS;
  const long Erow = (long)S->method->SchurC->NS;
  const long Ecol = (long)S->method->SchurC->NI;
  char *msg = calloc(MSG_SIZE1,1);
  IsNull(msg);
  
  S->method->SchurC->B       = alloc_matrix(REG_SF,Brow,Bcol);
  S->method->SchurC->E_Trans = alloc_matrix(REG_SF,Ecol,Erow);
    
  jacobian_field_eq(patch,S->method->SchurC);
  jacobian_bc_eq(patch,S->method->SchurC);
  
  sprintf(msg,"{ Populate B and E ...\n"
              "} Populate B and E --> Done. ( Wall-Clock = %.0fs )\n",
                get_time_sec()-time1);
  return msg;
}

/* making B when the grid only has one patch.
// refer to the note on the very top. */
static void making_B_single_patch(Patch_T *const patch)
{
  Solving_Man_T *const S      = patch->solving_man;
  const Uint cf           = S->cf;
  fEquation_T *const jacobian_field_eq = S->jacobian_field_eq[cf];
  fEquation_T *const jacobian_bc_eq = S->jacobian_bc_eq[cf];
  const long Brow = (long)S->method->SchurC->NS;
  const long Bcol = (long)S->method->SchurC->NS;
  
  S->method->SchurC->B       = alloc_matrix(REG_SF,Brow,Bcol);
    
  jacobian_field_eq(patch,S->method->SchurC);
  jacobian_bc_eq(patch,S->method->SchurC);
  
}

/* having been made partial g's, now make g column */
static void make_g(Grid_T *const grid)
{
  Patch_T *patch;
  DDM_Schur_Complement_T *Schur;
  Sewing_T *sewing;
  Pair_T *pair1, *pair2;
  Uint *Imap;/* interface map */
  Uint *Smap;/* subface map */
  double *g,*pg1,*pg2;/* g and partial g's */
  Uint npair,NSubFP;
  Uint p,s,pr;
  
  FOR_ALL_PATCHES(p,grid)
  {
    patch  = grid->patch[p];
    Schur  = patch->solving_man->method->SchurC;
    sewing = Schur->sewing[p];
    
    if (!sewing)
      continue;
    
    npair  = sewing->npair;
    Imap   = Schur->Imap;
      
    Schur->g = alloc_double(Schur->NI);
    g = Schur->g;
      
    for (pr = 0; pr < npair; ++pr)
    {
      pair1  = sewing->pair[pr];
      pair2  = pair1->mirror;
      NSubFP = pair1->subface->np;
      Smap   = pair1->subface->id;
      pg1    = pair1->pg;
      pg2    = pair2->pg;
      
      for (s = 0; s < NSubFP; ++s)
      {
        Uint s_node = Smap[s];
        if (Imap[s_node] == UINT_MAX)
          continue;
        /* NOTE: because we've filled each pg's in order 
        // now we can add them with same indices; otherwise we couldn't. */
        g[Imap[s_node]] += pg1[s]+pg2[s];
      }
      
      free(pg1);
      free(pg2);
    }
    
  }
}

/* since Shur Complement needs:
// 1. specific labeling for grid 
// 2. boundary points counted only once
// 3. adjacency patch of each patch and how other 
// patches use this patch to impose their boundary conditions.
// this function provides these ingredients, among others.
*/
static void preparing_ingredients(Solve_Equations_T *const SolveEqs)
{
  Grid_T *const grid = get_grid_solve_equations(SolveEqs);
  
  Uint p;
  Uint count = 0;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    /* if this patch doesn't have Schur stucture */
    if (!patch->solving_man->method->Schur_Complement)
      count ++;
  }
  
  if (!count)/* if all structures are ready, exit */
    return;
  else if (count != 0 && count < grid->np)/* is some structures have been made and some not */
    Error0(" How come that some Schur structure are ready and some are not!\n");
  
  /* populating Schur sturct */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
      
    DDM_Schur_Complement_T *SchurC = calloc(1,sizeof(*SchurC));
    IsNull(SchurC);
    SchurC->map = malloc(patch->nn*sizeof(*SchurC->map));
    IsNull(SchurC->map);
    SchurC->inv = malloc(patch->nn*sizeof(*SchurC->inv));
    IsNull(SchurC->inv);  
    SchurC->patch = patch;
    SchurC->np    = patch->grid->np;
    patch->solving_man->method->Schur_Complement = 1;
    patch->solving_man->method->SchurC =  SchurC;
    
    /* making map and inv */
    make_map_and_inv(patch);
    /* populating sewing structure for boundary info */
    populate_sewing(patch);
  }
  
  /* adding more info and struct */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    mirror_pairs(patch);
    miscellany_in_sewing(patch);
    set_NSs_NIs(patch);
  }
  checks_and_constraints(grid);
  
}

/* conduct some tests and checkups */
static void checks_and_constraints(const Grid_T *const grid)
{
  DDM_Schur_Complement_T *Schur;
  Uint p,count;
  
  count = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (patch->outerB && patch->innerB)
      Error0("In this version of domain decompostion method it is assumed \n"
              "a patch can only have one outer boundary or inner bounday at a time;\n"
              "but, it seems this patch has both outer boundary and inner bounday at a time!\n");
    
    Schur = patch->solving_man->method->SchurC;
    if (!Schur->NI)
      count++;
    
  }
  
  if (count > 1)
    Error0("Disconnected Manifold: It seems that the grid has some gaps in it!\n"
    " At least two of the patches are disconnected.");
}

/* set NS_p, NI_p, NS_total and NI_total for each patch */
static void set_NSs_NIs(Patch_T *const patch)
{	
  Grid_T *const grid = patch->grid;
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  Uint np = grid->np;
  Uint NS_total = 0;
  Uint NI_total = 0;
  Uint *NS_p = calloc(np,sizeof(*NS_p));
  Uint *NI_p = calloc(np,sizeof(*NI_p));
  IsNull(NS_p);
  IsNull(NI_p);
  Uint p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch_p = grid->patch[p];
    
    NS_p[p]   = patch_p->solving_man->method->SchurC->NS;
    NS_total += NS_p[p];
    NI_p[p]   = patch_p->solving_man->method->SchurC->NI;
    NI_total += NI_p[p];
  }
  Schur->NS_p     = NS_p;
  Schur->NS_total = NS_total;
  Schur->NI_p     = NI_p;
  Schur->NI_total = NI_total;
}

/* since some of the info should be added after all other sewings have
// already been made, this function is needed. */
static void miscellany_in_sewing(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const S = patch->solving_man->method->SchurC;
  DDM_Schur_Complement_T *S2;
  Sewing_T *sewing;
  Grid_T *const grid = patch->grid;
  const Uint np  = grid->np;
  const Uint cp  = patch->pn;
  Uint p;
  
  /* go thru all sewings */
  for (p = 0; p < np; ++p)
  {
    sewing = S->sewing[p];
    
    if (!sewing)
      continue;
      
    else if (p == cp)
    {
      sewing->NS   = S->NS;
      sewing->NI   = S->NI;
      sewing->Oi   = S->Oi;
      /* note the following are soft copy */
      sewing->map  = S->map;
      sewing->inv  = S->inv;
      sewing->Imap = S->Imap;
      sewing->Iinv = S->Iinv;
    }
    else
    {
      S2 = grid->patch[p]->solving_man->method->SchurC;
      sewing->NS   = S2->NS;
      sewing->NI   = S2->NI;
      sewing->Oi   = S2->Oi;
      /* note the following are hard copy */
      sewing->map  = dup_UINT(S2->map,S2->NS+S2->NI);
      sewing->inv  = dup_UINT(S2->inv,S2->NS+S2->NI);
      sewing->Imap = dup_UINT(S2->Imap,S2->NS+S2->NI);
      sewing->Iinv = dup_UINT(S2->Iinv,S2->NI);
    }
  }
}

/* populating sewing struct in each patch using subfaces
// and geometry of grid. in effect, we duplicate subfaces in order
// for each patch can make boundary condition without causing racing
// condition for each thread. thus, each patch knows exactly what others
// and itself need for boundary conditions like interpolation and continuity.
*/
static void populate_sewing(Patch_T *const patch)
{
  Grid_T *const grid = patch->grid;
  DDM_Schur_Complement_T *const SchurC = 
                                patch->solving_man->method->SchurC;
  const Uint np = patch->grid->np;
  Sewing_T **sewing = 0;
  Uint p;
  
  sewing = calloc(np,sizeof(*sewing));
  IsNull(sewing);
  
  /* initialize sewings and pairs */
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch2 = grid->patch[p];
    
    if (patch2 == patch)
    {
      make_its_sewing(patch,sewing);
    }
    else
    {
      make_others_sewing(patch,patch2,sewing);
    }
  }
  
  SchurC->sewing  = sewing;
  SchurC->nsewing = np;
}

/* mirror pairs of different patches */
static void mirror_pairs(Patch_T *const patch)
{
  Grid_T *const grid = patch->grid;
  DDM_Schur_Complement_T *const SchurC = 
                                patch->solving_man->method->SchurC;
  Sewing_T *const sewing = SchurC->sewing[patch->pn];
  Uint p;
  
  /* if there is no sewing, happen for single patch */
  if (!sewing)
    return;
  
  for (p = 0; p < sewing->npair; ++p)
  {
    Pair_T *pair = sewing->pair[p];
    Uint ap = pair->subface->adjPatch;
    Sewing_T *sewing2 = 
          grid->patch[ap]->solving_man->method->SchurC->sewing[patch->pn];
    pair->mirror = find_pair_in_sewing(sewing2,pair->subface);
    assert(pair->mirror);
  }
}  

/* finding the pair for a given subface and sewing. it finds based on
// face number and subface number of the given subface
// ->return value: found pair, 0 otherwise.
*/
static Pair_T *find_pair_in_sewing(const Sewing_T *const sewing,const SubFace_T *const subface)
{
  Pair_T *pair = 0;
  Uint i;
  
  for (i = 0; i < sewing->npair; ++i)
  {
    SubFace_T *s = sewing->pair[i]->subface;
    
    if (s->patch != subface->patch)
      Error0("These subfaces are supposed to be on a same patches!\n");
      
    if (s->face == subface->face && s->sn == subface->sn)
    {
      pair = sewing->pair[i];
      break;
    }
  }
  
  return pair;
}

/* making all of sewings deduced form its own subfaces */
static void make_its_sewing(const Patch_T *const patch,Sewing_T **const sewing)
{
  const Uint nintfc = countf(patch->interface);
  const Uint p = patch->pn;
  Uint intfc;

  /* loop over all interfaces */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    Uint nsfc = interface->ns;
    Uint sfc;
    
    /* loop over all subfaces */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      
      if (!subface->exterF)
      {
        Error0(INCOMPLETE_FUNC);
      }
      else if (subface->innerB)
      {
        continue;
      }
      else if (subface->outerB)
      {
        continue;
      }
      else
      {
        if (!sewing[p])
        {
          sewing[p] = alloc_sewing();
          sewing[p]->patchN = p;
        }
          
        populate_pair(sewing[p],subface,ITS);
      }
      
    }/* end of for (sfc = 0; sfc < nsfc; ++sfc) */
  }/* end of for (intfc = 0; intfc < nintfc; ++intfc) */
  
}

/* making all of sewings deduced from subfaces of patch2 */
static void make_others_sewing(const Patch_T *const patch,const Patch_T *const patch2,Sewing_T **const sewing)
{
  const Uint nintfc = countf(patch2->interface);
  const Uint p2 = patch2->pn;
  Uint intfc;

  /* loop over all interfaces */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch2->interface[intfc];
    Uint nsfc = interface->ns;
    Uint sfc;
    
    /* loop over all subfaces */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      
      if (!subface->exterF)
      {
        Error0(INCOMPLETE_FUNC);
      }
      else if (subface->innerB)
      {
        continue;
      }
      else if (subface->outerB)
      {
        continue;
      }
      else if (subface->adjPatch == patch->pn)
      {
        if (!sewing[p2])
        {
          sewing[p2] = alloc_sewing();
          sewing[p2]->patchN = p2;
        }
          
        populate_pair(sewing[p2],subface,OTHERS);
      }
   
    }/* end of for (sfc = 0; sfc < nsfc; ++sfc) */
  }/* end of for (intfc = 0; intfc < nintfc; ++intfc) */
} 

/* populate pair in sewing.
// if flag == ITS, it only connects subfaces pointers, 
// if flag == OTHERS, it duplicate the subface.
*/
static void populate_pair(Sewing_T *const sewing,SubFace_T *const subface,const DDM_SC_Flag_T flag)
{
  Grid_T *const grid = subface->patch->grid;
  Pair_T *const pair = calloc(1,sizeof(*pair));
  const Uint np = subface->np;
  Uint i;
  IsNull(pair);
  
  if (flag == OTHERS)
  {
    pair->subface = calloc(1,sizeof(*pair->subface));
    IsNull(pair->subface);
    copy_subface(pair->subface,subface);
  }
  else if (flag == ITS)
    pair->subface = subface;
  else
    Error0("Wrong flag.\n");
  
  
  /* if this subface needs normal vector */
  if (subface->df_dn)
  {
    double *N;
    Point_T point;
    
    pair->nv = calloc(subface->np,sizeof(*pair->nv));
    IsNull(pair->nv);
    
    point.patch = subface->patch;
    point.face  = subface->face;
    
    for (i = 0; i < np; ++i)
    {
      point.ind = subface->id[i];
      N = normal_vec(&point);
      pair->nv[i].N[0] = N[0];
      pair->nv[i].N[1] = N[1];
      pair->nv[i].N[2] = N[2];  
    }
    
  }
  
  /* if this interface is an interpolation one and is for others */
  if (!subface->copy && flag == OTHERS)
  {
    const Uint *const node = subface->id;
    const Patch_T *adjPatch = grid->patch[subface->adjPatch];
    double X[3];
    
    pair->ip = calloc(subface->np,sizeof(*pair->ip));
    IsNull(pair->ip);
    
    for (i = 0; i < np; ++i)
    {
      double *x = subface->patch->node[node[i]]->x;
      assert(X_of_x_precision(X,x,adjPatch,subface->precision_factor));
      
      pair->ip[i].X[0] = X[0];
      pair->ip[i].X[1] = X[1];
      pair->ip[i].X[2] = X[2];
    }
  }
  
  pair->sewing = sewing;
  pair->patchN = sewing->patchN;
  sewing->pair = 
    realloc(sewing->pair,(sewing->npair+1)*sizeof(*sewing->pair));
  IsNull(sewing->pair);
  sewing->pair[sewing->npair] = pair;
  sewing->npair++;
}

/* making map and inv map for Schur complement method.
// inner mesh nodes come first, then outer boundary nodes 
// and then finally inner boundar nodes.
*/
static void make_map_and_inv(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const SchurC = 
    patch->solving_man->method->SchurC;
  Uint *const map = SchurC->map;
  Uint *const inv = SchurC->inv;
  Uint *Imap = 0;
  Uint *Iinv = 0;
  const Uint nn = patch->nn;
  const Uint *n = patch->n;
  const Uint nintfc = countf(patch->interface);
  Uint i,j,j2, intfc;/* dummy indices and counters */
  
  /* keep tracking of points, 1 means counted, 0 means not */
  Uint *flag_point = calloc(nn,sizeof(*flag_point));
  IsNull(flag_point);
  
  if (patch->is_a_closed || patch->is_b_closed || patch->is_c_closed)
    Error0(INCOMPLETE_FUNC);
    
  /* filling inner points */
  j = 0;
  for (i = 0; i < nn; ++i)
  {
    if (!OnFace(n,i))
    {
      map[i] = j;
      inv[j] = i;
      flag_point[i] = 1;
      j++;
    }
  }
  
  /* set initial index for outer boundary points. */
  patch->solving_man->method->SchurC->Oi = j;
  
  /* filling boundary points */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    Uint nsfc = interface->ns;
    Uint sfc;
    
    /* loop over all subfaces to filling boundary */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      const Uint *const B = subface->id;
      
      if (subface->outerB || subface->innerB)/* if it reaches any kind of boundaries */
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
  
  /* all of subdomain points already considered so : */
  patch->solving_man->method->SchurC->NS = j;
  patch->solving_man->method->SchurC->NI = nn-j;
  
  j2 = 0;
  Imap = calloc(nn,sizeof(*Imap));
  IsNull(Imap);
  Iinv = calloc(nn-j,sizeof(*Iinv));
  IsNull(Iinv);
  
  /* make sure if it is given a point outside of its domain
  // it returns UINT_MAX. */
  for (i = 0; i < nn; ++i)
    Imap[i] = UINT_MAX;
  
  /* filling the other remaining points */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    Uint nsfc = interface->ns;
    Uint sfc;
    
    /* loop over all subfaces to filling boundary points */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      Uint *T = subface->id;
      
      if (!subface->exterF)/* if subface is internal */
      {
        Error0(INCOMPLETE_FUNC);
      }
      else if (subface->innerB)/* if there is inner boundary */
      {
        continue;
      }
      else if (subface->outerB)/* if it reaches outer boundary */
      {
        continue;
      }
      else
      {
        for (i = 0; i < subface->np; ++i)
          if (!flag_point[T[i]])
          {
            map[T[i]] = j;
            inv[j] = T[i];
            Imap[T[i]] = j2;
            Iinv[j2] = T[i];
            flag_point[T[i]] = 1;
            ++j;
            ++j2;
          }
      }
    }
  }
  
  if (j != nn)
  {
    Error0("Not all points are mapped.\n");
  }
  
  SchurC->Imap = Imap;
  SchurC->Iinv = Iinv;
  
  free(flag_point);
}

/* read the names of the fields to be solved
// in order given in the input file.
// the parameter is called solve_Order
// ->return value: name and number of fields. */
char **get_solving_field_name(const char *const solving_order,Uint *const nf)
{
  char COMMA = ',';
  char *par;
  char **field_name = 0;
  char *tok,*save = 0;
  *nf = 0;
  
  /* finding fields's name */
  par = dup_s(solving_order);/* par = f1,f2,... */
  tok = tok_s(par,COMMA,&save);/* tok = f1 */
  while(tok)
  {
    field_name = realloc(field_name,(*nf+1)*sizeof(*field_name));
    IsNull(field_name);
    field_name[*nf] = dup_s(tok);
    tok = tok_s(0,COMMA,&save);/* tok = f2 */
    (*nf)++;
  }
  free(par);
  
  return field_name;
}

/* set solving_man->cf, cf stands for current field */
static void set_solving_man_cf(Solve_Equations_T *const SolveEqs)
{
  Grid_T *const grid = get_grid_solve_equations(SolveEqs);
  const char *const field_name = SolveEqs->field_name;
  Uint p;
  
  assert(field_name);
    
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    char **fields = patch->solving_man->field_name;
    Uint nf = patch->solving_man->nf;
    patch->solving_man->cf = find_index_string(fields,nf,field_name);
  }
}

/* set solving_man->settings */
static void set_solving_man_settings(Solve_Equations_T *const SolveEqs)
{
  Grid_T *const grid = get_grid_solve_equations(SolveEqs);
  Uint p;
    
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    /* relaxation factor: */
    patch->solving_man->settings->relaxation_factor = 
      get_relaxation_factor_solve_equations(SolveEqs);
    
    patch->solving_man->settings->Frms_i  = DBL_MAX;
    patch->solving_man->settings->HFrms   = 0;
    patch->solving_man->settings->NHFrms  = 0;
    patch->solving_man->settings->solver_step  = 0;
    patch->solving_man->settings->last_sol = 0;
    patch->solving_man->settings->umfpack_size = SolveEqs->umfpack_size;
    patch->solving_man->settings->umfpack_refine = SolveEqs->umfpack_refine;
    
  }
}

/* making f column.
// THEARD SAFA.
*/
static void make_f(Patch_T *const patch)
{
  DDM_Schur_Complement_T *Schur = patch->solving_man->method->SchurC;
  
  Schur->f = alloc_double(Schur->NS);
  
  f_in_equation_part(patch);
  f_in_boundary_part(patch);
  
}

/* making partial of g which coming from this patch.
// for example if g = y1-y2, g partly coming from patch1
// and partly coming from patch2. when all of partial g's are made
// the g coloumn itself will be populated.
*/
static void make_partial_g(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const Schur = patch->solving_man->method->SchurC;
  Sewing_T **const sewing = Schur->sewing;
  const Uint np = Schur->np;
  Uint p;
  
  /* go thru all of sewings  */
  for (p = 0; p < np; ++p)
  {
    Uint pr;
    
    if (!sewing[p])
      continue;
      
    /* go thru all of pairs in each sewings */
    for (pr = 0; pr < sewing[p]->npair; ++pr)
    {
      Pair_T *const pair = sewing[p]->pair[pr];
      make_pg(patch,pair);/* make partial g pertinent to this pair */
    }
  }
}

/* make partial g pertinent to this pair */
static void make_pg(Patch_T *const patch, Pair_T *const pair)
{
  SubFace_T *const subface = pair->subface;
  
  pair->pg = alloc_double(subface->np);
  
  if (!subface->exterF)/* if subface is internal */
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (subface->innerB)/* if there is inner boundary */
  {
    Error0("Wrong subface:\n"
        "It isn't suppoed to have this subface here!\n");
  }
  else if (subface->outerB)/* if it reaches outer boundary */
  {
    Error0("Wrong subface:\n"
        "It isn't suppoed to have this subface here!\n");
  }
  else if (subface->touch)/* if two patches are in touch */
  {
    if (subface->copy)/* if it is collocated point */
    {
      pg_collocation(patch,pair);
    }
    else
    {
      pg_interpolation(patch,pair);
    }
  }
  else /* if there is an overlap case */
  {
    pg_interpolation(patch,pair);
  }
}

/* making pg whose part coming from collocation points,
// thus, we need copy the field or their normal derivatives.
*/
static void pg_collocation(Patch_T *const patch, Pair_T *const pair)
{
  SubFace_T *const subface = pair->subface;
  const Uint cf = patch->solving_man->cf;
  const char *const field_aliased = patch->solving_man->field_aliased[cf];
  Field_T *const f = patch->fields[Ind(field_aliased)];
  const Uint ppn = pair->patchN;
  const Uint NSubFP = subface->np;
  const Uint *node = 0;
  double *const pg = pair->pg; 
  double sign;
  Uint s,s_node;
  
  /* if this pg is for the same patch */
  if (patch->pn == ppn)
  {
    node = subface->id;
  }
  else
  {
    node = subface->adjid;
  }
    
  if (subface->df_dn)
  {
    if (patch->pn == ppn) sign = -1;/* -n.y1'+n.y2'=0 */
    else                  sign =  1;/* -n.y2'+n.y1'=0 */

    double *f_x, *f_y, *f_z;
    
    f_x = Partial_Derivative(f,"x");
    f_y = Partial_Derivative(f,"y");
    f_z = Partial_Derivative(f,"z");
    
    for (s = 0; s < NSubFP; ++s)
    {
      double *N = pair->nv[s].N;
      s_node = node[s];
      pg[s] += sign*(N[0]*f_x[s_node]+N[1]*f_y[s_node]+N[2]*f_z[s_node]);
    }
    
    free(f_x);
    free(f_y);
    free(f_z);
  }/* end of if (subface->df_dn) */
  else
  {
    if (patch->pn == ppn) sign =  1;/* y1-y2=0 */
    else		  sign = -1;/* y2-y1=0 */
    
    for (s = 0; s < NSubFP; ++s)
    {
      s_node = node[s];
      pg[s] += sign*f->v[s_node];
    }
  }
}

/* making pg whose part coming from interpolation points,
// thus, we need copy the field or their normal derivatives
// at interpolation points.
*/
static void pg_interpolation(Patch_T *const patch, Pair_T *const pair)
{
  SubFace_T *const subface = pair->subface;
  const Uint cf = patch->solving_man->cf;
  const char *const field_aliased = patch->solving_man->field_aliased[cf];
  Field_T *const f = patch->fields[Ind(field_aliased)];
  const Uint np = subface->np;
  const Uint ppn = pair->patchN;
  double *const pg = pair->pg; 
  double sign;
  Uint i;
  
  /* if normal derivatives of fields should be continuous */
  if (subface->df_dn)
  {
    Patch_T tmp_patch = make_temp_patch(patch);
    Field_T *f_x = add_field("df_dx","(3dim)",&tmp_patch,NO),
            *f_y = add_field("df_dy","(3dim)",&tmp_patch,NO), 
            *f_z = add_field("df_dz","(3dim)",&tmp_patch,NO);
    double *N;
    
    f_x->v = Partial_Derivative(f,"x");
    f_y->v = Partial_Derivative(f,"y");
    f_z->v = Partial_Derivative(f,"z");
    
    /* -(Nx.y1_x+Ny.y1_y+Nz.y1_z)+Nx.interpolation(y2_x)+Ny.interpolation(y2_y)+Nz.interpolation(y2_z) = 0*/
    if (patch->pn == ppn)
    {
      sign = -1;
      Uint *node = subface->id;
      
      for (i = 0; i < np; ++i)
      {
        N = pair->nv[i].N;
        pg[i] += sign*(
                  N[0]*f_x->v[node[i]]
                  +
                  N[1]*f_y->v[node[i]]
                  +
                  N[2]*f_z->v[node[i]]
                  );
      }
    }
    /* -(Nx.y2_x+Ny.y2_y+Nz.y2_z)+Nx.interpolation(y1_x)+Ny.interpolation(y1_y)+Nz.interpolation(y1_z) = 0 */
    else
    {
      sign = 1;
      double *X;
      Interpolation_T *interp_x = init_interpolation();
      Interpolation_T *interp_y = init_interpolation();
      Interpolation_T *interp_z = init_interpolation();
      interp_x->field = f_x;
      interp_y->field = f_y;
      interp_z->field = f_z;
      
      fill_interpolation_flags(interp_x,patch,subface);
      fill_interpolation_flags(interp_y,patch,subface);
      fill_interpolation_flags(interp_z,patch,subface);
      plan_interpolation(interp_x);
      plan_interpolation(interp_y);
      plan_interpolation(interp_z);
      
      for (i = 0; i < np; ++i)
      {
        X = pair->ip[i].X;
        N = pair->nv[i].N;
        
        interp_x->X = X[0];
        interp_x->Y = X[1];
        interp_x->Z = X[2];
        
        interp_y->X = X[0];
        interp_y->Y = X[1];
        interp_y->Z = X[2];
        
        interp_z->X = X[0];
        interp_z->Y = X[1];
        interp_z->Z = X[2];

        pg[i] += sign*(
                  N[0]*execute_interpolation(interp_x)
                  +
                  N[1]*execute_interpolation(interp_y)
                  +
                  N[2]*execute_interpolation(interp_z)
                  );
      }
      free_interpolation(interp_x);
      free_interpolation(interp_y);
      free_interpolation(interp_z);
      
    }
    
    remove_field(f_x);
    remove_field(f_y);
    remove_field(f_z);
    free_temp_patch(&tmp_patch);
  }/* end of if (subface->df_dn) */
  /* if field should be continuous */
  else
  {
    /* y1-interpolation(y2) = 0 */
    if (patch->pn == ppn)
    {
      sign = 1;
      Uint *node = subface->id;
      
      for (i = 0; i < np; ++i)
      {
        pg[i] += sign*f->v[node[i]];
      }
    }
    /* y2-interpolation(y1) = 0 */
    else
    {
      sign = -1;
      Interpolation_T *interp = init_interpolation();
      interp->field = f;
      fill_interpolation_flags(interp,patch,subface);
      plan_interpolation(interp);

      for (i = 0; i < np; ++i)
      {
        interp->X = pair->ip[i].X[0];
        interp->Y = pair->ip[i].X[1];
        interp->Z = pair->ip[i].X[2];
        
        pg[i] += sign*execute_interpolation(interp);
      }
      
      free_interpolation(interp);
    }
  }
}

/* filling flags of iterpolation based on the given subface. */
static void fill_interpolation_flags(Interpolation_T *const it,Patch_T *const patch,const SubFace_T *const sf)
{
  if (sf->sameX)
  {
    it->YZ_dir_flag = 1;
    it->I = const_index_of_face(patch,sf);
  }
  else if (sf->sameY)
  {
    it->XZ_dir_flag = 1;
    it->J = const_index_of_face(patch,sf);
  }
  else if (sf->sameZ)
  {
    it->XY_dir_flag = 1;
    it->K = const_index_of_face(patch,sf);
  }
  else if (!sf->sameX && !sf->sameY && !sf->sameZ)
  {
    it->XYZ_dir_flag = 1;
  }
  else if (sf->sameX && sf->sameY && sf->sameZ)
  {
    Error0("How come to have all the sameX,Y,Z flags of subface be the same.\n"
    "It means there is something wrong at finding of adjacent patches.\n");
  }

}

/* it finds the equation of a face(plane) in index format.
// note: we need this function for case like interpolation and etc.
// thus the sf subface is the subface of juxtapose patch which toches the patch,
// so sf->adjFace is a face in patch.
// more explanation:
// given subface, if it is a juxtapose kind (touch), 
// it returns the coordinate index of the surface of neighbor patch 
// in which this interface touches it.
// e.g. if two patches are juxtapose at X[0] = 4.5 in patch A and in this patch
// it is correspond to the index I in patch->node[i_j_k_to_ijk(n,I,*,*,)]->X[0] so this function
// returns I in case the subface of patch B which touches the mentioned interface is give.
// ->return value: coordinate index of plane X = const. if they won't touch it gives UINT_MAX. */
// ->return value: constant index(coords) of a given face, if not found UINT_MAX. */
static Uint const_index_of_face(Patch_T *const patch,const SubFace_T *const sf)
{
  const Uint f = sf->adjFace;
  const Uint *const n = patch->n;
  Uint C = UINT_MAX;/* constant value */
  
  if (sf->touch)
  {
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
        Error0(NO_OPTION);
    }
  }
  
  return C;
}


/* calculating the part of f coming from equation */
static void f_in_equation_part(Patch_T *const patch)
{
  Solving_Man_T *const S      = patch->solving_man;
  fEquation_T *const field_eq = S->field_eq[S->cf];
  
  field_eq(patch,S->method->SchurC);
}

/* calculating the part of f coming from boundary points */
static void f_in_boundary_part(Patch_T *const patch)
{
  const Uint nintfc        = countf(patch->interface);
  Solving_Man_T *const S   = patch->solving_man;
  const Uint cf            = S->cf;
  fEquation_T *const bc_eq = S->bc_eq[cf];
  Boundary_Condition_T bc;
  Uint intfc;
  bc.patch   = patch;
  
  /* loop over all interfaces */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    Uint nsfc = interface->ns;
    Uint sfc;
    
    /* loop over all subfaces and look for boundary */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      
      if (subface->outerB || subface->innerB)/* if it reaches any kind of boundaries */
      {
        bc.subface = subface;
        bc.node    = subface->id;
        bc.nn      = subface->np;
        bc_eq(&bc,S->method->SchurC);
      }
    }
  }
}



/* calculate root mean square of F, in Jx=-F, for the whole grid when
// it has only one patch and set it in patch->solving_man->Frms */
static void calculate_residual_single_patch(Patch_T *const patch)
{
  DDM_Schur_Complement_T *S = patch->solving_man->method->SchurC;
  double *f = S->f;
  double sqr = dot(S->NS,f,f);
  
  printf("\nResidual:\n");
  patch->solving_man->Frms = sqrt(sqr/S->NS);
  printf("-------->%s = %e\n", patch->name,patch->solving_man->Frms);
}

/* calculate root mean square of F, in Jx=-F, for the whole grid. 
// and set it in patch->solving_man->settings->Frms_i for single patch */
static void set_solving_man_settings_Frms_i_single_patch(Patch_T *const patch)
{
  DDM_Schur_Complement_T *S = patch->solving_man->method->SchurC;
  double *f = S->f;
  double sqr = dot(S->NS,f,f);
  double *HFrms = alloc_double(1);
  
  patch->solving_man->settings->Frms_i = sqrt(sqr/S->NS);
  HFrms[0] = patch->solving_man->settings->Frms_i;
  patch->solving_man->settings->HFrms  = HFrms;
  patch->solving_man->settings->NHFrms = 1;
}

/* calculate root mean square of F, in Jx=-F, for the whole grid. 
// and set it in patch->solving_man->settings->Frms_i */
static void set_solving_man_settings_Frms_i(Grid_T *const grid)
{
  const Uint npatch = grid->np;
  Uint p;
  
  DDM_SCHUR_OpenMP(omp parallel for)
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    DDM_Schur_Complement_T *S = patch->solving_man->method->SchurC;
    double *f = S->f;
    double *g = S->g;
    double sqr1 = dot(S->NS,f,f);
    double sqr2 = dot(S->NI,g,g);
    double sqrs = sqr1+sqr2;
    patch->solving_man->settings->Frms_i = sqrt(sqrs/(S->NS+S->NI));
  }
}

// set the current step of solver in patch->solving_man->settings->solver_step 
// it is used for stop criteria */
static void set_solving_man_settings_solver_step(Grid_T *const grid,const int current_step)
{
  const Uint npatch = grid->np;
  Uint p;
  
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    patch->solving_man->settings->solver_step = current_step;
  }
}

/* calculate root mean square of F, in Jx=-F, for the whole grid. 
// and set it in patch->solving_man->Frms */
static void calculate_residual(Grid_T *const grid)
{
  const Uint npatch = grid->np;
  double *HFrms = alloc_double(npatch);
  Uint NHFrms;
  double *extd = 0;
  Uint p,i;
  
  DDM_SCHUR_OpenMP(omp parallel for)
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    DDM_Schur_Complement_T *S = patch->solving_man->method->SchurC;
    double *f = S->f;
    double *g = S->g;
    double sqr1 = dot(S->NS,f,f);
    double sqr2 = dot(S->NI,g,g);
    double sqrs = sqr1+sqr2;
    patch->solving_man->Frms = sqrt(sqrs/(S->NS+S->NI));
    HFrms[p]                 = patch->solving_man->Frms;
  }
  
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    NHFrms = patch->solving_man->settings->NHFrms;
    
    extd = patch->solving_man->settings->HFrms;
    extd = realloc(extd,(NHFrms+1)*sizeof(*extd));
    IsNull(extd);
    extd[NHFrms] = HFrms[p];
    patch->solving_man->settings->HFrms = extd;
    patch->solving_man->settings->NHFrms++;
  }
  free(HFrms);
  HFrms = 0;
  
  /* print residual */
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Uint cf = patch->solving_man->cf;
    const char *field_name = patch->solving_man->field_name[cf];
    
    printf("\nResidual History:\n");
    printf("|---> %s equation:\n",field_name);
    printf("      |---> %s:\n", patch->name);
    HFrms  = patch->solving_man->settings->HFrms;
    NHFrms = patch->solving_man->settings->NHFrms;
    for (i = 0; i < NHFrms; ++i)
      printf("%*s|---> Newton step %d: %e\n",12," ",i,HFrms[i]);
    printf("\n");
    fflush(stdout);
  }
}

/* free patch->solving_man->settings */
static void free_solving_man_settings(Grid_T *const grid)
{
  const Uint npatch = grid->np;
  Uint p;
  
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    Free(patch->solving_man->settings->HFrms);
    patch->solving_man->settings->HFrms  = 0;
    patch->solving_man->settings->NHFrms = 0;
    
    Free(patch->solving_man->settings->last_sol);
    patch->solving_man->settings->last_sol = 0;
  }
   
}

/* figure out if a given node is located on a face or not,
// base on number of points in each direction.
// ->return value: if is on face 1, otherwise 0.
*/
static Uint OnFace(const Uint *const n, const Uint p)
{
  Uint i,j,k;
  
  ijk_to_i_j_k(p,n,&i,&j,&k);
  
  if (i == n[0]-1 || i == 0)  return 1;
  if (j == n[1]-1 || j == 0)  return 1;
  if (k == n[2]-1 || k == 0)  return 1;
  
  return 0;
}

/* solving Bx = f when grid only has only one patch. */
static void solve_Bx_f(Patch_T *const patch)
{
  DDM_Schur_Complement_T *Schur = patch->solving_man->method->SchurC;
  const Uint NS = Schur->NS;
  double *f = Schur->f;
  double *x = alloc_double(NS);
  Matrix_T *B = Schur->B;
  Matrix_T *B_ccs = cast_matrix_ccs(B);
  Umfpack_T *umfpack = init_umfpack();
  
  free_matrix(B);
  
  umfpack->a = B_ccs;
  umfpack->b = f;
  umfpack->x = x;
  umfpack->refinement = patch->solving_man->settings->umfpack_refine;
  direct_solver_umfpack_di(umfpack);
  
  free_matrix(B_ccs);
  free(f);
  free_umfpack(umfpack);
  Schur->x = x;
}

/* verifying the jacobian of equations, if they have been written correctly.
// by jacobian we mean matrix J in Jx = -F in Newton method. 
// it also could be used for testing of jacobian made by Shur method. */
void test_Jacobian_of_equations(Solve_Equations_T *const SolveEqs)
{
  char **field_name = 0;/* name of all fields to be solved */
  Uint nf = 0;/* number of all fields */
  Uint f;/* dummy index */
  char status_str[100];
  int status;
  
  /* read order of fields to be solved from input */
  field_name = get_solving_field_name(SolveEqs->solving_order,&nf);
  
  /* solving fields in order */
  for (f = 0; f < nf; ++f)
  {
    printf("Verifying the Jacobian of equation for '%s' ...\n",field_name[f]);
    
    /* set the name of the field we are solving it */
    SolveEqs->field_name = field_name[f];
    
    /* set solving_man->cf */
    set_solving_man_cf(SolveEqs);
    
    /* set solving_man->settings */
    set_solving_man_settings(SolveEqs);
    
    /* picking up labeling, mapping etc. */
    preparing_ingredients(SolveEqs);

    /* set the name of the field we are solving it */
    SolveEqs->field_name = field_name[f];

    status = Jwritten_vs_Jequation(SolveEqs);
    
    if (status == TEST_SUCCESSFUL)
      sprintf(status_str,":)");
    else if (status == TEST_UNSUCCESSFUL)
      sprintf(status_str,":(");
    else
      Error0(NO_OPTION);
      
    printf("Verifying the Jacobian of equation for '%s' --> %s\n",field_name[f],status_str);
  }
  
  /* free */
  free_2d_mem(field_name,nf);
  free_solve_equations(SolveEqs);
}

/* compare J written by user and J computed from (eq(F+dF)-eq(F))/dF.
// ->return value: status which is TEST_UNSUCCESSFUL or TEST_SUCCESSFUL */
static int Jwritten_vs_Jequation(Solve_Equations_T *const SolveEqs)
{
  Grid_T *const grid = get_grid_solve_equations(SolveEqs);
  Matrix_T *J_Schur,*J_Reg;
  int status;
  
  /* if number grid only has one patch */
  if (grid->np == 1)
    Error0(INCOMPLETE_FUNC);
  
  J_Schur = making_J_Schur_Method(SolveEqs);
  J_Reg   = making_J_Old_Fashion( SolveEqs);
  
  status = compare_Js(grid,J_Reg,J_Schur);
  free_matrix(J_Reg);
  free_matrix(J_Schur);
  
  return status;
}

/* comparing entries of J_Schur and J_Reg
// ->return value: TEST_UNSUCCESSFUL or TEST_SUCCESSFUL */
static int compare_Js(Grid_T *const grid,const Matrix_T *const J_Reg,const Matrix_T *const J_Schur)
{
  const Uint dim = total_nodes_grid(grid);
  const double ERR = 1e-9;
  double **const J_s = J_Schur->reg->A;
  double **const J_r = J_Reg->reg->A;
  Flag_T flg = NONE;
  Uint i,j;
  
  for (i = 0; i < dim; ++i)
  {
    for (j = 0; j < dim; ++j)
      if (GRT(ABSd(J_s[i][j]-J_r[i][j]),ERR))
      {
        printf("J_Schur = %g, J_Reg = %g, diff = %g\n",
                  J_s[i][j],J_r[i][j],J_s[i][j]-J_r[i][j]);
                  
        flg = FOUND;
      }
  }
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;

}

/* free(Schur->f) free(Schur->g) */
static void free_schur_f_g(Grid_T *const grid)
{
  DDM_Schur_Complement_T *Schur;
  Uint p;
  
  /* free f */
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Schur = patch->solving_man->method->SchurC;
    free(Schur->f);
    free(Schur->g);
  }
}

/* making J by varying of F in Jx = -F.
// ->return value: J in regular order. */
static Matrix_T *making_J_Old_Fashion(Solve_Equations_T *const SolveEqs)
{
  Grid_T *const grid = get_grid_solve_equations(SolveEqs);
  const Uint dim = total_nodes_grid(grid);
  Matrix_T *J_Reg = alloc_matrix(REG_SF,dim,dim);
  double **const J = J_Reg->reg->A;
  const double CONST = 1.;
  const Uint npatch = grid->np;
   
  double *F1,*F2;
  Uint R;/* reference */
  Uint p,pn,ijk,df;
  
  /* making F1 = F(f) */
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    make_f(patch);
    make_partial_g(patch);
  }
  make_g(grid);
  
  F1 = make_col_F(grid);
  free_schur_f_g(grid);
  
  /* making F2 = F(f+df) */
  R = 0;
  for (pn = 0; pn < npatch; ++pn)
  {
    Patch_T *patch2 = grid->patch[pn];
    // NOTE: potentially might be a problem when we have an alias field
    Field_T *f = patch2->fields[LookUpField(SolveEqs->field_name,patch2)];
    double EPS = CONST/patch2->nn;
    
    for (df = 0; df < patch2->nn; ++df)
    {
      f->v[df] += EPS;
      free_coeffs(f);
      
      if (SolveEqs->FieldUpdate)/* if any FieldUpdate set */
        SolveEqs->FieldUpdate(patch2,SolveEqs->field_name);
      
      for (p = 0; p < npatch; ++p)
      {
        Patch_T *patch = grid->patch[p];
        make_f(patch);
        make_partial_g(patch);
      }
      make_g(grid);
      
      F2 = make_col_F(grid);
      free_schur_f_g(grid);
      
      f->v[df] -= EPS;
      free_coeffs(f);
      
      for (ijk = 0; ijk < dim; ++ijk)
        J[ijk][df+R] = (F2[ijk]-F1[ijk])/EPS;
      
      free(F2);
    }
    R += patch2->nn;
  }/* end of for (pn = 0; pn < npatch; ++pn) */
  
  free(F1);
  return J_Reg;
}

/* making the whole col F in Jx = -F in Newton method.
// ->return value: F */
static double *make_col_F(Grid_T *const grid)
{
  const Uint dim = total_nodes_grid(grid);
  double *F = alloc_double(dim);
  DDM_Schur_Complement_T *Schur;
  double *f,*g;
  Uint R;/* reference */
  const Uint *inv,*Iinv;
  Uint NS,NI;
  Uint p,i,j;
  
  R = 0;
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Schur = patch->solving_man->method->SchurC;
    NS = Schur->NS;
    NI = Schur->NI;
    inv = Schur->inv;
    Iinv = Schur->Iinv;
    f = Schur->f;
    g = Schur->g;
    
    for (i = 0; i < NS; ++i)
      F[R+inv[i]] = f[i];
    
    for (j = 0; j < NI; ++j)
      F[R+Iinv[j]] = g[j];
      
    R += patch->nn;
  }
  
  return F;
}

/* using Schur Complement method, makes J.
// also it returns it back to normal order not new labeled one which
// used in Schur complement method.
// ->return value: J in regular order. */
static Matrix_T *making_J_Schur_Method(Solve_Equations_T *const SolveEqs)
{
  Grid_T *const grid = get_grid_solve_equations(SolveEqs);
  const Uint dim = total_nodes_grid(grid);
  Matrix_T *J_Schur = alloc_matrix(REG_SF,dim,dim);
  double **const J = J_Schur->reg->A; 
  DDM_Schur_Complement_T *Schur;
  Matrix_T *B,*Et,**F,**C;
  const Uint npatch = grid->np;
  const Uint *inv,*Iinv,*NI_p;
  Uint NS,NI;
  Uint R;/* reference */
  Uint p;
  
  DDM_SCHUR_OpenMP(omp parallel for)
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    making_B_and_E(patch);
    making_F_and_C_Regular(patch);
  }
  
  R = 0;
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Uint i,j,k;
    
    Schur = patch->solving_man->method->SchurC;
    NS    = Schur->NS;
    NI    = Schur->NI;
    inv   = Schur->inv;
    Iinv  = Schur->Iinv;
    B     = Schur->B;
    Et    = Schur->E_Trans;
    F     = Schur->F;
    C     = Schur->C;
    NI_p  = Schur->NI_p;
    
    for (i = 0; i < NS; ++i)
      for (j = 0; j < NS; ++j)
        J[R+inv[i]][R+inv[j]] = B->reg->A[i][j];
        
    free_matrix(B);
    
    for (i = 0; i < NS; ++i)
      for (j = 0; j < NI; ++j)
        J[R+inv[i]][R+Iinv[j]] = Et->reg->A[j][i];/* E is transpose */
        
    free_matrix(Et);
    
    Uint r = 0;/* reference */
    for (k = 0; k < npatch; ++k)
    {
      DDM_Schur_Complement_T *Schur2;
      const Uint *Iinv2;
      Schur2 = grid->patch[k]->solving_man->method->SchurC;
      Iinv2 = Schur2->Iinv;
      
      if (F[k])
      {
        for (i = 0; i < NI_p[k]; ++i)
          for (j = 0; j < NS; ++j)
            J[r+Iinv2[i]][R+inv[j]] = F[k]->reg->A[i][j];
        
        free_matrix(F[k]);
      }
      if (C[k])
      {
        for (i = 0; i < NI_p[k]; ++i)
          for (j = 0; j < NI; ++j)
            J[r+Iinv2[i]][R+Iinv[j]] = C[k]->reg->A[i][j];
        
        free_matrix(C[k]);
      }
      r += grid->patch[k]->nn;
    }
    R += patch->nn;
  }
  
  return J_Schur;
}

/* given the grid and the name of the field being solved,
// calculate the residual of the equation on the grid 
// and save it in the the field whose name is save_name.
// it assumes the Schur struct and equation structs are filled. */
void calculate_equation_residual(Solve_Equations_T *const SolveEqs)
{
  FUNC_TIC
  
  Grid_T *grid;
  char **field_name = 0;/* name of all fields to be solved */
  const char *const fsuffix = SolveEqs->residual_suffix;
  Uint nf = 0;/* number of all fields */
  Uint f;/* dummy index */
  
  /* read order of fields to be solved from input */
  field_name = get_solving_field_name(SolveEqs->solving_order,&nf);
  
  /* synchronize pools with the default grid */
  sync_patch_pools(SolveEqs->grid,SolveEqs);
  
  /* solving fields in order */
  for (f = 0; f < nf; ++f)
  {
    /* set the name of the field we are solving it */
    SolveEqs->field_name = field_name[f];
    
    /* get the computational grid */
    grid = get_grid_solve_equations(SolveEqs);
    
    /* set solving_man->cf */
    set_solving_man_cf(SolveEqs);
    
    /* set solving_man->settings */
    set_solving_man_settings(SolveEqs);
    
    /* picking up labeling, mapping etc. */
    preparing_ingredients(SolveEqs);
    
    /* if number grid only has one patch */
    if (grid->np == 1)
    {
      Patch_T *patch = grid->patch[0];
      DDM_Schur_Complement_T *Schur = patch->solving_man->method->SchurC;
      char field_res[MSG_SIZE1];
      sprintf(field_res,"%s%s",field_name[f],fsuffix);
      empty_field(patch->fields[Ind(field_res)]);
      patch->fields[Ind(field_res)]->v = alloc_double(patch->nn);
        
      double *const res = patch->fields[Ind(field_res)]->v;
      const Uint NS = Schur->NS;
      const Uint *const inv = Schur->inv;
      Uint s,s_node;
      
      make_f(patch);/* making f */
      for (s = 0; s < NS; ++s)
      {
        s_node = inv[s];
        res[s_node] = Schur->f[s];
      }
      free(Schur->f);/* free{f} */
    }
    else/* multi-domains grid */
    {
      const Uint npatch = grid->np;
      Uint p;
      
      DDM_SCHUR_OpenMP(omp parallel for)
      for (p = 0; p < npatch; ++p)
      {
        Patch_T *patch = grid->patch[p];
        make_f(patch);
        make_partial_g(patch);
      }
      make_g(grid);/* free pg */
      
      DDM_SCHUR_OpenMP(omp parallel for)
      for (p = 0; p < npatch; ++p)
      {
        Patch_T *patch = grid->patch[p];
        DDM_Schur_Complement_T *Schur = patch->solving_man->method->SchurC;
        char field_res[MSG_SIZE1];
        sprintf(field_res,"%s%s",field_name[f],fsuffix);
        empty_field(patch->fields[Ind(field_res)]);
        patch->fields[Ind(field_res)]->v = alloc_double(patch->nn);
        
        double *const res = patch->fields[Ind(field_res)]->v;
        const Uint NS = Schur->NS;
        const Uint NI = Schur->NI;
        const Uint *const inv  = Schur->inv;
        const Uint *const Iinv = Schur->Iinv;
        Uint s,s_node,i,i_node;
        
        for (s = 0; s < NS; ++s)
        {
          s_node      = inv[s];
          res[s_node] = Schur->f[s];
        }
        for (i = 0; i < NI; ++i)
        {
          i_node      = Iinv[i];
          res[i_node] = Schur->g[i];
        }
      }
      free_schur_f_g(grid);/* free {f,g} */
    }
    sync_patch_pools(grid,SolveEqs);
  }/* end of for (f = 0; f < nf; ++f) */
  
  /* free names */
  free_2d_mem(field_name,nf);
  
  FUNC_TOC
}

/* print an introcuction to Schur complement method */
static void pr_intro_ddm_schur_complement(void)
{
  pr_line_custom('*');
  const char *const pr_message =
  
  "\n* Quick Introduction to Schur Complement Domain Decomposition Method *\n\n"

  "  We want to solve the system of equations A z = b on the given grid,\n"
  "for example, Jx = -F equations in the Newton-Raphson method.\n"
  "The system of equations A z = b are reordered as follows:\n\n"


  "               --  --\n"
  "              | B  E |  |x|   |f| -> x's are inner points of a patch    \n"
  "  A z = b =>  |      |  | | = | |                                       \n"
  "              | F  C |  |y|   |g| -> y's are interface points of a patch\n"
  "               --  --\n"
  "\n"   
  "          =>  B x + E y = f\n"
  "          =>  F x + C y = g\n"
  "\n"
  "          =>                        x = B^(-1) * ( f - E y )\n"
  "          =>  ( C - F * B^(-1)* E ) y = g - F * B^(-1) * f\n"
  "\n"
  "  Algorithm:\n"
  "\n" 
  "  1. Solve:   BE'= E and Bf'= f where E' = B^(-1)*E and f' = B^(-1)*f\n"
  "  2. Compute: g' = g - Ff'\n"
  "  3. Compute: S  = C - FE'\n"
  "  4. Solve:   Sy = g'\n"
  "  5. Compute: x  = f'-E'y\n";

  printf("%s\n",pr_message);
  pr_line_custom('*');
}

/* calloc one sewing.
// ->return value = memory for 1 sewing struct. */
Sewing_T *alloc_sewing(void)
{
  Sewing_T *sewing = calloc(1, sizeof(*sewing));
  IsNull(sewing);
  
  return sewing;
}

/* free thoroughly patch->interface */
void free_patch_SolMan_method_Schur(Patch_T *const patch)
{
  if (!patch->solving_man)
    return;
  //if (!patch->solving_man->method)
  //  return;
  
  DDM_Schur_Complement_T *s = patch->solving_man->method->SchurC;
  if (!s)
    return;
    
  Sewing_T **se             = s->sewing;
  Pair_T **p                = 0;
  Uint i,j;
  
  /* note, those matrices and double populated during solver
  // won't needed to be freed */
  Free(s->map);
  Free(s->inv);
  Free(s->Imap);
  Free(s->Iinv);
  Free(s->NS_p);
  Free(s->NI_p);
  
  for (i = 0; i < s->nsewing; ++i)
  {
    p = 0;
    if (se[i])
    {
      if (se[i]->patchN != patch->pn)
      {
        Free(se[i]->map);
        Free(se[i]->inv);
        Free(se[i]->Imap);
        Free(se[i]->Iinv);
      }
      p = se[i]->pair;
      for (j = 0; j < se[i]->npair; ++j)
      {
        Free(p[j]->ip);
        Free(p[j]->nv);
        
        /* because only for the following patches we allocate memory */
        if (se[i]->patchN != patch->pn)
        {
          Free(p[j]->subface->flags_str);
          Free(p[j]->subface->id);
          Free(p[j]->subface->adjid);
          Free(p[j]->subface);  
          p[j]->subface = 0;
        }
        free(p[j]);
      }
      Free(p);
    }
    Free(se[i]);
  }
  Free(se);
  Free(s);
  
  patch->solving_man->method->SchurC = 0;
}


