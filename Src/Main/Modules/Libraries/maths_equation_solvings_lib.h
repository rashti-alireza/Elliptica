void make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_dfs_df_values(Grid_T *const grid);
void test_dInterp_a_df(Grid_T *const grid);
void *init_eq(void);
void add_eq(sEquation_T ***const data_base, fEquation_T *const eq,const char *const name);
void initialize_solving_man(Grid_T *const grid,sEquation_T **const field_eq,sEquation_T **const bc_eq,sEquation_T **const jacobian_field_eq, sEquation_T **const jacobian_bc_eq);
int solve_eqs(Grid_T *const grid);
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
double read_matrix_entry_ccs(Matrix_T *const m, const long r,const long c);
fJs_T *get_j_reader(const Matrix_T *const m);
void test_solve_ddm_schur_complement(Grid_T *const grid);


/* defining some macros to improve the readability and simplicity */

/* macros for jacobian of equations */
#define DDM_SCHUR_JACOBIAN_EQ_DECLARE \
  Patch_T *const patch  = vp1;\
  DDM_Schur_Complement_T *const S = vp2;\
  double **const B = S->B->reg->A;\
  double **E_Trans;\
  const unsigned *const node = S->inv;\
  const unsigned Ni = S->Oi;/* number of inner mesh nodes */\
  const unsigned Nj = S->NS;/* number of inner mesh+outer-boundary + inner-boundary nodes */\
  const unsigned K0 = S->NS;/* number of inner mesh+outer-boundary + inner-boundary nodes */\
  const unsigned Nk = patch->nn;/* total number of nodes */\
  unsigned i,j,k;

/* macro for B part of jacobian */
#define DDM_SCHUR_JACOBIAN_EQ_Bpart_OPEN \
  for (i = 0; i < Ni; ++i)\
  {\
    ijk = node[i];\
    for (j = 0; j < Nj; ++j)\
    {\
      lmn = node[j];\

#define DDM_SCHUR_JACOBIAN_EQ_Bpart_CLOSE \
    }/* end of for (i = 0; i < Ni; ++i) */\
  }/* end of for (j = 0; j < Nj; ++j) */

/* macros for E part of jacobian */
#define DDM_SCHUR_JACOBIAN_EQ_Epart_OPEN \
  if (S->NI)/* if there is any interface points then E is needed */\
  {\
    E_Trans = S->E_Trans->reg->A;\
    for (k = K0; k < Nk; ++k)\
    {\
      lmn = node[k];\
      j = k-K0;\
      for (i = 0; i < Ni; ++i)\
      {\
        ijk = node[i];

#define DDM_SCHUR_JACOBIAN_EQ_Epart_CLOSE \
     }/* end of for (i = 0; i < Ni; ++i) */\
    }/* end of for (k = K0; k < Nk; ++k) */\
  }/* end of if (S->NI) */


/* macros for jacobian of boundary condition equations */
#define DDM_SCHUR_JACOBIAN_BC_DECLARE \
  Patch_T *const patch  = vp1;\
  DDM_Schur_Complement_T *const S = vp2;\
  double **const B = S->B->reg->A;\
  double **E_Trans;\
  const unsigned *const node = S->inv;\
  const unsigned I0 = S->Oi;/* number of inner mesh nodes */\
  const unsigned Ni = S->NS;/* number of inner mesh+outer-boundary + inner-boundary nodes */\
  const unsigned Nj = S->NS;/* number of inner mesh+outer-boundary + inner-boundary nodes */\
  const unsigned K0 = S->NS;/* number of inner mesh+outer-boundary + inner-boundary nodes */\
  const unsigned Nk = patch->nn;/* total number of nodes */\
  unsigned i,j,k;

/* macro for B part of outer boundary jacobian */
#define DDM_SCHUR_JACOBIAN_BC_Bpart_OPEN \
  for (i = I0; i < Ni; ++i)\
  {\
    ijk = node[i];\
    for (j = 0; j < Nj; ++j)\
    {\
      lmn = node[j];\

#define DDM_SCHUR_JACOBIAN_BC_Bpart_CLOSE \
    }/* end of for (i = I0; i < Ni; ++i) */\
  }/* end of for (j = 0; j < Nj; ++j) */

/* macros for E part of jacobian */
#define DDM_SCHUR_JACOBIAN_BC_Epart_OPEN \
  if (S->NI)/* if there is any interface points then E is needed */\
  {\
    E_Trans = S->E_Trans->reg->A;\
    for (k = K0; k < Nk; ++k)\
    {\
      lmn = node[k];\
      j = k-K0;\
      for (i = I0; i < Ni; ++i)\
      {\
        ijk = node[i];

#define DDM_SCHUR_JACOBIAN_BC_Epart_CLOSE \
     }/* end of for (i = I0; i < Ni; ++i) */\
    }/* end of for (k = K0; k < Nk; ++k) */\
  }/* end of if (S->NI) */


/* macro for equation */
#define DDM_SCHUR_EQ_DECLARE \
  Patch_T *const patch = vp1;\
  DDM_Schur_Complement_T *const S = vp2;\
  double *const F = S->f;\
  const unsigned *const node  = S->inv;/* inverse map to node */\
  const unsigned N = S->Oi;/* number of inner mesh nodes */\
  unsigned n;
  
#define DDM_SCHUR_EQ_OPEN \
  for (n = 0; n < N; ++n)\
  {\
    ijk  = node[n];


#define DDM_SCHUR_EQ_CLOSE }

/* macro for boundary condition */
#define DDM_SCHUR_BC_DECLARE \
  Boundary_Condition_T *const bc = vp1;\
  DDM_Schur_Complement_T *const S = vp2;\
  double *const F      = S->f;\
  unsigned *const map  = S->map;\
  Patch_T *const patch = bc->patch;\
  const unsigned *const node = bc->node;/* nodes at boundary */\
  const unsigned N = bc->nn;/* number of nodes at boundary */\
  unsigned n;

#define DDM_SCHUR_BC_OPEN \
  for (n = 0; n < N; ++n)\
  {\
    ijk  = node[n];
    
#define DDM_SCHUR_BC_CLOSE }

