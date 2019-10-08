enum ROOT_FINDER_enum
{
  ROOT_FINDER_UNDEF/* undefined */,
  ROOT_FINDER_OK/* root was found successfully */,
  ROOT_FINDER_EXTREMA/* it stuck in an extrema */,
  ROOT_FINDER_MAX_ITER/* exceeds from maximum number of iteration */,
  ROOT_FINDER_NO_IMPROVEMENT/* it could not improve it more */,
  ROOT_FINDER_INTERRUPTED/* it was interrupted by a condition by user */,
  ROOT_FINDER_NAN/* the residual gets nan */
};

/* struct for root finder routine */
typedef struct ROOT_FINDER_T
{
  const char *type;/* type of root finder */
  const char *description;/* if might give some description for the root finder */
  double residual;/* residual of the function from zero */
  double tolerance;/* tolerance for f(x) = 0, 
                   // if |f^{iter+1}(x)-f^{iter}(x)| < tol, the root finder stops */
  unsigned n;/* number of variables (or equations) that make f = 0, 
             // e.g in {f1(x1,x2) = 0,f2(x1,x2) = 0, n is 2 */
  unsigned MaxIter;/* maximum iteration */
  const double *x_gss;/* initial guess */
  double *x_sol;/* solution of f(x) = 0 */
  unsigned FD_Left : 1;/* if 1 it uses finite difference with Left side stencil */
  unsigned FD_Right: 1;/* if 1 it uses finite difference with Right side stencil */
  void *params;/* parameters needed for evaluation of f(x) */ 
  /* f(x1,x2,...) = 0, params is supposed to refere to whatever is needed for evaluation of f */
  // note: since it might be systems of equations like {f1=0,f2=0,...} I used pointer to pointer function */
  double (**f)(void *params,const double *const x);
  /* df/dx^{dir}, params is the parameters are used for evalution of df_dx,
  // x is the dependent variables and dir is the direction of derivative */
  double (**df_dx)(void *params,const double *const x,const unsigned dir);
  double *(*root_finder_func)(struct ROOT_FINDER_T *const root);
  enum ROOT_FINDER_enum exit_status;/* exit status of root finder */
  int interrupt;/* if interrupt != 0, the root finder is interrupted.
                // for example, this controls if during search of root, 
                // root finder exceeds the domain of function. note, this 
                // must be set by the user at the equation function f(x). */
  unsigned verbose: 1;/* if 1, prints every step of root finding */
}Root_Finder_T;

/* solve equation struct that is passed to the solver.
// it may contain various functions and parameters to control and
// execute different tasks. */
typedef struct SOLVE_EQUATIONS_T
{
  Grid_T *grid;/* default grid, if no specific grid for particular
               // field has been specified, it uses this grid which
               // is given at the time of initialization. */
  const char *field_name;/* the name of the field that is being solved now */
  const char *solving_order;/* field name separated with comma to be solved,
                      // e.g. "phi,psi' means solve for first phi 
                      // and then psi. */
                      
  double relaxation_factor;/* in relaxation scheme we have :
                           // X_new = A*X'+(A-1)X_old, where
                           // A is the relaxation_factor, 
                           // X' is the solution found by the solver. 
                           // this factor can be set for each field separately
                           // and if no info is given, it is equal to 1, which
                           // means no relaxation. */
  /* some fields need their own grid, called sgrid (Special GRID) here. 
  // e.g. phi in Euler's equations is solved only in NS not the whole grid */
  struct
  {
    char *name;/* name of the field with special grid, e.g. phi */
    Grid_T *sgrid;/* e.g. the grid composed of NS patches for phi */
  }**Sgrid;/* the end of this struct determined by Null */
  
  /* instructions for updating field and its derivative according to the
  // field name particulare task for updating is done. note, if 
  // it has not been assigned it won't be execute. */
  void (*FieldUpdate)(Patch_T *const patch,const char *const name);
  
  /* instructions for updating the sources after the field has been 
  // solved on the whole grid and its derivative according to the given
  // name of the field. 
  // note, if it has not been assigned it won't be executed.*/
  void (*SourceUpdate)(Grid_T *const grid,const char *const name);
  
  /* this is the function specifies the stop criteria of the solver
  // if 1 it means continue, 0 means stop. if no function defined 
  // the default function is made using 
  // Solving_Residual and Solving_Max_Number_of_Solver_Step parameter */
  int (*StopCriteria)(Grid_T *const grid,const char *const name);
  
}Solve_Equations_T;

void print_root_finder_exit_status(const Root_Finder_T *const root);
Root_Finder_T *init_root_finder(const unsigned n);
double *execute_root_finder(Root_Finder_T *const root);
void plan_root_finder(Root_Finder_T *const root);
void free_root_finder(Root_Finder_T *root);
int solve_eqs(Solve_Equations_T *const SolveEqs);
void free_solve_equations(Solve_Equations_T *solve);
double get_relaxation_factor_solve_equations(Solve_Equations_T *const solve);
Solve_Equations_T *init_solve_equations(Grid_T *const grid);
Grid_T *get_grid_solve_equations(Solve_Equations_T *const solve);
void add_special_grid_solve_equations(Grid_T *const grid,const char *const name, Solve_Equations_T *const solve);
void make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_dfs_df_values(Grid_T *const grid);
void test_dInterp_a_df(Grid_T *const grid);
void *init_eq(void);
void add_eq(sEquation_T ***const data_base, fEquation_T *const eq,const char *const name);
void initialize_solving_man(Grid_T *const grid,sEquation_T **const field_eq,sEquation_T **const bc_eq,sEquation_T **const jacobian_field_eq, sEquation_T **const jacobian_bc_eq);
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
double read_matrix_entry_ccs(Matrix_T *const m, const long r,const long c);
fJs_T *get_j_reader(const Matrix_T *const m);
void test_solve_ddm_schur_complement(Grid_T *const grid);
void test_root_finders(Grid_T *const grid);

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

