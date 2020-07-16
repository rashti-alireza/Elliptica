#ifndef maths_equation_solvings_LIB_H
#define maths_equation_solvings_LIB_H

#include "maths_linear_algebra_lib.h"

#define MAX_STR 400

/* forward declaration structures */
struct FIELD_T;
struct MATRIX_T;


/* a general prototype to embrace various types of equations */
typedef void *fEquation_T(void *vp1,void *vp2);

/* elements of Jacobian for equations like dfxx_df etc. */
typedef double fJs_T(struct MATRIX_T *const m,const long i,const long j);

/* equation stucture */
typedef struct sEQUATION_T
{
  char name[MAX_STR];
  fEquation_T *eq;/* the equation needed to be satisfied */
}sEquation_T;

/* different quantities giving info abour pairing used in Schur complement */
typedef struct PAIR_T
{
  struct SEWING_T *sewing;/* refers to its sewing */
  double *pg;/* partial g that comping from this pair*/
  SubFace_T *subface;/* differet subfaces due to patch[pn] that
                     // is related to the current patch that equations
                     // are being set up 
                     */
  struct PAIR_T *mirror;/* the pair that is mirror of itself but
                        // from the other patch. */
  unsigned patchN;/* patch number which is equal to its sewing number */
  struct/* interpolation points;general coords of points
        // needed for interpolation subfaces */
  {
    double X[3];
  }*ip;
  struct/* normal vector at the subface */
  {
    double N[3];
  }*nv;
  
}Pair_T;

/* boundary information and how different patches are sown.
// this struct is made specially for having a good concurrency.
*/
typedef struct SEWING_T
{
  Pair_T **pair;
  unsigned patchN;/* patch number which is equal to its sewing number */
  unsigned npair;/* number of pairs */
  /* the following are the quantities that 
  // patch[patchN]->method->SchurC has.
  // it's used for purpose of concurrency and avoing race condition
  // bewteen pairs of different patches. more definition of each quantity
  // refer to SchurC strcut. */
  unsigned NS;
  unsigned NI;
  unsigned Oi;
  unsigned *map;
  unsigned *inv;
  unsigned *Imap;
  unsigned *Iinv;
}Sewing_T;

/* ingredients needed for mapping, labeling and etc for
// domain decomposition schur complement method
*/
typedef struct DDM_SCHUR_COMPLEMENT_T
{
  struct PATCH_T *patch;/* refers to its patch itself */
  /* regular means L(n,i,j,k) */
  unsigned *map;/* map: regular -> relabeled. ex: map[2] = 5 */
  unsigned *inv;/* inv: relabeled -> regular. ex: inv[5] = 2 */
  unsigned *Imap;/* interface point map, if it is given a point
                 // outside of its domain, it returns UINT_MAX. */
  unsigned *Iinv;/* interface point inverse map */
  unsigned NS;/* Number of subdomain points i.e. 
              // number of inner points + outerboundar points (NO) =>
              // total nodes - NS = number of interface points. Note:
              // outerboundary points excluded from interface points.
              */
  unsigned NI;/* total number of interface points, if 0, it means there
              // is no interface for this patch, for example when you
              // only have one single patch, all sides of the patch
              // are outerbounday so no interface with other patches. */
  unsigned Oi;/* initial index of outer boundary points at new label.
              // e.g. if NS = 10 and the last 3 points are 
              // outer boundary points then Oi = 7. 
              // furthermore, if there is no any outer boundary points 
              // then Oi = NS. */
  
/* namings:
   |B E||x|   |f|
   |F C||y| = |g|
*/
  double *f;
  double *f_prime;
  double *F_by_f_prime;
  double *g;
  double *x;
  double *y;
  struct MATRIX_T *B;
  struct MATRIX_T *E_Trans;/* NOE: this is TRANSPOSE of E */
  struct MATRIX_T *E_Trans_prime;/* NOTE: it is E' of E_Trnas. */
  struct MATRIX_T *F_by_E_prime;/* it is made in CCS format */
  struct MATRIX_T **F;
  struct MATRIX_T **C;
  struct MATRIX_T *C_ccs;/* combining all of the C's into one CCS format matrix */
  Sewing_T **sewing;/* sewing[patch_number] */
  unsigned nsewing;/* number of sewings which is = number of patches */
  unsigned np;/* total number of patches */
  unsigned *NS_p;/* SchurC->NS for each patch p */
  unsigned NS_total;/* summation of all NS_p */
  unsigned *NI_p;/* SchurC->NI for each patch p */
  unsigned NI_total;/* summation of all NI_p */
  
}DDM_Schur_Complement_T;

/* solving management */
typedef struct SOLVING_MAN_T
{
  struct PATCH_T *patch;/* refers to its patch itself */
  char **field_name;/* field to be solved */
  unsigned nf;/* number of fields */
  unsigned cf;/* current field; index of the field the is being solved */
  double Frms;/* the current residual(rms) of F in, Jx=-F for this field 
              // at this patch. note: it's initialized to DBL_MAX. */
  fEquation_T **field_eq;/* the equation needed to be satisfied */
  fEquation_T **bc_eq;/* the B.C. needed to be satisfied */
  fEquation_T **jacobian_field_eq;/* jacobian for field equations */
  fEquation_T **jacobian_bc_eq;/* jacobian for B.C. equations */
  struct/* jacobian elements */
  {
    char type[MAX_STR];
    struct MATRIX_T *J;
  }**jacobian;
  unsigned nj;/* number of jacobian */
  
  struct/* various method to solve */
  {
    /* type of method */
    unsigned Schur_Complement: 1;/*1 if schur complement, 0 otherwise*/
    DDM_Schur_Complement_T *SchurC;
  }method[1];
  
  struct/* settings and options for solver */
  {
    double relaxation_factor;/* relaxation factor in relaxation scheme */
    double Frms_i;/* the very beginning Frms (see Frms above for definition)
                 // which this field has at the its entrance to the solver. */
    double *HFrms;/* history of all Frms start form 0 to NFrms */
    double *last_sol;/* it is back up of last solution, 
                     // so in case the residula goes up, it uses this value. */
    unsigned NHFrms;/* number of HFrms */
    int solver_step;/* number of steps have been taken by solver till now. starting from 0 */
    int umfpack_size;/* 0 = di, otherwise dl (default is 0) */
    double umfpack_refine;/* max iter. refinement step, default is the default of UMFPACK which is 2 */
  }settings[1];
}Solving_Man_T;

/* equation solver */
typedef int fEquation_Solver_T(void *vp);

/* general function for variation of various kind of interpolation with 
// respect to the field. */
typedef double fdInterp_dfs_T(Patch_T *const patch,const double *const X,const unsigned df,const unsigned plane);


/* boundary condition struct */
typedef struct BOUNDARY_CONDITION_T
{
  Patch_T *patch;/* patch that has this boundary */
  SubFace_T *subface;/* the subface located at interesting boundary */
  struct FIELD_T *field;/* the field this B.C.s to be imposed */
  unsigned cn;/* collection number */
  unsigned *node;/* nodes's index at the boundary, i.e node[i] = node number used in the patch */
  unsigned nn;/* number of nodes */
}Boundary_Condition_T;


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
  unsigned eq_number;/* current equation number, there are cases ,e.g. PDE, 
                     // that the equations are the same but they are 
                     // at different point, this could help to populate the
                     // root->f with one function but the funcation is evaluated
                     // at different points. */
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
  double a_bisect;/* note: f(x) must change sign for x in [a,b]. */
  double b_bisect;/* note: f(x) must change sign for x in [a,b]. */
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
  
  int umfpack_size;/* (0 = di) otherwise long, default is 0 */
  double umfpack_refine;/* max iter. refinement step, default is the default of UMFPACK which is 2 */
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

void calculate_equation_residual(Solve_Equations_T *const SolveEqs);
char **get_solving_field_name(const char *const solving_order,unsigned *const nf);
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
struct MATRIX_T *get_j_matrix(const Patch_T *const patch,const char *type);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
double read_matrix_entry_ccs(struct MATRIX_T *const m, const long r,const long c);
fJs_T *get_j_reader(const struct MATRIX_T *const m);
void test_Jacobian_of_equations(Solve_Equations_T *const SolveEqs);
void test_root_finders(Grid_T *const grid);
fdInterp_dfs_T *get_dInterp_df(const Patch_T *const patch,const SubFace_T *const sf,const char *const dir);
Sewing_T *alloc_sewing(void);
void free_db_eqs(sEquation_T **db);
void free_patch_SolMan_jacobian(Patch_T *const patch);
void free_patch_SolMan_method_Schur(Patch_T *const patch);


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


#endif



