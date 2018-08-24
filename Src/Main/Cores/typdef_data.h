/*
// Alireza Rashti
// May 2018
*/

/* *******************************************
// constants:
// *******************************************
*/

#define __1MAX_STR_LEN1__ 400

/* *******************************************
// enum:
// =====
//
// NOTE: DON'T CHANGE THE INITIALIZATION NUMBER 
// OF ENUMS. THESE CONVENTIONS ARE USED THROUGHOUT
// THE CODE. CHANGING THEM MIGHT LEAD TO SEGMENTAION
// FALUT.
// *******************************************
*/

/* common flags*/
typedef enum FLAG_T
{
  UNDEFINED = 0,
  NONE,
  NO,
  YES,
  READY,
  INUSE,
  FOUND,
  CLEAN,
  BRUTE_FORCE,
  FATAL
}Flag_T;

/* collocation */
typedef enum COLLOCATION_T
{
  UNDEFINED_COLLOCATION = 0,
  EquiSpaced,
  Chebyshev_Extrema
}Collocation_T;

/* types of basis enum */
typedef enum BASIS_T
{
  UNDEFINED_BASIS = 0/* undefined basis  */,
  No_BASIS/* when no basis is used */,
  Chebyshev_Tn_BASIS/* first kind Chebyshev basis T_n*/
}Basis_T;

/* coordinate system */
typedef enum COORD_T
{
  UNDEFINED_COORD = 0,
  Cartesian,
  CubedSphere
}Coord_T;

/* print flags */
typedef enum PRINT_T
{
  UNDEFINED_PRINT,
  PRINT_PARAMETERS,
  PRINT_COORDS,
  PRINT_INTERFACES
}Print_T;

/* face (interface) number */
typedef enum FACE_T
{
  I_0 = 0,
  I_n0,
  J_0,
  J_n1,
  K_0,
  K_n2,
  TOT_FACE = 6/* it's assumed we have 6 faces for each patch.
              // one might change this number with caveat.
              */
}Face_T;

/* various directions for derivative and Jacobian transformation
// NOTE: DON'T CHANGE THE NUMERATION.
*/
typedef enum DD_T
{
  _N0_ = 0/* for normalized 0-coord [-1,1] */,
  _N1_/* for normalized 1-coord [-1,1] */,
  _N2_ /* for normalized 2-coord [-1,1] */,
  _x_/* for Carteisian 0-coord */,
  _y_/* for Carteisian 1-coord */,
  _z_/* for Carteisian 2-coord */,
  _a_/* for Curvilinear 0-coord */,
  _b_/* for Curvilinear 1-coord */,
  _c_/* for Curvilinear 2-coord */,
  N_DD_T/* number of total DD_T*/,
  UNDEFINED_DIR/* for referring to undefined Dd_T */,
  _3_ = 3/* number three for tracking dimensions for different quantities */
 
}Dd_T;

/* matrix storage format */
typedef enum MATRIX_SF_T
{
  REG_SF/* regular storage format */,
  CCS_SF/* compressed column storage format */,
  CRS_SF/* compressed row storage format */,
  TRI_SF/* triplet storage format */,
  UNDEF_SF/* undefined */
}Matrix_SF_T;

/* *******************************************
// parameter:
// *******************************************
*/

/* parameters */
typedef struct PARAMETER_T
{
  /* syntax is expected to be lv = rv */
  char *lv;/* letf value */
  char *rv;/* right value */
}Parameter_T;

/* *******************************************
// project:
// *******************************************
*/

/* functions */
typedef int ProjFunc(void);

/* projects */
typedef struct PROJECT_T
{
  char *name;/* name of project */
  char *des;/* description */
  ProjFunc *func;/* project function */
}Project_T;

/* *******************************************
// grid, patch, solve, equations etc.
// *******************************************
*/

/* node*/
typedef struct NODE_T
{
  double x[3];/* for Cartesian value x,y,z */
  double *X;/* for general curvilinear value a,b,c */
}Node_T;

/* point */
typedef struct POINT_T
{
  unsigned ind;/* linear index in which node[ind] refers to this point */
  double N[3];/* normal vector */
  double *x;/* points to some triplet values for coords of a point */
  struct PATCH_T *patch;/* refers to its patch */
  struct POINT_T *adjPoint;/* adjacent point */
  unsigned face    ;/* the interface in which this point located */
  unsigned adjFace ;/* adjacent face used in interpolation */
  unsigned adjPatch;/* adjacent patch used in interpolation */
  unsigned sameX  : 1;/* 1 if addjacent face is on X = const */
  unsigned sameY  : 1;/* 1 if addjacent face is on Y = const */
  unsigned sameZ  : 1;/* 1 if addjacent face is on Z = const */
  unsigned touch  : 1;/* touch state 1, overlap state 0 */
  unsigned copy   : 1;/* copy state 1, interpolation state 0 */
  unsigned exterF : 1;/* external interface 1, internal 0
                      // external means it can reach other interface */
  unsigned outerB : 1;/* if it is outer boundary of grid 
                      // and needs boundary condition */
  unsigned innerB : 1;/* if it is inner boundary of grid 
                      // and needs boundary condition */
  unsigned houseK : 1;/* house-keeping purposes like if it is already
                      // counted 1 otherwise 0 among others */
}Point_T;

/* a subset of point on an interface */
typedef struct SUBFACE_T
{
  struct PATCH_T *patch;/* refers to its patch */
  char *flags_str ;/* encodes all of flags info in string format */
  unsigned sn     ;/* its subface number */
  unsigned adjsn  ;/* adjacent surface number */
  unsigned np     ;/* number of points this surface has */
  unsigned *id    ;/* id of points this subface made of; it refers to node number*/
  unsigned *adjid ;/* id of adjacent point of each point, their index must be matched 
                   // e.g. adjacent point of id[ind1]=? is adjid[ind1]=? */
  unsigned face    ;/* the interface in which this point located */
  unsigned adjFace ;/* adjacent face used in interpolation or copy */
  unsigned adjPatch;/* adjacent patch used in interpolation or copy */
  unsigned df_dn  : 1;/* 1 if d(field)/dn is set at interface, 0 otherwise */
  unsigned sameX  : 1;/* 1 if addjacent face is on X = const */
  unsigned sameY  : 1;/* 1 if addjacent face is on Y = const */
  unsigned sameZ  : 1;/* 1 if addjacent face is on Z = const */
  unsigned touch  : 1;/* touch state 1, overlap state 0 */
  unsigned copy   : 1;/* copy state 1, interpolation state 0 */
  unsigned exterF : 1;/* external interface 1, internal 0
                      // external means it can reach other interfaces */
  unsigned outerB : 1;/* if it is outer boundary of grid 
                      // and needs boundary condition */
  unsigned innerB : 1;/* if it is inner boundary of grid 
                      // and needs boundary condition */
}SubFace_T;

/* interface (face) */
typedef struct Interface_T
{
  struct PATCH_T *patch;/* refers to its patch */
  unsigned np;/* number of points in Point_T */
  unsigned fn;/* its interface number */
  unsigned ns;/* number of subfaces */
  Point_T **point;/* points on the interface */
  SubFace_T **subface;/* subset of points on this interface with same flags */
}Interface_T;

/* Jacobian transformation between different coords  */
typedef struct JACOBIAN_TRANS_T
{
  double (*j)(const struct PATCH_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);/* function for transformation */
  double *dX_dx[3][3];/* saving some transformation to save time for dX[0..2]/dx[0..2] */
  double *dx_dX[3][3];/* saving some transformation to save time dx[0..2]/dX[0..2] */
}JacobianTrans_T;

/* Field structure */
typedef struct FIELD_T
{
  char *name;/* name of field */
  double *v;/* values on each node on patch */
  double *v2;/* other values. if this field has two kinds of value:
             // e.g. fields in spectral expansion needs both 
             // coefficients of expansion and values of field on each node.
             // so for field with expansion this v2 refers to coefficents
             // value.
             */
  char *attr;/* attributes of fields like its dimension 
             // or other essential info.
             */
  char *info;/* each field might need more info or attributes 
             // which will save here and they are temporary.
             // e.g. the info about the coeffs of field. this info 
             // is dynamic.
             */
  struct PATCH_T *patch;/* refers to its patch which this field defined */
  void *v_ptr;/* refers to various quantities based on need */
}Field_T;

/* matrix */
typedef struct MATRIX_T
{
  unsigned reg_f: 1;/* regular format */
  unsigned tri_f: 1;/* 1 if tripet storage format, 0 otherwise*/
  unsigned ccs_f: 1;/* 1 if compressed column storage format, 0 otherwise */
  unsigned crs_f: 1;/* 1 if compressed row storage format,0 otherwise */
  long unsigned row;
  long unsigned col;
  struct/* triplet storage format */
  {
    long *row;
    long *col;
    double *a;
  }tri[1];
  struct/* compressed column storage format */
  {
    long *Ap;
    long *Ai;
    double *Ax;
  }ccs[1];
  struct/* compressed row storage format */
  {
    long *Ap;
    long *Aj;
    double *Ax;
  }crs[1];
  struct/* regular storage format */
  {
    double **A;
  }reg[1];
}Matrix_T;

/* jacobian for equations */
typedef struct JACOBIAN_EQ_T
{
  Matrix_T *matrix;
}Jacobian_Eq_T;

/* function for making Jacobian for equations */
typedef void fJacobian_Eq_T(Jacobian_Eq_T *const jac);

/* equation solver */
typedef int fEquation_Solver_T(void *vp);

/* equation that field obeys */
typedef void fEquation_T(void *vp,double *const F);

/* elements of Jacobian for equations like J_xx etc. */
typedef double fJs_T(Matrix_T *const m,const unsigned i,const unsigned j);

/* solve collects the entites needed to solve fields 
// based on eq's and b.c.'s and etc.
*/
typedef struct SOLVE_T
{
  char **f_name;/* f_name[#]="name of field #" */
  unsigned nf;/* number of fields */
  unsigned *f_occupy;/* refers to starting memory in 
                     // which the field number f occupies;
                     // so b[f_occupy[f]] to b[f_occupy[f]+nn] 
                     // are all memories in b which are used by 
                     // field with id f, nn is total nodes of a patch
                     */ 
  double *b;/* in  a.x = b */
  double *x;/* in  a.x = b */
  Matrix_T *a;/* in a.x = b */
  Field_T **field;/* fields to be found to satisfy equations and B.C.s */
  fEquation_T **field_eq;/* the equation needed to be satisfied */
  fEquation_T **bc_eq;/* the B.C needed to be satisfied */
  fJacobian_Eq_T **jac_eq;/* jacobian for equations */
  fEquation_Solver_T *solver;/* solver for ax = b */
}Solve_T;

/* equation stucture */
typedef struct sEQUATION_T
{
  char name[__1MAX_STR_LEN1__];
  fEquation_T *eq;/* the equation needed to be satisfied */
}sEquation_T;

/* solution managing */
typedef struct SOLUTION_MAN_T
{
  char **field_name;/* field to be solved */
  unsigned nf;/* number of fields */
  unsigned ns;/* number of solutions structure */
  Solve_T **solve;/* group fields to be solved together
                  // and in order appeared in input file
                  // based on eq's & b.c.'s.
                  // so seeking_field[0] is the first
                  // group of fields to be solved
                  // seeking_field[1] is the next and so on.
                  */
  struct/* jacobian elements */
  {
    char type[__1MAX_STR_LEN1__];
    Matrix_T *J;
  }**jacobian;
  unsigned nj;/* number of jacobian */
  fJs_T *j_func;
}Solution_Man_T;

/* patch */
typedef struct PATCH_T
{
  struct GIRD_T *grid;/* refers to its grid */
  char *name;/* patch name */
  Coord_T coordsys;/* coord sys used in this patch */
  Collocation_T collocation[3];/* type of collocation in each direction */
  Basis_T basis[3];/* the type of basis for functions used in this patch 
                   // each refers to basis in that specific direction.
                   // e.g. basis[2] = Chebyshev_FirstKind_BASIS means
                   // in c direction it uses that basis.
                   */
  unsigned n[3];/* number of points (nodes) in each direction */
  unsigned nn;/* number of nodes in this patch */
  unsigned pn;/* its patch number i.e. patch[pn] = patch */
  unsigned nc;/* node counter, sum of all nodes in previous patches */
  unsigned nfld;/* number of fields */
  double c[3];/* center */
  double s[3];/* size like length, width and height */
  double min[3];/* minimum of each direction like x_min = min[0] */
  double max[3];/* maximum of each direction like b_max = max[1] */
  Node_T **node;/* node info */
  Interface_T **interface;/* interface info  */
  JacobianTrans_T *JacobianT;/* Jacobian transformation between the coords */
  Field_T **pool;/* pool of fields, 
                 // notation: pool[Ind("Phi_f")] refers 
                 // to field phi. note: Ind is macro.
                 // one can access to values of field on this 
                 // patch like pool[Ind("phi_f")]->v */
  Solution_Man_T *solution_man;/* solution managing */
  unsigned innerB:1;/* if this patch has inner boundary 1 otherwise 0 */
}Patch_T;

/* grid */
typedef struct GIRD_T
{
  char *kind;/* type of grid which refers how we cover the grid */
  Flag_T status;/* INUSE or READY */
  Patch_T **patch;/* covering patch */
  unsigned np;/* number of patches on grid */
  unsigned nn;/* total number of nodes on grid */
}Grid_T;

/* *******************************************
// some typedef are common and used in different functions 
// *******************************************
*/

/* needle which should be found.
// it is a general purpose structure could be used by different functions
// like,point_finder function
*/
typedef struct NEEDLE_T
{
  const double *x;
  Grid_T *grid;/* the grid which is used */
  unsigned *guess;/* these are guess patches searched firstly */
  unsigned *in;/* force it to look only inside these patches. */
          /* notation: in[?] = patch number */
  unsigned *ex;/* force it to not look inside these patches */
           /* notation: ex[?] = patch number */
  unsigned *ans;/* the answers found which pointing to patch number */
  unsigned Nans;/* number of answer */
  unsigned Nin;/* number of included patches */
  unsigned Nex;/* number of excluded patches */
  unsigned Ng;/* numbef of guess patches */
}Needle_T;

/* print field */
typedef struct PR_FIELD_T
{
  const Grid_T *grid;
  const Patch_T *patch;
  const char *par;
  const char *folder;
  int cycle;
  double time;
  unsigned ng;/* number of group */
  void *group;/* points to a group for printing */
  void *opt_patch;/* points to options for patch */
  void *opt_field;/* points to options for field */
  void *vptr;/* general pointer for different purposes */
  void *a;/* a in double or float */
  void *b;/* b in double or float */
  void *c;/* c in double or float */
  void *v;/* v in double or float */
  void *file;/* file */
  void *file2;/* file */
}Pr_Field_T;

/* interpolation struct used in interpolation function */
typedef struct INTERPOLATION_T
{
  Field_T *field;/* interesting field for interpolation */
  Point_T *point;/* info about the point that interpolation takes place */
}Interpolation_T;

/* boundary condition struct */
typedef struct BOUNDARY_CONDITION_T
{
  SubFace_T *subface;/* the subface located at interesting boundary */
  Field_T *field;/* the field this B.C.s to be imposed */
  Solve_T *solve;/* used in interpolation and copy among others */
  unsigned cn;/* collection number */
  unsigned *node;/* nodes's index at the boundary, i.e node[i] = node number used in the patch */
  unsigned nn;/* number of nodes */
}Boundary_Condition_T;

/* *******************************************
// structures for solvers
// *******************************************
*/

/* umfpack direct solver */
typedef struct UMFPACK_T
{
  int row;/* number of row in matrix A */
  int col;/* number of column in matrix A */
  /* Ap, Ai and Ax refer to compressed column storage format for matrix A */
  long *Ap_long;
  long *Ai_long;
  int  *Ap;
  int  *Ai;
  double *Ax;
  double *b;/* in Ax=b */
  double *x;/* in Ax=b */
}UmfPack_T;