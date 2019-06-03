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
  NOT_READY,
  INUSE,
  FOUND,
  CLEAN,
  BRUTE_FORCE,
  FATAL,
  INITIALIZE,
  UP = 0,
  DOWN = 1,
  LEFT = 2,
  RIGHT = 3,
  BACK = 4,
  FRONT = 5,
  NS_T_CS,/* NS type in cubed spherical */
  SR_T_CS,/* surrounding type in cubed spherical */
  OT_T1_CS,/* outermost type1 in cubed spherical */
  OT_T2_CS,/* outermost type2 in cubed spherical */
  NOT_INITIALIZE
}Flag_T;

/* collocation */
typedef enum COLLOCATION_T
{
  UNDEFINED_COLLOCATION = 0,
  EquiSpaced,
  Chebyshev_Extrema,
  Chebyshev_Nodes
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
  Spherical,
  ProjectiveHemisphereUp,
  ProjectiveHemisphereDown,
  StereographicSphereLeft,
  StereographicSphereRight,
  CubedSpherical
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
  CCS_L_SF/* long compressed column storage format */,
  CRS_L_SF/* long compressed row storage format */,
  TRI_L_SF/* long triplet storage format */,
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
  char *lv;/* letf value its name*/
  char *rv;/* right value string */
  double rv_double;/* right value double */
  double *rv_array;/* right value array */
  unsigned rv_n;/* right value unsigned or dimension of field */
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
  unsigned adjPatch;/* adjacent patch used in interpolation default is UINT_MAX */
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
  unsigned adjsn  ;/* adjacent subface number */
  unsigned np     ;/* number of points this surface has */
  unsigned *id    ;/* id of points this subface made of; it refers to node number*/
  unsigned *adjid ;/* id of adjacent point of each point, their index must be matched 
                   // e.g. adjacent point of id[ind1]=? is adjid[ind1]=? */
  unsigned face    ;/* the interface in which this point located */
  unsigned adjFace ;/* adjacent face used in interpolation or copy, for outerB or innerB they are UINT_MAX */
  unsigned adjPatch;/* adjacent patch used in interpolation or copy for outerB or innerB they are UINT_MAX */
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
  double (*j)(struct PATCH_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);/* function for transformation */
  double (*dN0_dx)(struct PATCH_T *const patch,const double *const X);/* specifically calculating dN0/dx at X */
  double (*dN0_dy)(struct PATCH_T *const patch,const double *const X);/* specifically calculating dN0/dy at X */
  double (*dN0_dz)(struct PATCH_T *const patch,const double *const X);/* specifically calculating dN0/dz at X */
  double (*dN1_dx)(struct PATCH_T *const patch,const double *const X);/* specifically calculating dN1/dx at X */
  double (*dN1_dy)(struct PATCH_T *const patch,const double *const X);/* specifically calculating dN1/dy at X */
  double (*dN1_dz)(struct PATCH_T *const patch,const double *const X);/* specifically calculating dN1/dz at X */
  double (*dN2_dx)(struct PATCH_T *const patch,const double *const X);/* specifically calculating dN2/dx at X */
  double (*dN2_dy)(struct PATCH_T *const patch,const double *const X);/* specifically calculating dN2/dy at X */
  double (*dN2_dz)(struct PATCH_T *const patch,const double *const X);/* specifically calculating dN2/dz at X */
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
  unsigned tri_l_f: 1;/* 1 if long tripet storage format, 0 otherwise*/
  unsigned ccs_l_f: 1;/* 1 if long compressed column storage format, 0 otherwise */
  unsigned crs_l_f: 1;/* 1 if long compressed row storage format,0 otherwise */
  long row;
  long col;
  struct/* triplet storage format */
  {
    int *row;
    int *col;
    double *a;
  }tri[1];
  struct/* compressed column storage format */
  {
    int *Ap;
    int *Ai;
    double *Ax;
  }ccs[1];
  struct/* compressed row storage format */
  {
    int *Ap;
    int *Aj;
    double *Ax;
  }crs[1];
  struct/* regular storage format */
  {
    double **A;
  }reg[1];
  struct/* triplet storage format */
  {
    long *row;
    long *col;
    double *a;
  }tri_long[1];
  struct/* compressed column storage format */
  {
    long *Ap;
    long *Ai;
    double *Ax;
  }ccs_long[1];
  struct/* compressed row storage format */
  {
    long *Ap;
    long *Aj;
    double *Ax;
  }crs_long[1];
}Matrix_T;

/* equation solver */
typedef int fEquation_Solver_T(void *vp);

/* a general prototype to embrace various types of equations */
typedef void *fEquation_T(void *vp1,void *vp2);

/* elements of Jacobian for equations like dfxx_df etc. */
typedef double fJs_T(Matrix_T *const m,const long i,const long j);

/* equation stucture */
typedef struct sEQUATION_T
{
  char name[__1MAX_STR_LEN1__];
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
  Matrix_T *B;
  Matrix_T *E_Trans;/* NOE: this is TRANSPOSE of E */
  Matrix_T *E_Trans_prime;/* NOTE: it is E' of E_Trnas. */
  Matrix_T *F_by_E_prime;/* it is made in CCS format */
  Matrix_T **F;
  Matrix_T **C;
  Matrix_T *C_ccs;/* combining all of the C's into one CCS format matrix */
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
  double Frms;/* the residual(rms) of F in, Jx=-F for this field. 
              // note: it's initialized to DBL_MAX. */
  fEquation_T **field_eq;/* the equation needed to be satisfied */
  fEquation_T **bc_eq;/* the B.C. needed to be satisfied */
  fEquation_T **jacobian_field_eq;/* jacobian for field equations */
  fEquation_T **jacobian_bc_eq;/* jacobian for B.C. equations */
  struct/* jacobian elements */
  {
    char type[__1MAX_STR_LEN1__];
    Matrix_T *J;
  }**jacobian;
  unsigned nj;/* number of jacobian */
  
  struct/* various method to solve */
  {
    /* type of method */
    unsigned Schur_Complement: 1;/*1 if schur complement, 0 otherwise*/
    DDM_Schur_Complement_T *SchurC;
  }method[1];
  
}Solving_Man_T;

/* patch */
typedef struct PATCH_T
{
  struct GIRD_T *grid;/* refers to its grid */
  char *name;/* patch name */
  Coord_T coordsys;/* coord sys used in this patch */
  struct
  {
   struct
   {
    double R1;/* smaller R */
    double R2;/* bigger R */
    Field_T *R1_f;/* smaller R field */
    Field_T *R2_f;/* bigger R field */
    Field_T *dR1_dx;/* dR1/dx */
    Field_T *dR1_dy;/* dR1/dy */
    Field_T *dR1_dz;/* dR1/dz */
    Field_T *dR2_dx;/* dR2/dx */
    Field_T *dR2_dy;/* dR2/dy */
    Field_T *dR2_dz;/* dR2/dz */
   }ProjectiveCoord[1];
   struct
   {
    Flag_T side;/* the side of this cubed coord, up, down, etc. */
    Flag_T type;/* type of cubed spherical, NS, SR, OT */
    Field_T *R1_f;/* cubed spherical surface function (small). Note: it's always positive for all sides. */
    Field_T *R2_f;/* cubed spherical surface function (big). Note: it's always positive for all sides. */
    Field_T *dR1_dx;/* dR1/dx */
    Field_T *dR1_dy;/* dR1/dy */
    Field_T *dR1_dz;/* dR1/dz */
    Field_T *dR2_dx;/* dR2/dx */
    Field_T *dR2_dy;/* dR2/dy */
    Field_T *dR2_dz;/* dR2/dz */
    double R1;/* small radius of outermost patches. Note: it's always positive for all sides. */
    double R2;/* big radius of outermost patches. Note: it's always positive for all sides.*/
    double xc1;/* const xc value for those patches that have square (small)*/
    double xc2;/* const xc value for those patches that have square (big)*/
   }CubedSphericalCoord[1];
  }CoordSysInfo[1];
  Collocation_T collocation[3];/* type of collocation in each direction */
  Basis_T basis[3];/* the type of basis for functions used in this patch 
                   // each refers to basis in that specific direction.
                   // e.g. basis[2] = Chebyshev_FirstKind_BASIS means
                   // in c direction it uses that basis.
                   */
  unsigned n[3];/* number of points (nodes) in each direction */
  unsigned nn;/* number of nodes in this patch */
  unsigned pn;/* its patch number i.e. patch[pn] = patch */
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
  Solving_Man_T *solving_man;/* solving managing */
  unsigned innerB:1;/* if this patch has inner boundary 1 otherwise 0 */
  unsigned is_a_closed: 1;/* if coordinate a is periodic or closed 1, otherwise 0 */
  unsigned is_b_closed: 1;/* if coordinate b is periodic or closed 1, otherwise 0 */
  unsigned is_c_closed: 1;/* if coordinate c is periodic or closed 1, otherwise 0 */
}Patch_T;

/* grid */
typedef struct GIRD_T
{
  char *kind;/* type of grid which refers how we cover the grid */
  Patch_T **patch;/* covering patch */
  unsigned gn;/* grid number */
  unsigned np;/* number of patches on grid */
  unsigned nn;/* total number of nodes on grid */
}Grid_T;

/* jacobian for equations */
typedef struct JACOBIAN_EQ_T
{
  Field_T **field;
  unsigned nf;/* number of fields */
}Jacobian_Eq_T;


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

struct INTERPOLATION_T;

/* interpolation function typedef */
typedef double fInterpolation_T(struct INTERPOLATION_T *const interp_s);

/* interpolation struct used in interpolation function */
typedef struct INTERPOLATION_T
{
  Field_T *field;/* interesting field for interpolation */
  fInterpolation_T *interpolation_func;/* interpolation function */
  double X,Y,Z;/* where interpolant calculated. 
               // MUST be provided in coords sys. used by patch.
               */
  unsigned X_dir_flag   : 1;/* 1-D interpolation in X direction */
  unsigned Y_dir_flag   : 1;/* 1-D interpolation in Y direction */
  unsigned Z_dir_flag   : 1;/* 1-D interpolation in Z direction */
  unsigned XY_dir_flag  : 1;/* 2-D interpolation in X&Y direction */
  unsigned XZ_dir_flag  : 1;/* 2-D interpolation in X&Z direction */
  unsigned YZ_dir_flag  : 1;/* 2-D interpolation in Y&Z direction */
  unsigned XYZ_dir_flag : 1;/* 3-D interpolation in X&Y&Z direction */
  unsigned I;/* the index held constant in case of interpolation in 1-D and 2-D */
  unsigned J;/* the index held constant in case of interpolation in 1-D and 2-D */
  unsigned K;/* the index held constant in case of interpolation in 1-D and 2-D */
}Interpolation_T;

/* boundary condition struct */
typedef struct BOUNDARY_CONDITION_T
{
  Patch_T *patch;/* patch that has this boundary */
  SubFace_T *subface;/* the subface located at interesting boundary */
  Field_T *field;/* the field this B.C.s to be imposed */
  unsigned cn;/* collection number */
  unsigned *node;/* nodes's index at the boundary, i.e node[i] = node number used in the patch */
  unsigned nn;/* number of nodes */
}Boundary_Condition_T;

/* *******************************************
// structure for various tasks and projects
// *******************************************
*/
typedef struct TOV_PROJECT_T
{
 double m;/* NS mass */
}TOV_T;

/* umfpack direct solver */
typedef struct UMFPACK_T
{
  const char *description;
  Matrix_T *a;/* a in a.x = b */
  double *b;/* in ax=b */
  double *x;/* in ax=b */
  double **xs;/* x for series solving Ax[i]=b[i], i runs 0,...,ns */
  double **bs;/* b for series solving Ax[i]=b[i], i runs 0,...,ns */
  unsigned ns;/* number of series */
}UmfPack_T;
