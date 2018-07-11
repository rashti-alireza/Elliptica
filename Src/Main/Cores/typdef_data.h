/*
// Alireza Rashti
// May 2018
*/

/* *******************************************
// enum:
// *******************************************
*/

/* common flags*/
typedef enum FLAG_T
{
  YES,
  NO,
  NONE,
  READY,
  INUSE,
  FOUND,
  BRUTE_FORCE
}Flag_T;

/* collocation */
typedef enum COLLOCATION_T
{
  EquiSpaced,
  Chebyshev_Zero,
  NDEF_COLLOCATION
}Collocation_T;

/* coordinate system */
typedef enum COORD_T
{
  Cartesian,
  CubedSphere
}Coord_T;

/* print flags */
typedef enum PRINT_T
{
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
  TOT_FACE = 6
}Face_T;

/* types of basis enum */
typedef enum BASIS_T
{
  No_Basis = 0,/* when no basis is used */
  Chebyshev_FirstKind_Basis/* first kind Chebyshev basis */
}Basis_T;

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
// grid and related: 
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
  unsigned adjFace ;/* adjacent face used in interpolation */
  unsigned adjPatch;/* adjacent patch used in interpolation */
  unsigned Dn_Df  : 1;/* 1 if Dn_Dfield is set at interface, 0 otherwise */
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

/* patch */
typedef struct PATCH_T
{
  struct GIRD_T *grid;/* refers to its grid */
  char *name;/* box name */
  char *coordsys;/* coord sys used in this patch */
  Collocation_T collocation;/* type of collocation in this patch */
  Basis_T basis;/* the type of basis for functions used in this patch */
  unsigned n[3];/* number of point in each direction */
  unsigned pn;/* its patch number i.e. patch[pn] = patch */
  double c[3];/* center */
  double s[3];/* size like length, width and height */
  double min[3];/* minimum of each direction like x_min = min[0] */
  double max[3];/* maximum of each direction like b_max = max[1] */
  Node_T **node;/* node info */
  Interface_T **interface;/* interface info  */
  unsigned innerB:1;/* if this patch has inner boundary 1 otherwise 0 */
}Patch_T;

/* field */
typedef struct FIELD_T
{
  char *name;/* its name like alpha or psi */
  double *value;/* its value on each grid point */
  double *coeffs;/* coefficients of basis if needed */
}Field_T;

/* grid */
typedef struct GIRD_T
{
  char *kind;/* type of grid which refers how we cover the grid */
  Flag_T status;/* INUSE or READY */
  Patch_T **patch;/* covering patch */
  unsigned nn;/* total number of nodes on grid */
  Field_T  **field;/* fields */
  unsigned nf;/* number of fields */
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
  double *x;
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
