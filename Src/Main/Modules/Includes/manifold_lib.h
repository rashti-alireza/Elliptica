#ifndef manifold_LIB_H
#define manifold_LIB_H
#include "elliptica_system_lib.h"

/* prefix for patch name */
#define PATCH_NAME_P_      "grid"

/* prefix for patch name used in sprint */
#define PATCH_NAME_PRT_P_  PATCH_NAME_P_"%u_"

/* for now only two objects */
#define NPARAMS_GRID_CHAR (2)

/* to keep the previous notation these macros defined 
// with default precision factor 1. */
#define X_of_x(XXX,xxx,ppp) (X_of_x_precision((XXX),(xxx),(ppp),1.0))
#define x_of_X(xxx,XXX,ppp) (x_of_X_precision((xxx),(XXX),(ppp),1.0))

/* forward declaration structures */
struct FIELD_T;
struct SOLVING_MAN_T;
struct POINTSET_T;

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
  CubedSpherical
}Coord_T;

/* grid kind */
typedef enum GRID_KIND_T
{
  Grid_UNDEFINED = 0,
  Grid_SplitCubedSpherical_SNS,
  Grid_SplitCubedSpherical_SBH,
  Grid_SplitCubedSpherical_NSNS,
  Grid_SplitCubedSpherical_BHNS,
  Grid_SplitCubedSpherical_BHBH,
  Grid_CubedSpherical_BHNS,
  Grid_CubedSpherical_NSNS,
  Grid_Box
}Grid_Kind_T;

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
              // one might change this number with caveat. */
}Face_T;

/* various directions for derivative and Jacobian transformation
// NOTE: DON'T CHANGE THE NUMERATION.*/
typedef enum DD_T
{
  _N0_ = 0/* for normalized 0-coord [-1,1] */,
  _N1_/* for normalized 1-coord [-1,1] */,
  _N2_ /* for normalized 2-coord [-1,1] */,
  _x_/* for Carteisian 0-coord */,
  _y_/* for Carteisian 1-coord */,
  _z_/* for Carteisian 2-coord */,
  _a_/* for Curvilinear 0-coord (coordinate patch) */,
  _b_/* for Curvilinear 1-coord (coordinate patch) */,
  _c_/* for Curvilinear 2-coord (coordinate patch) */,
  N_DD_T/* number of total DD_T*/,
  UNDEFINED_DIR/* for referring to undefined Dd_T */,
  _3_ = 3/* number three for tracking dimensions for different quantities */
}Dd_T;

/* node*/
typedef struct NODE_T
{
  double x[3];/* for Cartesian value x,y,z */
  double *X;/* for general curvilinear value a,b,c 
            // it's the coordinate the patch uses. */
}Node_T;

/* point */
typedef struct POINT_T
{
  Uint ind;/* linear index in which node[ind] refers to this point */
  double N[3];/* normal vector */
  double *x;/* points to some triplet values for coords of a point */
  struct PATCH_T *patch;/* refers to its patch */
  struct POINT_T *adjPoint;/* adjacent point */
  Uint face    ;/* the interface in which this point located */
  Uint adjFace ;/* adjacent face used in interpolation */
  Uint adjPatch;/* adjacent patch used in interpolation default is UINT_MAX */
  Uint adjInd;/* adjacent point index correspond to adjPatch node[adjInd] */
  Uint adjIndF;/* adjacent point index correspond to adjFace point[adjIndF] */
  Uint IsOnEdge:1;/* if on edge 1 otherwise 0 */
  Uint sameX  : 1;/* 1 if addjacent face is on X = const */
  Uint sameY  : 1;/* 1 if addjacent face is on Y = const */
  Uint sameZ  : 1;/* 1 if addjacent face is on Z = const */
  Uint touch  : 1;/* touch state 1, overlap state 0 */
  Uint copy   : 1;/* copy state 1, interpolation state 0 */
  Uint exterF : 1;/* external interface 1, internal 0
                      // external means it can reach other interface */
  Uint outerB : 1;/* if it is outer boundary of grid 
                      // and needs boundary condition */
  Uint innerB : 1;/* if it is inner boundary of grid 
                      // and needs boundary condition */
  Uint houseK : 1;/* house-keeping purposes like if it is already
                      // counted 1 otherwise 0 among others */
}Point_T;

/* a subset of point on an interface */
typedef struct SUBFACE_T
{
  struct PATCH_T *patch;/* refers to its patch */
  char *flags_str ;/* encodes all of flags info in string format */
  Uint sn     ;/* its subface number */
  Uint adjsn  ;/* adjacent subface number */
  Uint np     ;/* number of points this surface has */
  Uint *id    ;/* id of points this subface made of; it refers to node number*/
  Uint *adjid ;/* id of adjacent point of each point, their index must be matched 
                   // e.g. adjacent point of id[ind1]=? is adjid[ind1]=? */
  Uint face    ;/* the interface in which this point located */
  Uint adjFace ;/* adjacent face used in interpolation or copy, for outerB or innerB they are UINT_MAX */
  Uint adjPatch;/* adjacent patch used in interpolation or copy for outerB or innerB they are UINT_MAX */
  Uint df_dn  : 1;/* 1 if d(field)/dn is set at interface, 0 otherwise */
  Uint sameX  : 1;/* 1 if addjacent face is on X = const */
  Uint sameY  : 1;/* 1 if addjacent face is on Y = const */
  Uint sameZ  : 1;/* 1 if addjacent face is on Z = const */
  Uint touch  : 1;/* touch state 1, overlap state 0 */
  Uint copy   : 1;/* copy state 1, interpolation state 0 */
  Uint exterF : 1;/* external interface 1, internal 0
                      // external means it can reach other interfaces */
  Uint outerB : 1;/* if it is outer boundary of grid 
                      // and needs boundary condition */
  Uint innerB : 1;/* if it is inner boundary of grid 
                      // and needs boundary condition */
  double precision_factor;/* precision factor for X_of_x_precision */
}SubFace_T;

/* interface (face) */
typedef struct INTERFACE_T
{
  struct PATCH_T *patch;/* refers to its patch */
  Uint np;/* number of points in Point_T */
  Uint fn;/* its interface number */
  Uint ns;/* number of subfaces */
  Point_T **point;/* points on the interface */
  SubFace_T **subface;/* subset of points on this interface with same flags */
  struct POINTSET_T **innerP;/* all points on the interface but edge */
  struct POINTSET_T **edgeP;/* all edge points on the interface  */
  double centerN[3];/* unit outward normal at the center of this face. */
  double centerx[3];/* x-coords of center of this face. */
  Uint df_dn:1;/* Drichlet = 0, Neumann = 1 */
  Uint df_dn_set:1;/* if df_dn is set = 1, otherwise 0. */
  Uint innerB:1;/* if has innerB is 1, otherwise 0. */
  Uint outerB:1;/* if has outerB is 1, otherwise 0. */
}Interface_T;

struct Collocation_s
{
  double min;
  double max;
  Uint n;
  double stp;
  double a;
  double b;
  Collocation_T c;
};

/* Jacobian transformation between different coords  */
typedef struct JACOBIAN_TRANS_T
{
  double (*j)(struct PATCH_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p);/* function for transformation */
  double *dX_dx[3][3];/* coords Jacobian for dX[0...2]/dx[0...2]|[ijk]. 
                        // notation: 0 = x/X, 1 = y/Y, and 2 = z/Z. */
  double *d2X_dxdy[3][6];/* coords Jacobian for d^2X[0...2]/dxdy[0...5]|[ijk].
                          // notation: 0 = x/X, 1 = y/Y, and 2 = z/Z.
                          // for [0...5] => xx, xy, xz, yy, yz, and zz. */
  double dN_dX[3];/* coords Jacobian for dN/dX, where N is the normalized coords, 
                  // i.e., [-1,1]. Since we ASSUMED dN/dX = constant then 
                  // for example dN_dX[1]= dN^1/dY. 
                  // notation: 0 = X, 1 = Y, 2 = Z. */
  /* split cubed spherical stuffs */
  struct
  {
   double sign;/* sign */
   Uint iper,jper,kper;/* permutation */
  }SCS[1];
  
}JacobianTrans_T;

/* patch */
typedef struct PATCH_T
{
  struct GRID_T *grid;/* refers to its grid */
  char *name;/* patch name */
  Coord_T coordsys;/* coord sys used in this patch */
  struct
  {
   char region[999];/* region this patch covers, parentheses separated */
   struct
   {
    Flag_T side;/* the side of this cubed coord, up, down, etc. */
    Flag_T type;/* type of cubed spherical, NS, SR, OT */
    struct FIELD_T *R1_f;/* cubed spherical surface function (small). Note: it's always positive for all sides. */
    struct FIELD_T *R2_f;/* cubed spherical surface function (big). Note: it's always positive for all sides. */
    struct FIELD_T *dR1_dx;/* dR1/dx */
    struct FIELD_T *dR1_dy;/* dR1/dy */
    struct FIELD_T *dR1_dz;/* dR1/dz */
    struct FIELD_T *dR2_dx;/* dR2/dx */
    struct FIELD_T *dR2_dy;/* dR2/dy */
    struct FIELD_T *dR2_dz;/* dR2/dz */
    double R1;/* small radius of outermost patches. Note: it's always positive for all sides. */
    double R2;/* big radius of outermost patches. Note: it's always positive for all sides.*/
    double xc1;/* const xc value for those patches that have cube (small)*/
    double xc2;/* const xc value for those patches that have cube (big)*/
   }CubedSphericalCoord[1];
  }CoordSysInfo[1];
  Collocation_T collocation[3];/* type of collocation in each direction */
  Basis_T basis[3];/* the type of basis for functions used in this patch 
                   // each refers to basis in that specific direction.
                   // e.g. basis[2] = Chebyshev_FirstKind_BASIS means
                   // in c direction it uses that basis.
                   */
  Uint n[3];/* number of points (nodes) in each direction */
  Uint nn;/* number of nodes in this patch */
  Uint nsplit[3];/* number of splits taken place for this type */
  Uint pn;/* its patch number i.e. patch[pn] = patch */
  Uint nfld;/* number of fields */
  double c[3];/* center */
  double s[3];/* size like length, width and height */
  double min[3];/* minimum of each direction like x_min = min[0] */
  double max[3];/* maximum of each direction like b_max = max[1] */
  Node_T **node;/* node info */
  Interface_T **interface;/* interface info  */
  JacobianTrans_T *JacobianT;/* Jacobian transformation between the coords */
  struct FIELD_T **fields;/* all fields for this patch, 
                 // notation: fields[Ind("Phi_f")] refers 
                 // to field phi. note: Ind is macro.
                 // one can access to values of field on this 
                 // patch like fields[Ind("phi_f")]->v */
  struct SOLVING_MAN_T *solving_man;/* solving management */
  Uint innerB:1;/* if this patch has inner boundary 1 otherwise 0 */
  Uint outerB:1;/* if this patch has outer boundary 1 otherwise 0 */
  Uint is_a_closed: 1;/* if coordinate a is periodic or closed 1, otherwise 0 */
  Uint is_b_closed: 1;/* if coordinate b is periodic or closed 1, otherwise 0 */
  Uint is_c_closed: 1;/* if coordinate c is periodic or closed 1, otherwise 0 */
}Patch_T;

/* grid */
typedef struct GRID_T
{
  Grid_Kind_T kind;/* type of grid which refers how we cover the grid */
  Patch_T **patch;/* covering patch */
  Uint gn;/* grid number */
  Uint np;/* number of patches on grid */
}Grid_T;

/* jacobian for equations */
typedef struct JACOBIAN_EQ_T
{
  struct FIELD_T **field;
  Uint nf;/* number of fields */
}Jacobian_Eq_T;

/* needle which should be found.
// it is a general purpose structure could be used by different functions
// like,point_finder function */
typedef struct NEEDLE_T
{
  const double *x;/* Cartesian coords */
  double precision_factor;/* precision factor for X_of_x_precision
                          // default = 1. . */
  Grid_T *grid;/* the grid which is used */
  Uint *guess;/* these are guess patches searched firstly */
  Uint *in;/* force it to look only inside these patches. */
          /* notation: in[?] = patch number */
  Uint *ex;/* force it to not look inside these patches */
           /* notation: ex[?] = patch number */
  Uint *ans;/* the answers found which pointing to patch number */
  Uint Nans;/* number of answer */
  Uint Nin;/* number of included patches */
  Uint Nex;/* number of excluded patches */
  Uint Ng;/* numbef of guess patches */
}Needle_T;

/* characteristics info of grid */
typedef struct GRID_CHAR_T
{
 Grid_T *grid;/* the new pristine grid */
 double S;/* separation between the objects in a binary system or 
          // the size of box around the single object. */
 struct/* for each object */
 {
  const char *obj;/* BH or NS */
  const char *dir;/* left or right (must be lower case) */
  double *relClm;/* Re(Clm) at Ylm expansion of surface function */
  double *imgClm;/* Im(Clm) at Ylm expansion of surface function */
  Uint lmax;/* lmax in Ylm expansion */
  double l,w,h;/* length(x), width(y) and hight(z) of the central cube */
  double r_min,r_max;/* min and max of surface function */
  Uint occupied:1;/* 1, if it occupied already, otherwise 0. */
 }params[NPARAMS_GRID_CHAR][1];
 
}Grid_Char_T;

double point_value(const Uint i, const struct Collocation_s *const coll_s);
void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const colloc,const Uint dir);
int make_patches(Grid_T *const grid);
int realize_interfaces(Grid_T *const grid);
int make_JacobianT(Grid_T *const grid);
int X_of_x_precision(double *const X,const double *const x,const Patch_T *const patch,const double precision_factor);
int x_of_X_precision(double *const x,const double *const X,const Patch_T *const patch,const double precision_factor);
double x_coord(const Uint i,const Patch_T *const patch);
double y_coord(const Uint i,const Patch_T *const patch);
double z_coord(const Uint i,const Patch_T *const patch);
double X_coord(const Uint i,const Patch_T *const patch);
double Y_coord(const Uint i,const Patch_T *const patch);
double Z_coord(const Uint i,const Patch_T *const patch);
double dq2_dq1(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p);
double JT_Cartesian_patch(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p);
double *normal_vec(Point_T *const point);
double General2ChebyshevExtrema(const double X,const Uint dir,const Patch_T *const patch);
void grid_characteristics_example(Grid_T *const grid);
double R_interpolation_CS(struct FIELD_T *const R,const double *const X);
void SignAndIndex_permutation_CubedSphere(const Flag_T side,Uint *const a,Uint *const b,Uint *const c,double *const s);
void test_CubedSpherical_Coordinates(Grid_T *const grid);
void test_dq_dN(Grid_T *const grid);
void needle_ex(Needle_T *const needle,const Patch_T *const patch);
void needle_in(Needle_T *const needle,const Patch_T *const patch);
void point_finder(Needle_T *const needle);
void populate_right_BH(Grid_T *const grid,const Uint pn);
void populate_right_BH_central_box(Grid_T *const grid,const Uint pn);
int make_nodes(Grid_T *const grid);
void alloc_nodes(Grid_T *const grid);
void alloc_interface(Patch_T *const patch);
void *alloc_point(const Uint s);
void free_points(Grid_T *const grid);
void free_patch_interface(Patch_T *const patch);
void free_needle(Needle_T *needle);
void *alloc_needle(void);
Patch_T make_temp_patch(const Patch_T *const patch);
void free_temp_patch(Patch_T *const patch);
void *alloc_grid(void);
void alloc_patches(Grid_T *const grid);
void free_grid(Grid_T *grid);
void free_patch(Patch_T *patch);
void free_grid_db(void);
void set_params_of_split_cubed_spherical_grid(Grid_Char_T *const grid_char);
void theta_phi_of_XY_CS(double *const theta,double *const phi,const double *const X,const Flag_T side);

void 
populate_CS_patch_SplitCS
  (
  Grid_T *const grid,
  const char *const obj0,/* NS, BH or etc. */
  const Flag_T dir0/* LEFT or RIGHT or CENTER or NONE */
  );

void 
populate_box_patch_SplitCS
  (
  Grid_T *const grid,
  const char *const obj0,/* filling_box,central_box. */
  const Flag_T dir0,/* direction */
  const char *const region/* covering region */
  );

int 
IsItCovering
  (
  const Patch_T *const patch,/* the patch */
  const char *const region/* BH/NS etc. see the list above */
  );
  
  

Patch_T **
collect_patches
  (
  Grid_T *const grid,/* the grid */
  const char *const region,/* see the list in IsItCovering function */
  Uint *const Np/* number of patches found */
  );

Grid_Char_T *init_grid_char(Grid_T *const new_grid);
void free_grid_char(Grid_Char_T *g);
Grid_Kind_T set_grid_kind(const char *const grid_kind);

void 
find_XYZ_and_patch_of_theta_phi_CS
 (
 double *const X/* found X,Y,Z Note: X[2] must be filled 
                // to determine the surface */,
 Patch_T **const ppatch,/* found patch */
 const double *const center/* center of S2 in general is not patch->c */,
 const double theta/* given theta */,
 const double phi/* given phi */,
 Patch_T **const patches,/* search among these patches */
 const Uint Np/* number of patches */
 );

Patch_T *x_in_which_patch(const double x[3],Patch_T **const patches,
                          const Uint Np);
Patch_T *X_in_which_patch(const double x[3],Patch_T **const patches,
                          const Uint Np);

void find_theta_phi_of_XYZ_CS(double *const theta,double *const phi,
                              const double *const X,const Flag_T side);                          

Patch_T **
collect_patches_regex
  (
  Grid_T *const grid,/* the grid */
  const char *const regex,/* regex */
  Uint *const Np/* number of patches found */
  );


void free_grid_params(const Grid_T *const grid);
Patch_T *x_in_which_patch_force(const double x[3],Patch_T **const patches,
                                const Uint Np,double *const X);
void carryover_interfaces(Grid_T *const new_grid,Grid_T *const old_grid);


#endif





