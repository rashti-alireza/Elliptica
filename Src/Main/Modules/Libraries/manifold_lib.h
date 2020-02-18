#ifndef manifold_lib
#define manifold_lib



/* interrelation structures */
struct FIELD_T;
struct SOLVING_MAN_T;


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
  _a_/* for Curvilinear 0-coord */,
  _b_/* for Curvilinear 1-coord */,
  _c_/* for Curvilinear 2-coord */,
  N_DD_T/* number of total DD_T*/,
  UNDEFINED_DIR/* for referring to undefined Dd_T */,
  _3_ = 3/* number three for tracking dimensions for different quantities */
}Dd_T;

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
typedef struct INTERFACE_T
{
  struct PATCH_T *patch;/* refers to its patch */
  unsigned np;/* number of points in Point_T */
  unsigned fn;/* its interface number */
  unsigned ns;/* number of subfaces */
  Point_T **point;/* points on the interface */
  SubFace_T **subface;/* subset of points on this interface with same flags */
}Interface_T;

struct Collocation_s
{
  double min;
  double max;
  unsigned n;
  double stp;
  double a;
  double b;
  Collocation_T c;
};

/* Jacobian transformation between different coords  */
typedef struct JACOBIAN_TRANS_T
{
  double (*j)(struct PATCH_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);/* function for transformation */
  double *dX_dx[3][3];/* saving some transformation to save time for dX[0..2]/dx[0..2] */
  double *dx_dX[3][3];/* saving some transformation to save time dx[0..2]/dX[0..2] */
}JacobianTrans_T;

/* patch */
typedef struct PATCH_T
{
  struct GRID_T *grid;/* refers to its grid */
  char *name;/* patch name */
  Coord_T coordsys;/* coord sys used in this patch */
  struct
  {
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
  struct FIELD_T **pool;/* pool of fields, 
                 // notation: pool[Ind("Phi_f")] refers 
                 // to field phi. note: Ind is macro.
                 // one can access to values of field on this 
                 // patch like pool[Ind("phi_f")]->v */
  struct SOLVING_MAN_T *solving_man;/* solving management */
  unsigned innerB:1;/* if this patch has inner boundary 1 otherwise 0 */
  unsigned outerB:1;/* if this patch has outer boundary 1 otherwise 0 */
  unsigned is_a_closed: 1;/* if coordinate a is periodic or closed 1, otherwise 0 */
  unsigned is_b_closed: 1;/* if coordinate b is periodic or closed 1, otherwise 0 */
  unsigned is_c_closed: 1;/* if coordinate c is periodic or closed 1, otherwise 0 */
}Patch_T;

/* grid */
typedef struct GRID_T
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
  struct FIELD_T **field;
  unsigned nf;/* number of fields */
}Jacobian_Eq_T;

/* needle which should be found.
// it is a general purpose structure could be used by different functions
// like,point_finder function */
typedef struct NEEDLE_T
{
  const double *x;/* Cartesian coords */
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


double point_value(const unsigned i, const struct Collocation_s *const coll_s);
void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const colloc,const unsigned dir);
int make_patches(Grid_T *const grid);
int realize_geometry(Grid_T *const grid);
int X_of_x(double *const X,const double *const x,const Patch_T *const patch);
int x_of_X(double *const x,const double *const X,const Patch_T *const patch);
double x_coord(const unsigned i,const Patch_T *const patch);
double y_coord(const unsigned i,const Patch_T *const patch);
double z_coord(const unsigned i,const Patch_T *const patch);
double X_coord(const unsigned i,const Patch_T *const patch);
double Y_coord(const unsigned i,const Patch_T *const patch);
double Z_coord(const unsigned i,const Patch_T *const patch);
double dq2_dq1(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_Cartesian_patch(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double *normal_vec(Point_T *const point);
double General2ChebyshevExtrema(const double X,const unsigned dir,const Patch_T *const patch);
void grid_characteristics_example(Grid_T *const grid);
double R_interpolation_CS(struct FIELD_T *const R,const double *const X);
void SignAndIndex_permutation_CubedSphere(const Flag_T side,unsigned *const a,unsigned *const b,unsigned *const c,double *const s);
void test_CubedSpherical_Coordinates(Grid_T *const grid);
void test_dq_dN(Grid_T *const grid);
void needle_ex(Needle_T *const needle,const Patch_T *const patch);
void needle_in(Needle_T *const needle,const Patch_T *const patch);
void point_finder(Needle_T *const needle);
void populate_right_BH(Grid_T *const grid,const unsigned pn);
void populate_right_BH_central_box(Grid_T *const grid,const unsigned pn);
int make_nodes(Grid_T *const grid);
void alloc_nodes(Grid_T *const grid);
void alloc_interface(Patch_T *const patch);
void *alloc_point(const unsigned s);
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




#endif





