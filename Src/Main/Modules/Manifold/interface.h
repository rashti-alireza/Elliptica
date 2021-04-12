#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"

/* handy error sprintf, note: patch and f are required  */
#define err_spr_adj(s,x,adjpatch,adjface)  \
sprintf(s,\
  "At '%s' on %s for \n"\
  " point (%g,%g,%g) no adjacent patch found!\n"\
  " Possible adjacent patch is:\n '%s' on %s.\n"\
  " Please consider to increase the resolution.",\
  patch->name,FaceName[f],(x)[0],(x)[1],(x)[2],\
  adjpatch->name,FaceName[adjface]);

#define err_spr(s,x)  \
sprintf(s,\
  "At '%s' on %s for \n"\
  " point (%g,%g,%g) no adjacent patch found!\n"\
  " Please consider to increase the resolution.",\
  patch->name,FaceName[f],(x)[0],(x)[1],(x)[2]);


/* handy error sprintf for warning about normal match  */
#define warn_normal_spr(s,patch,adj_patch,adj_face)  \
   sprintf(s,"~> Warning!\nfor '%s' on %s and \nadjacent '%s' on %s:\n"\
   "normal vectors are not align: angle = %g degree(s).\n"\
   "Please consider to increase the resolution.",\
   patch->name,FaceName[f],(adj_patch ? adj_patch->name : "?"),\
   (adj_face < NFaces ? FaceName[adj_face]: "?"), \
    acos(N1dotN2)*180/M_PI);

/* type point */
enum Type
{
  INNER,
  EDGE
};

/* Face and normal properties for adjPnt */
struct Face_S
{
  Uint on_f;/* if on this face 1; otherwise 0 */
  Uint FitFlg;/* 1 if (N1dotN2 == -1); otherwise 0 */
  Uint OrthFlg;/* 1 if (N1dotN2==0); othewise 0 */
  double N2[3];/*normal at a this face*/
  double N1dotN2;/* dot product of normals(adjPnt,point) */
};

/* adjacent info */
typedef struct ADJPOINT_T
{
  Uint p;/* adjacent patch */
  Uint FaceFlg;/* equals 1 if found on an interface; 0 otherwise */
  Uint node;/* node refers to index of adjacent point if any */
  Uint on_c;/* 1 if it is collocation, 0 otherwise */
  struct Face_S fs[TOT_FACE];
  Uint InterpFace;/* face number for interpolation */
  Uint CopyFace;/* face number for copy */
}AdjPoint_T;

/* points to be studied for realizing of geometry */
typedef struct POINTSET_T
{
  Point_T     *Pnt;/* the point under study */
  AdjPoint_T *adjPnt;/* its adjacent points */
  Uint NadjPnt;/* number of adjacent point */
  Uint idFit;/* id of adjPnt for fittest case */
  Uint idOrth;/* id of adjPnt for Orthogonal cases */
  Uint idInterp;/* id of adjPnt for interpolation case */
  Uint overlap;/* 1 if overlaps, 0 otherwise */
  enum Type type;
}PointSet_T;

/* subface struct for setting some flags like df_dn and pairing purposes */
typedef struct SUBF_T
{
  SubFace_T *subf;/* subface */
  Uint n;/* its number in a chain, ex: chain[n] = this struct */
  Uint set  : 1;/* set 1 unset 0 */
  Uint df_dn: 1;
  struct SUBF_T *next;/* the next subface in which this subface connected */
  struct SUBF_T *prev;/* the previous subface in which itconnected to this subface */
}Subf_T;

void carryover_interfaces(Grid_T *const new_grid,Grid_T *const old_grid);
static void fill_basics(Patch_T *const patch);
static void fill_geometry(Grid_T * const grid,Uint **const point_flag);
static void FindInnerB_Cartesian_coord(Patch_T *const patch);
static void FindExterF_Cartesian_coord(Patch_T *const patch);
static void init_Points(const Interface_T *const interface,PointSet_T ***const innP,PointSet_T ***const edgP);
static void set_min_max_sum(const Uint *const n,const Uint f,Uint *const im,Uint *const iM,Uint *const jm,Uint *const jM,Uint *const km,Uint *const kM,Uint *const sum);
static void free_PointSet(PointSet_T **pnt);
static void alloc_PointSet(const Uint N,PointSet_T ***const pnt);
static int realize_adj(PointSet_T **const Pnt,Uint **const point_flag);
static void find_adjPnt(PointSet_T *const Pnt);
static void fill_adjPnt(PointSet_T *const pnt,const Uint N);
static void analyze_adjPnt(PointSet_T *const Pnt,Uint **const point_flag);
static void add_adjPnt(PointSet_T *const pnt,const Uint *const p, const Uint np);
static void normal_vec_Cartesian_coord(Point_T *const point);
static void set_sameXYZ(Point_T *const p,const Uint f);
void tangent(const Point_T *const pnt,double *const N);
static Uint NumPoint(const Interface_T *const interface,const enum Type type);
static Uint L2(const Uint *const n,const Uint f, const Uint i, const Uint j, const Uint k);
static int realize_neighbor(Patch_T *const patch,Uint **const point_flag);
static int IsOverlap(PointSet_T *const Pnt);
static int IsNormalFit(PointSet_T *const Pnt);
static int IsOrthOutBndry(PointSet_T *const Pnt);
static int IsOrthInnBndry(PointSet_T *const Pnt);
static int IsOutBndry(PointSet_T *const Pnt);
static int IsInnBndry(PointSet_T *const Pnt);
static int IsInterpolation(PointSet_T *const Pnt);
static int IsOnSubface(const Point_T *const pnt, const char *const lead);
static int ReachBnd(PointSet_T *const Pnt,const Uint p,const Uint f);
static void test_subfaces(const Grid_T *const grid);
static char *inspect_flags(const Point_T *const pnt);
static void add_to_subface(const Point_T *const pnt,const char *const lead);
static void add_point(SubFace_T *const subface,const Point_T *const pnt);
static Point_T *get_p2(const PointSet_T *const Pnt);
static SubFace_T *find_subface(const SubFace_T *const sub);
void point_finder(Needle_T *const needle);
double *normal_vec(Point_T *const point);
void needle_ex(Needle_T *const needle,const Patch_T *const patch);
void needle_in(Needle_T *const needle,const Patch_T *const patch);
Uint find_node(const double *const x, const Patch_T *const patch,Flag_T *const flg);
Uint node_onFace(const double *const x, const Uint f,const Patch_T *const patch);
int realize_interfaces(Grid_T *const grid);
void make_normal_outward(Point_T *const point);
static void misc(Grid_T *const grid);
static int IsMatchedOtherInnerSubface(PointSet_T *const Pnt);
static void normal_vec_CS_coord(Point_T *const point);
static void FindInnerB_CS_coord(Patch_T *const patch);
static void FindExterF_CS_coord(Patch_T *const patch);
static void set_one_Dirichlet_BC(Interface_T **const face);
static void set_df_dn_and_pair(Grid_T *const grid);
static Subf_T *add_ring(Subf_T ***chain);
static Subf_T **compose_the_chain(SubFace_T *const subf1);
static void set_df_dn(Subf_T *const ring,const Uint df_dn);
void alloc_nodes(Grid_T *const grid);
void alloc_interface(Patch_T *const patch);
void *alloc_point(const Uint s);
void free_points(Grid_T *const grid);
void free_patch_interface(Patch_T *const patch);
static void ri_split_cubed_spherical(Grid_T *const grid);
static void ri_general_method(Grid_T *const grid);
static void check_houseK(Patch_T *const patch);
static void flush_houseK(Patch_T *const patch);
static void fill_N(Patch_T *const patch);
static void add_to_subface_scs(Point_T *const pnt);
static void add_point_scs(SubFace_T *const subface,const Point_T *const pnt);
static void set_subfaces_scs(Grid_T *const grid,Patch_T *const patch);

static void 
find_adjacent_scs
  (
  Grid_T *const grid,
  Patch_T *const patch,
  Uint *const point_flag
  );

static void pair_subfaces_and_set_bc(Grid_T *const grid);
static Uint counter_n_adjacent_faces(const Interface_T *const face);

