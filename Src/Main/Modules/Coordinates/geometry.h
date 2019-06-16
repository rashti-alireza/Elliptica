#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"


/* type point */
enum Type
{
  INNER,
  EDGE
};

/* Face and normal properties for adjPnt */
struct Face_S
{
  unsigned on_f;/* if on this face 1; otherwise 0 */
  unsigned FitFlg;/* 1 if (N1dotN2 == -1); otherwise 0 */
  unsigned OrthFlg;/* 1 if (N1dotN2==0); othewise 0 */
  double N2[3];/*normal at a this face*/
  double N1dotN2;/* dot product of normals(adjPnt,point) */
};

/* adjacent info */
typedef struct ADJPOINT_T
{
  unsigned p;/* adjacent patch */
  unsigned FaceFlg;/* equals 1 if found on an interface; 0 otherwise */
  unsigned node;/* node refers to index of adjacent point if any */
  unsigned on_c;/* 1 if it is collocation, 0 otherwise */
  struct Face_S fs[TOT_FACE];
  unsigned InterpFace;/* face number for interpolation */
  unsigned CopyFace;/* face number for copy */
}AdjPoint_T;

/* points to be studied for realizing of geometry */
typedef struct POINTSET_T
{
  Point_T     *Pnt;/* the point under study */
  AdjPoint_T *adjPnt;/* its adjacent points */
  unsigned NadjPnt;/* number of adjacent point */
  unsigned idFit;/* id of adjPnt for fittest case */
  unsigned idOrth;/* id of adjPnt for Orthogonal cases */
  unsigned idInterp;/* id of adjPnt for interpolation case */
  unsigned overlap;/* 1 if overlaps, 0 otherwise */
  enum Type type;
}PointSet_T;

/* subface struct for setting some flags like df_dn and pairing purposes */
typedef struct SUBF_T
{
  SubFace_T *subf;/* subface */
  unsigned n;/* its number in a chain, ex: chain[n] = this struct */
  unsigned set  : 1;/* set 1 unset 0 */
  unsigned df_dn: 1;
  struct SUBF_T *next;/* the next subface in which this subface connected */
  struct SUBF_T *prev;/* the previous subface in which itconnected to this subface */
}Subf_T;

static void fill_basics(Patch_T *const patch);
static void fill_N(Patch_T *const patch);
static void fill_geometry(Grid_T * const grid);
static void FindInnerB_Cartesian_coord(Patch_T *const patch);
static void FindExterF_Cartesian_coord(Patch_T *const patch);
static void init_Points(const Interface_T *const interface,PointSet_T ***const innP,PointSet_T ***const edgP);
static void set_min_max_sum(const unsigned *const n,const unsigned f,unsigned *const im,unsigned *const iM,unsigned *const jm,unsigned *const jM,unsigned *const km,unsigned *const kM,unsigned *const sum);
static void free_PointSet(PointSet_T **pnt);
static void alloc_PointSet(const unsigned N,PointSet_T ***const pnt);
static int realize_adj(PointSet_T **const Pnt);
static void find_adjPnt(PointSet_T *const Pnt);
static void fill_adjPnt(PointSet_T *const pnt,const unsigned N);
static void analyze_adjPnt(PointSet_T *const Pnt);
static void add_adjPnt(PointSet_T *const pnt,const unsigned *const p, const unsigned np);
static void normal_vec_Cartesian_coord(Point_T *const point);
static void set_sameXYZ(Point_T *const p,const unsigned f);
void tangent(const Point_T *const pnt,double *const N);
static unsigned NumPoint(const Interface_T *const interface,const enum Type type);
static unsigned L2(const unsigned *const n,const unsigned f, const unsigned i, const unsigned j, const unsigned k);
static int realize_neighbor(Patch_T *const patch);
static int IsOverlap(PointSet_T *const Pnt);
static int IsNormalFit(PointSet_T *const Pnt);
static int IsOrthOutBndry(PointSet_T *const Pnt);
static int IsOutBndry(PointSet_T *const Pnt);
static int IsInnBndry(PointSet_T *const Pnt);
static int IsInterpolation(PointSet_T *const Pnt);
static int IsOnSubface(const Point_T *const pnt, const char *const lead);
static int ReachBnd(PointSet_T *const Pnt,const unsigned p,const unsigned f);
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
void flush_houseK(Patch_T *const patch);
void check_houseK(const Patch_T *const patch);
unsigned find_node(const double *const x, const Patch_T *const patch,Flag_T *const flg);
unsigned node_onFace(const double *const x, const unsigned f,const Patch_T *const patch);
int realize_geometry(Grid_T *const grid);
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
static void set_df_dn(Subf_T *const ring,const unsigned df_dn);

