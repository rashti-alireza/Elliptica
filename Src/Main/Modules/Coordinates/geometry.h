#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"


/* type point */
enum Type
{
  INNER,
  EDGE
};

/* adjacent info */
typedef struct ADJACENT_T
{
  int p;// adjacent patch
  int f[TOT_FACE];// 0 if not located on an interface, positive otherwise
  int node;// node refers to index of adjacent point if any
  double N2[TOT_FACE][3];// normal of this adjPnt
  double N1dotN2[TOT_FACE];// dot product of normals(adjPnt,point)
  int FaceFlg;// equals 1 if found on an interface 0 otherwise
  int FitF;// equals face number if (N1dotN2 == -1); otherwise -1
}Adjacent_T;

/* points to be studied for realizing of geometry */
typedef struct POINTSET_T
{
  Point_T    *point;// the point under study
  Adjacent_T *adjPnt;// its adjacent points
  int NadjPnt;// number of adjacent point
  int FitN;// fittest number for adjPnt
}PointSet_T;

static void fill_basics(Patch_T *patch);
static void fill_N(Patch_T *patch);
static void fill_geometry(Grid_T *grid);
static void FindInnerB_Cartesian_coord(Patch_T *patch);
static void FindExterF_Cartesian_coord(Patch_T *patch);
static void init_Points(Interface_T *interface,PointSet_T ***innP,PointSet_T ***edgP);
static void set_min_max_sum(int *n,int f,int *im,int *iM,int *jm,int *jM,int *km,int *kM,int *sum);
static void free_PointSet(PointSet_T **pnt);
static void alloc_PointSet(int N,PointSet_T ***pnt);
static void realize_adj(PointSet_T **Pnt,enum Type type);
static void find_adjPnt(PointSet_T *Pnt,enum Type type);
static void find_adjNode(PointSet_T *pnt,const int N);
static void analyze_adjPnt(PointSet_T *Pnt,enum Type type);
static void add_adjPnt(PointSet_T *pnt,int *p, int np);
static void normal_vec_Cartesian_coord(Point_T *point);
static void tangent(Point_T *pnt,double *N);
static int NumPoint(Interface_T *interface,enum Type type);
static int L2(int *n,int f, int i, int j, int k);
static int  realize_neighbor(Patch_T *patch);
static int IsOverlap(PointSet_T *Pnt);
static int IsNormalFit(PointSet_T *Pnt);
static Point_T *get_p2(PointSet_T *Pnt);
void point_finder(Needle_T *needle);
double *normal_vec(Point_T *point);
void needle_ex(Needle_T *needle,Patch_T *patch);
void needle_in(Needle_T *needle,Patch_T *patch);
void flush_houseK(Patch_T *patch);
int find_node(double *x, Patch_T *patch);
