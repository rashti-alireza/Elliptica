#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"

/* face number */
enum Face
{
  I_0 = 0,
  I_n0,
  J_0,
  J_n1,
  K_0,
  K_n2,
  TOT_FACE
};

/* type point */
enum Type
{
  INNER,
  EDGE
};

typedef struct ADJACENT_T
{
  int p;// adjacent patch
  int f;// adjacent face
}Adjacent_T;

/* points to be studied for realizing of geometry */
typedef struct POINTSET_T
{
  Point_T *point;
  Adjacent_T *adj;
  int Np;// number of point
  int Nadj;// number of adjacent
  enum Type type;// inner or edge point type
}PointSet_T;

static void fill_basics(Patch_T *patch);
static void fill_N(Patch_T *patch);
static void fill_geometry(Grid_T *grid);
static void FindInnerB_Cartesian_coord(Patch_T *patch);
static void FindExterF_Cartesian_coord(Patch_T *patch);
static void init_Points(Interface_T *interface,PointSet_T ***innP,PointSet_T ***edgP);
static int L2(int *n,int f, int i, int j, int k);
static void set_min_max_sum(int *n,int f,int *im,int *iM,int *jm,int *jM,int *km,int *kM,int *sum);
static void free_PointSet(PointSet_T **pnt);
static void *alloc_PointSet(int N);
static int NumPoint(Interface_T *interface,enum Type type);
static int  RealizeNeighbor(Patch_T *patch);
double *normal_vec(Point_T *point);
static void normal_vec_Cartesian_coord(Point_T *point);
