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
  K_n2
};

typedef struct POINTSET_T
{
  Point_T *point;
  int *adjFace;// adjFace[0]...adjFace[nf-1] says 
               // the face number of adjacent faces which each 
               // of thess points can reach
  int *adjPatch;// adjPatch is like adjFace but for patch
  int nf;// number of face
  int np;// number of point
  int npatch;// number of patch
}PointSet_T;

static void fill_basics(Patch_T *patch);
static void fill_N(Patch_T *patch);
static void fill_geometry(Grid_T *grid);
static void FindInnerB_Cartesian_coord(Patch_T *patch);
static void FindExterF_Cartesian_coord(Patch_T *patch);
static void FindNeighbor(Patch_T *patch);
double *normal_vec(Point_T *point);
static void normal_vec_Cartesian_coord(Point_T *point);
