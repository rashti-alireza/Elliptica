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

static void fill_basics(Patch_T *const patch);
static void fill_N(Patch_T *const patch);
static void fill_geometry(Grid_T *grid);
static void FindInnerB_Cartesian_coord(Patch_T *patch);
static void FindExterF_Cartesian_coord(Patch_T *patch);
//static void FindNeighbor_Cartesian_grid(Patch_T *patch);
double *normal_vec(Point_T *point);
static void normal_vec_Cartesian_coord(Point_T *point);
