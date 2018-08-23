#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "maths_calculus_lib.h"
#include "coordinates_lib.h"
#include "maths_approximation_lib.h"
#include "memory_managing_lib.h"

#define OMP_PARALLEL_PATCH(x) _Pragma ( #x )

typedef int fSolve_T (Grid_T *const grid);

int solve_eqs(Grid_T *const grid);
static int parallel_patch_method (Grid_T *const grid);
static int bndry_outerB_ppm(Boundary_Condition_T *const bc);
static int bndry_copy_ppm(Boundary_Condition_T *const bc);
static int bndry_interpolate_ppm(Boundary_Condition_T *const bc);
static int b_in_ax_b_whole_ppm(Patch_T *const patch,const unsigned cn);
static int b_in_ax_b_bndry_ppm(Patch_T *const patch,const unsigned cn);
static double *normal_vec_curvilinear(Point_T *const point);
