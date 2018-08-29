#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "maths_calculus_lib.h"
#include "coordinates_lib.h"
#include "maths_approximation_lib.h"
#include "memory_managing_lib.h"
#include "maths_solvers_lib.h"
#include "maths_general_lib.h"

#define PARALLEL_PATCH_METHOD_OpenMP(x) _Pragma ( #x )
static const double _Dirac_Delta_[2] = {0.0,1.0};/* dirac's delta */
#define Dirac_Delta(x,y) _Dirac_Delta_[(x)==(y)]

int parallel_patch_method (Grid_T *const grid);
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
static int solve_ax_b_ppm(Patch_T *const patch,const unsigned cn);
static int b_bndry_outerB_ppm(Boundary_Condition_T *const bc);
static int b_bndry_copy_ppm(Boundary_Condition_T *const bc);
static int b_bndry_interpolate_ppm(Boundary_Condition_T *const bc);
static int a_bndry_outerB_ppm(Boundary_Condition_T *const bc);
static int a_bndry_copy_ppm(Boundary_Condition_T *const bc);
static int a_bndry_interpolate_ppm(Boundary_Condition_T *const bc);
static int b_in_ax_b_whole_ppm(Patch_T *const patch,const unsigned cn);
static int b_in_ax_b_bndry_ppm(Patch_T *const patch,const unsigned cn);
static int a_in_ax_b_whole_ppm(Patch_T *const patch,const unsigned cn);
static int a_in_ax_b_bndry_ppm(Patch_T *const patch,const unsigned cn);
static double *normal_vec_curvilinear(Point_T *const point);
