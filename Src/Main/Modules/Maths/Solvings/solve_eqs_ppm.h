#include "core_lib.h"
#include "utilities_lib.h"
#include "error_handling_lib.h"
#include "maths_calculus_lib.h"
#include "coordinates_lib.h"
#include "maths_approximation_lib.h"
#include "memory_managing_lib.h"
#include "maths_solvers_lib.h"
#include "maths_general_lib.h"
#include "maths_linear_algebra_lib.h"

#define PARALLEL_PATCH_METHOD_OpenMP(x) _Pragma ( #x )
static const double _Dirac_Delta_[2] = {0.0,1.0};/* dirac's delta */
#define Dirac_Delta(x,y) _Dirac_Delta_[(x)==(y)]

int parallel_patch_method (Grid_T *const grid);
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
static void initialize_ppm(Grid_T *const grid);
static void copy_initial_values_ppm(Grid_T *const grid);
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
static void fill_interpolation_flags(Interpolation_T *const it,const SubFace_T *const sf);
static unsigned const_index_on_AdjFace(const SubFace_T *const sf);
fJs_T *get_j_reader(const Matrix_T *const m);
static void update_fields_ppm(Grid_T *const grid);
static Flag_T check_residual(const Grid_T *const grid,const double res_input);
