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

#define DDM_SCHUR_COMPLEMENT_OpenMP(x) _Pragma ( #x )


int ddm_schur_complement(Grid_T *const grid);
static void preparing_ingredients(Grid_T *const grid);
static void make_map_and_inv(Patch_T *const patch);
static Flag_T check_residual(const Grid_T *const grid,const double res_input);
static void make_f(Patch_T *const patch);
static char **read_fields_in_order(unsigned *const nf);
static int solve_field(Grid_T *const grid);
static void set_cf(Grid_T *const grid,const char *const field_name);
static void f_in_equation_part(Patch_T *const patch);
static void f_in_outerboundary_part(Patch_T *const patch);
static void make_g_partial(Patch_T *const patch);



