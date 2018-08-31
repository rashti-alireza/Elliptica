#include "core_lib.h"
#include "coordinates_lib.h"
#include "macros_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_analytic_lib.h"
#include "maths_calculus_lib.h"

#define DO 1
#define NOT_DO 0

/* types of derivatives; new one "must" be added to one before the last */
enum FUNC_E
{
  FUNC = 0,
  FUNC_x = 1,/* this number must to be changed */
  FUNC_y,
  FUNC_z,
  FUNC_xx,
  FUNC_yy,
  FUNC_zz,
  FUNC_xy,
  FUNC_xz,
  FUNC_yz,
  FUNC_xyz,
  N_FUNC
};

int derivative_tests(Grid_T *const grid);
int interpolation_tests(Grid_T *const grid);
Interpolation_T *init_interpolation(void);
double execute_interpolation(Interpolation_T *const interp_s);
void plan_interpolation(Interpolation_T *const interp_s);
static Flag_T read_F(sFunc_Patch2Pdouble_T **const F,sFunc_Patch2Pdouble_T **const DataBase_func,const enum FUNC_E fn);
static void enum2strcat(enum FUNC_E e,char *const f_derivative);
static void enum2str(enum FUNC_E e,char *const str);
static Flag_T compare_derivative(const char *const name,const double *const numc,const double *const anac,const Patch_T *const patch,const char *const path);
static double *make_random_number(const unsigned N,const double min,const double max);
static int interpolation_tests_x(Field_T *const field,const double *const x,const unsigned N);
static int interpolation_tests_y(Field_T *const field,const double *const y,const unsigned N);
static int interpolation_tests_z(Field_T *const field,const double *const z,const unsigned N);
static int interpolation_tests_xy(Field_T *const field,const double *const x,const double *const y,const unsigned Nx,const unsigned Ny);
static int interpolation_tests_xz(Field_T *const field,const double *const x,const double *const z,const unsigned Nx,const unsigned Nz);
static int interpolation_tests_yz(Field_T *const field,const double *const y,const double *const z,const unsigned Ny,const unsigned Nz);
static int interpolation_tests_xyz(Field_T *const field,const double *const x,const double *const y,const double *const z,const unsigned Nx,const unsigned Ny,const unsigned Nz);


