#include "core_lib.h"
#include "coordinates_lib.h"
#include "macros_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_analytic_lib.h"

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

int DerivativeTest(Grid_T *const grid);
static void ChebyshevFirstKindBasis_DerivativeTest(const Patch_T *const patch);
static Flag_T read_F(sFunc_Grid2Pdouble_T **const F,sFunc_Grid2Pdouble_T **const func,const enum FUNC_E fn);
static void enum2strcat(enum FUNC_E e,char *const f_derivative);
