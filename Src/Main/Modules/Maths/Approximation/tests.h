#include "core_lib.h"
#include "coordinates_lib.h"
#include "macros_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_analytic_lib.h"

enum FUNC_E
{
  FUNC = 0,
  FUNC_x = 1,
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

int DerivativeTest(const Grid_T *const grid)
static void ChebyshevFirstKindBasis_DerivativeTest(const Patch_T *const patch);
