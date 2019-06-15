#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "memory_managing_lib.h"
#include "coordinates_lib.h"

enum enum_dA_da
{
  da_dx = 0,
  da_dy,
  da_dz,
  db_dx,
  db_dy,
  db_dz,
  dc_dx,
  dc_dy,
  dc_dz,
  dA_da_UNDEFINED
};

/* returning value */
struct Ret_S
{
  char s0[20],s1[20],s2[20];
};

void make_keyword_parameter(struct Ret_S *const ret,const char *const box,const char *const needle);
enum enum_dA_da get_dA_da(const Dd_T q2_e, const Dd_T q1_e);
double dq2_dq1(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
