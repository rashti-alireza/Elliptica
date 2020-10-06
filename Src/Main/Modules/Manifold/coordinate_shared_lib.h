#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_special_functions_lib.h"
#include "maths_calculus_lib.h"

/* array size */
#define STR_SIZE1 (100)
#define STR_SIZE2 (200)
#define STR_SIZE3 (999)

/* some conventions: */
/* split cubed spherical suffix */
static const char *const scs_suffix = "X%uY%uZ%u";

/* prototype for surface function, ex: grid1_left_NS_sigmaU_back023 */
static const char *const scs_prototype_sigma = "grid%u_%s_%s_%s_%s%u%u%u";

/* surface function sprintf NOTE: d0,d1,d2 is required */
#define scs_wpar_sigma(Sstr,Ugridnumber,Sdirection,Sobject,Fside,Ssigma) \
  assert(sprintf(Sstr,scs_prototype_sigma,(Ugridnumber),Ssigma,\
         Sdirection,Sobject,StrSide[Fside],d0,d1,d2));

/* surface function stem */
#define SigmaU "sigmaU"
#define SigmaD "sigmaD"

/* center of patch with respect to the reference Cartesian coords. */
#define scs_wpar_center(Sstr,Ugridnumber,Sdirection,Sobject,Svar) \
  assert(sprintf(Sstr,"grid%u_%s_%s_center_%s",(Ugridnumber),\
         Sdirection,Sobject,Svar));

/* box lengths */
#define scs_wpar_box(Spar,Ugridnumber,Sdirection,Sobject,Svar) \
  assert(sprintf(Spar,"grid%u_%s_%s_%s%u%u%u",(Ugridnumber),\
         Sdirection,Sobject,Svar,d0,d1,d2));

/* sides, NOTE: the order is important, it MUST be like FLAG_T */
static const char *const StrSide[] = 
  {"up","down","left",
  "right","back","front",
  0};

/* enum for jacobian */
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
void alloc_patches_BBN_Split_CubedSpherical_grid(Grid_T *const grid);



