#ifndef manifold_header_LIB_H
#define manifold_header_LIB_H

#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_special_functions_lib.h"
#include "maths_calculus_lib.h"

/* internal parameter prefix 
// PLEASE keep it capitalized. */
#define P_ "MNFLD_"

/* array size */
#define STR_SIZE1 (100)
#define STR_SIZE2 (200)
#define STR_SIZE3 (999)

/* some conventions: */
/* split cubed spherical suffix */
#define SCS_suffix "_X%uY%uZ%u"

/* surface function stem */
#define SigmaU "sigmaU"
#define SigmaD "sigmaD"

/* par format for surface function.
// NOTE: d0,d1,d2,grid,dir,obj and side are required.
// ex: grid1_sigmaU_left_NS_backX0Y2Z3 */
#define SCS_par_sigma(par,sigma) \
  sprintf(par,PATCH_NAME_PRT_P_"%s_%s_%s_%s"SCS_suffix,\
  grid->gn,sigma,dir,obj,StrSide[side],d0,d1,d2);

/* par format for center of patch with respect to the 
// reference Cartesian coords. 
// NOTE: grid,dir,obj and side are required.
// ex: grid1_left_BH_front_center_a */
#define SCS_par_CS_center(par,axis) \
  sprintf(par,PATCH_NAME_PRT_P_"%s_%s_center_%s_%s",\
  grid->gn,dir,obj,StrSide[side],axis);

/* par format for center of patch with respect to the 
// reference Cartesian coords. 
// NOTE: grid,dir,obj,d0,d1,d2 and side are required.
// ex: grid1_left_central_box_left_center_aX0Y0Z1 */
#define SCS_par_box_center(par,axis) \
  sprintf(par,PATCH_NAME_PRT_P_"%s_%s_center_%s_%s"SCS_suffix,\
  grid->gn,dir,obj,StrSide[side],axis,d0,d1,d2);

/* par format for box lengths.
// NOTE: d0,d1,d2,grid,dir and obj are required.
// ex: grid1_central_box_len_a_leftX0Y2Z3 */
#define SCS_par_box_length(par,axis) \
  sprintf(par,PATCH_NAME_PRT_P_"%s_len_%s_%s"SCS_suffix,\
  grid->gn,obj,axis,dir,d0,d1,d2);

/* par format for xc lengths.
// NOTE: d0,d1,d2,grid,dir,side and obj are required.
// ex: grid1_right_BH_around_up_xc2_X0Y2Z3 */
#define SCS_par_xc_length(par,xc) \
  sprintf(par,PATCH_NAME_PRT_P_"%s_%s_%s_%s"SCS_suffix,\
  grid->gn,dir,obj,StrSide[side],xc,d0,d1,d2);

/* par format for patch->min[index]
// NOTE: d0,d1,d2,grid,dir,obj and side are required.
// ex: grid1_left_NS_up_min2X0Y2Z3 */
#define SCS_par_min(par,index) \
  sprintf(par,PATCH_NAME_PRT_P_"%s_%s_%s_min%u"SCS_suffix,\
  grid->gn,dir,obj,StrSide[side],index,d0,d1,d2);

/* par format for patch->max[index].
// NOTE: d0,d1,d2,grid,dir,obj and side are required.
// ex: grid1_left_NS_up_max2X0Y2Z3 */
#define SCS_par_max(par,index) \
  sprintf(par,PATCH_NAME_PRT_P_"%s_%s_%s_max%u"SCS_suffix,\
  grid->gn,dir,obj,StrSide[side],index,d0,d1,d2);


/* par format for patch->name.
// NOTE: d0,d1,d2,grid,dir,obj and side are required.
// ex: grid1_left_NS_up_X0Y2Z3 */
#define SCS_par_name(par) \
  sprintf(par,PATCH_NAME_PRT_P_"%s_%s_%s"SCS_suffix,\
  grid->gn,dir,obj,StrSide[side],d0,d1,d2);


/* sides, NOTE: the order is important, it MUST be like FLAG_T */
static const char *const StrSide[] = 
  {"up","down","left",
  "right","back","front",
  "center",
  0};

/* allowed object: */
static const char *const SCS_ObjType[] = 
 {"NS","BH",
  "NS_around","BH_around",
  "NS1","BH1",
  "NS1_around","BH1_around",
  "NS2","BH2",
  "NS2_around","BH2_around",
  "central_box","filling_box",
  "outermost",0
 };

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
double dq2_dq1(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p);
void alloc_patches_Split_CubedSpherical_grid(Grid_T *const grid);
void set_object_name_split_CS(char *const obj,const char *const type);
Grid_Kind_T set_grid_kind(const char *const grid_kind);

#endif




