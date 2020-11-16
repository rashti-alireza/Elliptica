#include "coordinate_shared_lib.h"
#include "maths_calculus_lib.h"
#include "maths_spectral_methods_lib.h"
#include "fields_lib.h"

#define Power3(a) ((a)*Pow2(a))

void populate_central_NS_central_box(Grid_T *const grid,const unsigned pn);
void populate_right_box_sns(Grid_T *const grid,const unsigned pn); 
void fill_patches_SBH_CubedSpherical_grid(Grid_T *const grid);
void fill_patches_SNS_CubedSpherical_Box_grid(Grid_T *const grid);
void fill_patches_SNS_CubedSpherical_grid(Grid_T *const grid);
void fill_patches_BNS_CubedSpherical_grid(Grid_T *const grid);
void fill_patches_BBN_CubedSpherical_grid(Grid_T *const grid);
void SignAndIndex_permutation_CubedSphere(const Flag_T side,unsigned *const a,unsigned *const b,unsigned *const c,double *const s);
void make_nodes_CubedSpherical_coord(Patch_T *const patch);
void populate_left_NS_central_box(Grid_T *const grid,const unsigned pn);
void populate_right_NS_central_box(Grid_T *const grid,const unsigned pn);
void populate_filling_box_CubedSpherical(Grid_T *const grid,const unsigned pn,const Flag_T side);
void alloc_patches_SBH_CubedSpherical_grid(Grid_T *const grid);
void alloc_patches_BNS_CubedSpherical_grid(Grid_T *const grid);
void alloc_patches_BBN_CubedSpherical_grid(Grid_T *const grid);
void alloc_patches_SNS_CubedSpherical_Box_grid(Grid_T *const grid);
void make_JacobianT_CubedSpherical_coord(Patch_T *const patch);
double JT_NS_T_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_NS_T_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_NS_T_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_NS_T_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_NS_T_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_NS_T_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_SR_T_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_SR_T_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_SR_T_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_SR_T_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_SR_T_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_SR_T_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T1_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T1_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T1_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T1_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T1_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T1_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T2_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T2_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T2_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T2_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T2_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T2_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
static void populate_left_NS(Grid_T *const grid,const unsigned pn);
static void populate_left_NS_surrounding(Grid_T *const grid,const unsigned pn);
static void populate_right_NS(Grid_T *const grid,const unsigned pn);
static void populate_right_NS_surrounding(Grid_T *const grid,const unsigned pn);
static void populate_right_BH_surrounding(Grid_T *const grid,const unsigned pn);
static void populate_outermost(Grid_T *const grid,const unsigned pn,const unsigned i);
static void populate_central_NS(Grid_T *const grid,const unsigned pn);
static void populate_central_NS_surrounding(Grid_T *const grid,const unsigned pn);
static void populate_central_BH_surrounding(Grid_T *const grid,const unsigned pn);
void populate_right_BH(Grid_T *const grid,const unsigned pn);
void populate_right_BH_central_box(Grid_T *const grid,const unsigned pn);
void populate_filling_box(Grid_T *const grid,const unsigned pn);
static void R1_derivative(Patch_T *const patch);
static void R2_derivative(Patch_T *const patch);  
double R_interpolation_CS(Field_T *const R,const double *const X);
int x_of_X(double *const x,const double *const X,const Patch_T *const patch);
void test_CubedSpherical_Coordinates(Grid_T *const grid);
void alloc_patches_SNS_CubedSpherical_grid(Grid_T *const grid);
void alloc_patches_Split_CubedSpherical_grid(Grid_T *const grid);
void set_params_split_CS(Grid_Char_T *const grid_char);
void fill_patches_Split_CubedSpherical_NS_BH_grid(Grid_T *const grid);
void set_object_name_split_CS(char *const obj,const char *const type);
static void R12_derivatives_SCS(Patch_T *const patch);

double JT_OJ_T_SCS(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_OT_T_SCS(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);

void 
populate_CS_patch_SplitCS
  (
  Grid_T *const grid,
  const char *const obj0,/* NS, BH or etc. */
  const Flag_T dir0,/* LEFT or RIGHT or NONE */
  unsigned *const pn/* starting patch number,is increased for each add */
  );
