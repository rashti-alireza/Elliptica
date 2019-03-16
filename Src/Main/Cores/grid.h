#include "core_lib.h"
#include "error_handling_lib.h"
#include "memory_managing_lib.h"
#include "utilities_lib.h"

/* returning value */
struct Ret_S
{
  char s0[20],s1[20],s2[20];
};

int make_patches(Grid_T *const grid);
static void fill_patches(Grid_T *const grid);
static void fill_patches_Cartesian_grid(Grid_T *const grid);
static void fill_patches_BNS_Projective_grid(Grid_T *const grid);
static void populate_left_NS_central_box(Grid_T *const grid,const unsigned pn);
static void populate_left_NS_hemisphere_up(Grid_T *const grid,const unsigned pn);
static void populate_left_NS_hemisphere_down(Grid_T *const grid,const unsigned pn);
static void populate_left_NS_surrounding_up(Grid_T *const grid,const unsigned pn);
static void populate_left_NS_surrounding_down(Grid_T *const grid,const unsigned pn);
static void populate_right_NS_central_box(Grid_T *const grid,const unsigned pn);
static void populate_right_NS_hemisphere_up(Grid_T *const grid,const unsigned pn);
static void populate_right_NS_hemisphere_down(Grid_T *const grid,const unsigned pn);
static void populate_right_NS_surrounding_up(Grid_T *const grid,const unsigned pn);
static void populate_right_NS_surrounding_down(Grid_T *const grid,const unsigned pn);
static void populate_left_outermost(Grid_T *const grid,const unsigned pn,const unsigned outermost_n);
static void populate_right_outermost(Grid_T *const grid,const unsigned pn,const unsigned outermost_n);
static void make_keyword_parameter(struct Ret_S *const ret,const char *const box,const char *const needle);
int make_nodes(Grid_T *const grid);
int make_JacobianT(Grid_T *const grid);
void check_houseK(Patch_T *const patch);
void flush_houseK(Patch_T *const patch);
Patch_T make_temp_patch(const Patch_T *const patch);
void free_temp_patch(Patch_T *const patch);
