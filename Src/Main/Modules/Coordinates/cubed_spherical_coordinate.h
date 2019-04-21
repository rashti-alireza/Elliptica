#include "coordinate_shared_lib.h"

void fill_patches_BNS_CubedSpherical_grid(Grid_T *const grid);
void SignAndIndex_permutation_CubedSphere(const Flag_T side,unsigned *const a,unsigned *const b,unsigned *const c,double *const s);
void make_nodes_CubedSpherical_coord(Patch_T *const patch);
void populate_left_NS_central_box(Grid_T *const grid,const unsigned pn);
void populate_right_NS_central_box(Grid_T *const grid,const unsigned pn);
void populate_filling_box_CubedSpherical(Grid_T *const grid,const unsigned pn,const Flag_T side);
void alloc_patches_BNS_CubedSpherical_grid(Grid_T *const grid);
static void populate_left_NS(Grid_T *const grid,const unsigned pn);
static void populate_left_NS_surrounding(Grid_T *const grid,const unsigned pn);
static void populate_right_NS(Grid_T *const grid,const unsigned pn);
static void populate_right_NS_surrounding(Grid_T *const grid,const unsigned pn);
static void populate_outermost(Grid_T *const grid,const unsigned pn,const unsigned i);
void populate_filling_box(Grid_T *const grid,const unsigned pn);
  