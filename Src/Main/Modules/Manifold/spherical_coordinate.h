#include "manifold_header.h"

void make_nodes_Spherical_coord(Patch_T *const patch);
void fill_patches_BNS_Spherical_grid(Grid_T *const grid);
void alloc_patches_BNS_Spherical_grid(Grid_T *const grid);
static void populate_left_NS_sphere(Grid_T *const grid,const Uint pn);
static void populate_left_NS_around_sphere(Grid_T *const grid,const Uint pn);
static void populate_right_NS_sphere(Grid_T *const grid,const Uint pn);
static void populate_right_NS_around_sphere(Grid_T *const grid,const Uint pn);
