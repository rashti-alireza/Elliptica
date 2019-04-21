#include "coordinate_shared_lib.h"

void make_nodes_Spherical_coord(Patch_T *const patch);
void fill_patches_BNS_Spherical_grid(Grid_T *const grid);
void populate_left_outermost(Grid_T *const grid,const unsigned pn,const unsigned outermost_n);
void populate_right_outermost(Grid_T *const grid,const unsigned pn,const unsigned outermost_n);
void alloc_patches_BNS_Spherical_grid(Grid_T *const grid);
static void populate_left_NS_sphere(Grid_T *const grid,const unsigned pn);
static void populate_left_NS_surrounding_sphere(Grid_T *const grid,const unsigned pn);
static void populate_right_NS_sphere(Grid_T *const grid,const unsigned pn);
static void populate_right_NS_surrounding_sphere(Grid_T *const grid,const unsigned pn);
