#include "manifold_header.h"

void make_nodes_Cartesian_coord(Patch_T *const patch);
void make_JacobianT_Cartesian_coord(Patch_T *const patch);
void populate_left_NS_central_box(Grid_T *const grid,const Uint pn);
void populate_right_NS_central_box(Grid_T *const grid,const Uint pn);
void populate_filling_box_CubedSpherical(Grid_T *const grid,const Uint pn,const Flag_T side);
void alloc_patches_Cartesian_grid(Grid_T *const grid);
double JT_Cartesian_patch(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p);
void fill_patches_Cartesian_grid(Grid_T *const grid);
void populate_right_box_sns(Grid_T *const grid,const Uint pn);
void populate_central_NS_central_box(Grid_T *const grid,const Uint pn);
void populate_right_BH_central_box(Grid_T *const grid,const Uint pn);


void 
populate_box_patch_SplitCS
  (
  Grid_T *const grid,
  const char *const obj0,/* filling_box,central_box. */
  const Flag_T dir0,/* direction */
  const char *const region/* covering region */
  );
  
  
  
