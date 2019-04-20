#include "core_lib.h"
#include "error_handling_lib.h"
#include "memory_managing_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"

int make_patches(Grid_T *const grid);
static void fill_patches(Grid_T *const grid);
void fill_patches_Cartesian_grid(Grid_T *const grid);
void fill_patches_BNS_Projective_grid(Grid_T *const grid);
void fill_patches_BNS_CubedSpherical_grid(Grid_T *const grid);
void fill_patches_BNS_Spherical_grid(Grid_T *const grid);
int make_nodes(Grid_T *const grid);
int make_JacobianT(Grid_T *const grid);
void check_houseK(Patch_T *const patch);
void flush_houseK(Patch_T *const patch);
Patch_T make_temp_patch(const Patch_T *const patch);
void free_temp_patch(Patch_T *const patch);
