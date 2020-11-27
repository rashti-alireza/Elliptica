#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"
#include "maths_equation_solvings_lib.h"
#include "manifold_header.h"

/* string length */
#define STR_LEN1 (99)

extern Grid_T **grids_global;

int make_patches(Grid_T *const grid);
void free_grid_db(void);
static void fill_patches(Grid_T *const grid);
void fill_patches_SNS_CubedSpherical_grid(Grid_T *const grid);
void fill_patches_Cartesian_grid(Grid_T *const grid);
void fill_patches_BNS_CubedSpherical_grid(Grid_T *const grid);
void fill_patches_BBN_CubedSpherical_grid(Grid_T *const grid);
void fill_patches_BNS_Spherical_grid(Grid_T *const grid);
void fill_patches_SNS_CubedSpherical_Box_grid(Grid_T *const grid);
void fill_patches_SBH_CubedSpherical_grid(Grid_T *const grid);
void alloc_patches_Cartesian_grid(Grid_T *const grid);
void alloc_patches_BNS_Spherical_grid(Grid_T *const grid);
void alloc_patches_BNS_CubedSpherical_grid(Grid_T *const grid);
void alloc_patches_BBN_CubedSpherical_grid(Grid_T *const grid);
void alloc_patches_SNS_CubedSpherical_Box_grid(Grid_T *const grid);
void alloc_patches_SNS_CubedSpherical_grid(Grid_T *const grid);
void alloc_patches_SBH_CubedSpherical_grid(Grid_T *const grid);
void fill_patches_Split_CubedSpherical_grid(Grid_T *const grid);
int make_nodes(Grid_T *const grid);
int make_JacobianT(Grid_T *const grid);
Patch_T make_temp_patch(const Patch_T *const patch);
void free_temp_patch(Patch_T *const patch);
static char *coord_sys_str(const Patch_T *const patch,char *const str);
static char *bases_str(const Patch_T *const patch,char *const str);
static char *collocation_str(const Patch_T *const patch,char *const str);
void *alloc_grid(void);
void alloc_patches(Grid_T *const grid);
void free_grid(Grid_T *grid);
void free_patch(Patch_T *patch);
void free_grid_params(const Grid_T *const grid);

