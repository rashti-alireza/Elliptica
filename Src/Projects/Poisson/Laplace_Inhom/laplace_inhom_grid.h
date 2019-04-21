#include "core_lib.h"
#include "memory_managing_lib.h"
#include "coordinates_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"

struct R_max
{
    double R_max_l;
    double R_max_r;
};

Grid_T *Laplace_Inhom_make_grid(void);
static void grid_characteristics_Laplace_Inhome (Grid_T *const grid);
static void characteristics_Cartesian_grid(Grid_T *const grid);
static void characteristics_BNS_Projective_grid(Grid_T *const grid);
static void NS_radii_BNS_Projective_grid(Grid_T *const grid,void *vp);
static void characteristics_BNS_Spherical_grid(Grid_T *const grid);
static void NS_radii_BNS_Spherical_grid(Grid_T *const grid,void *vp);
static void characteristics_BNS_CubedSpherical_grid(Grid_T *const grid);
static void NS_surface_BNS_CubedSpherical_grid(Grid_T *const grid);
