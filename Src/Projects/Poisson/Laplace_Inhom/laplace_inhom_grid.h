#include "core_lib.h"
#include "memory_managing_lib.h"
#include "coordinates_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"

Grid_T *Laplace_Inhom_make_grid(void);
static void grid_characteristics_Laplace_Inhome (Grid_T *const grid);
static void characteristics_Cartesian_grid(Grid_T *const grid);
static void characteristics_BNS_Projective_grid(Grid_T *const grid);
static void make_field_of_NS_radius(Grid_T *const grid,double *const R_max_l,double *const R_max_r);

