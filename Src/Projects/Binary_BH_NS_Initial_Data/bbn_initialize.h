#include "core_lib.h"
#include "TOV_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"
#include "memory_managing_lib.h"
#include "utilities_lib.h"

Grid_T *bbn_initialize_next_grid(Grid_T *const grid_prev);
static Grid_T *TOV_KerrShild_approximation(void);
static Grid_T *creat_grid_TOV_KerrShild(const double R_NS_l,const double R_BH_r,const double a_BH);
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,const double R_NS_l,const double R_BH_r,const double a_BH);
static void init_field_TOV_plus_KerrSchild(Grid_T *const grid);
static unsigned IsNSPatch(const Patch_T *const patch);
static void create_fields(Grid_T *const grid);




