#include "core_lib.h"
#include "TOV_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"
#include "memory_managing_lib.h"
#include "utilities_lib.h"
#include "maths_approximation_lib.h"
#include "physics_EoS_lib.h"
#include "physics_StressEnergyTensor_lib.h"

Grid_T *bbn_initialize_next_grid(Grid_T *const grid_prev);
static Grid_T *TOV_KerrShild_approximation(void);
static Grid_T *creat_grid_TOV_KerrShild(const double R_NS_l,const double R_BH_r,const double a_BH);
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,const double R_NS_l,const double R_BH_r,const double a_BH);
static void init_field_TOV_plus_KerrSchild(Grid_T *const grid,const TOV_T *const tov, const double a_BH, const double M_BH);
void bbn_allocate_fields(Grid_T *const grid);
void bbn_partial_derivatives_fields(Grid_T *const grid);
void bbn_populate_free_data(Grid_T *const grid);
void bbn_update_psi10A_UiUj(Grid_T *const grid);



