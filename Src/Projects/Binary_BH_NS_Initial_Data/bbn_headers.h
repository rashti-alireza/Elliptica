#include "core_lib.h"
#include "TOV_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"
#include "memory_managing_lib.h"
#include "utilities_lib.h"
#include "maths_approximation_lib.h"
#include "physics_EoS_lib.h"
#include "physics_StressEnergyTensor_lib.h"

void bbn_update_psi10A_UiUj(Grid_T *const grid);
void bbn_populate_free_data(Grid_T *const grid);
void bbn_free_data_gammas(Grid_T *const grid);
void bbn_free_conformal_metric_derivatives(Patch_T *const patch);
void bbn_preparing_conformal_metric_derivatives(Patch_T *const patch);
void bbn_free_data_Gamma(Grid_T *const grid);
void bbn_free_data_Ricci(Grid_T *const grid);
void bbn_free_data_KS_trKij(Patch_T *const patch);
