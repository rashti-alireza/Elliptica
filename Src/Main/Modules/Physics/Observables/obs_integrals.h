#include "core_lib.h"
#include "physics_observables_lib.h"
#include "utilities_lib.h"
#include "manifold_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "fields_lib.h"
#include "obs_headers.h"

void obs_plan_obs_CS(Observable_T *obs);
void obs_free_obs_CS(Observable_T *obs);
void obs_populate_ADM_integrand_PdS_GdV_binary(const Observable_T *const obs);
void obs_populate_ADM_integrand_PdS_GdV_single(const Observable_T *const obs);
double obs_Kommar_mass(Observable_T *const obs);
double obs_ADM_mass(Observable_T *const obs);
void obs_define_spin_integral(double S[3],Grid_T *const grid,const char *const kind);
void obs_define_spin_JRP(double S[3],Grid_T *const grid,const char *const kind);
void obs_define_spin_akv(double S[3],Grid_T *const grid,const char *const kind);
static double ADM_momentum_x_BHNS_CS(Observable_T *const obs);
static double ADM_momentum_y_BHNS_CS(Observable_T *const obs);
static double ADM_momentum_z_BHNS_CS(Observable_T *const obs);
static double ADM_angular_momentum_x_BHNS_CS(Observable_T *const obs);
static double ADM_angular_momentum_y_BHNS_CS(Observable_T *const obs);
static double ADM_angular_momentum_z_BHNS_CS(Observable_T *const obs);
static void n_physical_metric_surrounding(struct items_S *const adm,const Dd_T dir);
static void n_conformal_metric_surrounding(struct items_S *const adm,const Dd_T dir);




