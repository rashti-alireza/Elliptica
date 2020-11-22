#include "obs_header.h"

void obs_plan(Observable_T *obs);
void obs_free(Observable_T *obs);
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
static void n_physical_metric_around(struct items_S *const adm,const Dd_T dir);
static void n_conformal_metric_around(struct items_S *const adm,const Dd_T dir);
void obs_Rc_BH(double Rc[3],Grid_T *const grid);





