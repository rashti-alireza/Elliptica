#include "obs_header.h"

void obs_plan(Observable_T *obs);
void obs_define_spin_integral(double S[3],Physics_T *const phys);
void obs_define_spin_JRP(double S[3],Physics_T *const phys);
void obs_define_spin_akv(double S[3],Physics_T *const phys);
void obs_Rc_BH(double Rc[3],Physics_T *const phys);
static double ADM_momentum_x_BHNS_CS(Observable_T *const obs);
static double ADM_momentum_y_BHNS_CS(Observable_T *const obs);
static double ADM_momentum_z_BHNS_CS(Observable_T *const obs);
static double ADM_angular_momentum_x_BHNS_CS(Observable_T *const obs);
static double ADM_angular_momentum_y_BHNS_CS(Observable_T *const obs);
static double ADM_angular_momentum_z_BHNS_CS(Observable_T *const obs);
static void n_physical_metric_around(struct items_S *const adm,const Dd_T dir);
static void n_conformal_metric_around(struct items_S *const adm,const Dd_T dir);





