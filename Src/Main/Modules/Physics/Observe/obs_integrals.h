#include "obs_header.h"

void obs_calculate(Observe_T *const obs);
static void define_spin_campanelli(Observe_T *const obs);
static void define_spin_JRP(Observe_T *const obs);
static void define_spin_akv(Observe_T *const obs);
static void Rc_BH(Observe_T *const obs);
static double ADM_momentum_x_BHNS_CS(Observe_T *const obs);
static double ADM_momentum_y_BHNS_CS(Observe_T *const obs);
static double ADM_momentum_z_BHNS_CS(Observe_T *const obs);
static double ADM_angular_momentum_x_BHNS_CS(Observe_T *const obs);
static double ADM_angular_momentum_y_BHNS_CS(Observe_T *const obs);
static double ADM_angular_momentum_z_BHNS_CS(Observe_T *const obs);
static void n_physical_metric_around(struct items_S *const adm,const Dd_T dir);
static void n_conformal_metric_around(struct items_S *const adm,const Dd_T dir);





