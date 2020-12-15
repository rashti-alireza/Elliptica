#include "obs_header.h"

/* handy comparison */
#define IFss(X)  if(strstr_i(obs->quantity,X))
#define IFsc(X)  if(strcmp_i(obs->quantity,X))

/* put it to 1 if you want \int{Gdv} */
#define VOLUME_INTEGRAL (1)

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
static void calc_ADM_mass(Observe_T *const obs);
static void calc_ADM_PJ(Observe_T *const obs);
static void calc_Kommar_mass(Observe_T *const obs);





