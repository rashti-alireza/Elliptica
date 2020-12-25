#include "obs_header.h"

/* handy comparison */
#define IFss(X)    if(strstr_i(obs->quantity,X))
#define IFsc(X)    if(strcmp_i(obs->quantity,X))

/* is X set to Y? note: X is prefixed with physics */
#define IsIt(X,Y)  Pcmps(MyParam(X),Y)

/* put it to 1 if you want \int{Gdv}
// this is mainly for test purposes. */
#define VOLUME_INTEGRAL (1)

void obs_calculate(Observe_T *const obs);
static void define_spin_campanelli(Observe_T *const obs);
static void define_spin_JRP(Observe_T *const obs);
static void define_spin_akv(Observe_T *const obs);
static void Rc_BH(Observe_T *const obs);
static void n_physical_metric_around(struct items_S *const adm,const Dd_T dir);
static void n_conformal_metric_around(struct items_S *const adm,const Dd_T dir);
static void calc_ADM_mass(Observe_T *const obs);
static void calc_ADM_PJ(Observe_T *const obs);
static void calc_Kommar_mass(Observe_T *const obs);
static double integral_ADM_PJ(Observe_T *const obs,
                              const char *const sP/* integrand for S */,
                              const char *const sG/* intergrand for V */);
static void calc_irreducible_BH_mass(Observe_T *const obs);


