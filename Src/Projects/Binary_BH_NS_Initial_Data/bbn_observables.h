#include "core_lib.h"
#include "physics_observables_lib.h"
#include "utilities_lib.h"
#include "manifold_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "fields_lib.h"
#include "bbn_headers.h"

void bbn_plan_PsJs_ADM_CS(Observable_T *obs);
void bbn_free_PsJs_ADM_CS(Observable_T *obs);
void bbn_populate_ADM_integrand_PdS_GdV(const Observable_T *const obs);
static double ADM_momentum_x_BBN_CS(Observable_T *const obs);
static double ADM_momentum_y_BBN_CS(Observable_T *const obs);
static double ADM_momentum_z_BBN_CS(Observable_T *const obs);
static double ADM_angular_momentum_x_BBN_CS(Observable_T *const obs);
static double ADM_angular_momentum_y_BBN_CS(Observable_T *const obs);
static double ADM_angular_momentum_z_BBN_CS(Observable_T *const obs);
static void populate_normal_surrounding(struct PsJs_ADM_S *const adm,const Dd_T dir);




