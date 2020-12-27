#include "obs_header.h"

/* handy comparison */
#define IFss(X)    if(strstr_i(obs->quantity,X))
#define IFsc(X)    if(strcmp_i(obs->quantity,X))

/* is X set to Y? note: X is prefixed with physics */
#define IsIt(X,Y)  (strcmp_i(Gets(X),Y))

/* handy macro for collecting patches1 and patches2 depending on split */
#define Collect_outermost_S_V_SplitCS \
  /* for 1 split */ \
  if (Pgeti("grid_SplitCS_Nsplit_c") == 1) \
  { \
    /* surface part */ \
    region   = "outermost_OB"; \
    patches2 = collect_patches(grid,region,&N2); \
  } \
  else \
  { \
    /* volume part */ \
    region   = "outermost_OB"; \
    patches1 = collect_patches(grid,region,&N1); \
    /* surface part */ \
    region   = "outermost_OB"; \
    patches2 = collect_patches(grid,region,&N2); \
  }

/* set outermost surface integral flags depending on split */
#define Set_outermost_integral_S_SplitCS \
  adm[n]->surface_integration_flg = 1; \
  adm[n]->Z_surface = 1; \
  /* for 1 split */ \
  if (Pgeti("grid_SplitCS_Nsplit_c") == 1) \
  { \
    adm[n]->K = patch->n[2]-1; \
  } \
  else \
  { \
    adm[n]->K = 0; \
  } \
  n_physical_metric_around(adm[n],_c_);


void obs_calculate(Observe_T *const obs);
static void define_spin_campanelli(Observe_T *const obs);
static void define_spin_JRP(Observe_T *const obs);
static void define_spin_akv(Observe_T *const obs);
static void Rc_BH(Observe_T *const obs);
static void n_physical_metric_around(struct items_S *const adm,const Dd_T dir);
static void n_conformal_metric_around(struct items_S *const adm,const Dd_T dir);
static void calc_ADM_mass(Observe_T *const obs);
static void calc_ADM_J(Observe_T *const obs);
static void calc_ADM_P(Observe_T *const obs);
static void calc_Kommar_mass(Observe_T *const obs);
static void calc_irreducible_BH_mass(Observe_T *const obs);
double obs_integral_SV (Observe_T *const obs,
                        const char *const sS/* integrand for S */,
                        const char *const sV/* intergrand for V */,
                        const char sign_sS/* [+/-] integral of S */,
                        const char sign_sV/* [+/-] integral of V */);

