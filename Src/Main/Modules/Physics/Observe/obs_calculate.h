#include "obs_header.h"
#include "physics_star_lib.h"


/* string len */
#define obs_STR_LEN (999)

/* handy comparison */
#define IFss(X)    if(strstr_i(obs->quantity,X))
#define IFsc(X)    if(strcmp_i(obs->quantity,X))

/* is method set to X? */
#define IsIt(X)  (strcmp_i(obs->method,X))

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
#define Set_outermost_integral_S_SplitCS(item) \
  item[n]->surface_integration_flg = 1; \
  item[n]->Z_surface = 1; \
  /* for 1 split */ \
  if (Pgeti("grid_SplitCS_Nsplit_c") == 1) \
  { \
    item[n]->K = patch->n[2]-1; \
  } \
  else \
  { \
    item[n]->K = 0; \
  }

/* handy error msg */
#define SET_MSG  \
 char obs_err_msg[obs_STR_LEN]; \
 sprintf(obs_err_msg,"In '%s' NO such '%s' or '%s' found!\n", \
         __func__,obs->quantity,obs->method);


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
static void calc_CM(Observe_T *const obs);
static void calc_spin(Observe_T *const obs);
static void calc_baryonic_mass(Observe_T *const obs);


