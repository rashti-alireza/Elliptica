#include "star_header.h"


typedef void fAdjustment_t (Physics_T *const phys);

int star_NS_idealfluid_gConf_find_Euler_const(Physics_T *const phys);
int star_NS_idealfluid_extrapolate_matter_fields(Physics_T *const phys);
int star_NS_idealfluid_gConf_force_balance(Physics_T *const phys);
static double Euler_eq_const_gConf_rootfinder_eq(void *params,const double *const x);
static void W_spin_vector_idealfluid(Patch_T *const patch,const double Omega_NS[3],const double C_NS[3]);
void star_W_spin_vector_idealfluid_update(Physics_T *const phys);
static void parse_adjust_parameter(const char *const par,char *adjust[3]);
static fAdjustment_t *get_func_force_balance_adjustment(const char *const adjust);
static void force_balance_ddx_x_CM(Physics_T *const phys);
static void force_balance_ddy_x_CM(Physics_T *const phys);
static void force_balance_ddz_x_CM(Physics_T *const phys);
static void force_balance_ddx_y_CM(Physics_T *const phys);
static void force_balance_ddy_y_CM(Physics_T *const phys);
static void force_balance_ddz_y_CM(Physics_T *const phys);
static void force_balance_ddx_Omega(Physics_T *const phys);
static void force_balance_ddy_Omega(Physics_T *const phys);
static void force_balance_ddz_Omega(Physics_T *const phys);
static void force_balance_ddCM_Omega(Physics_T *const phys);
static void force_balance_eq_root_finders(Physics_T *const phys,const int dir, const char *const par);


