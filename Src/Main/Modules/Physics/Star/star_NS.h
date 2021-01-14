#include "star_header.h"


typedef void fAdjustment_t (Physics_T *const phys);

int star_NS_ifluid_gConf_find_EulerC_fix_baryon_mass(Physics_T *const phys);
int star_NS_idealfluid_extrapolate_matter_fields(Physics_T *const phys);
int star_NS_idealfluid_gConf_force_balance(Physics_T *const phys);
static double Euler_eq_const_gConf_rootfinder_eq(void *params,const double *const x);
static void W_spin_vector_idealfluid(Patch_T *const patch,const double Omega_NS[3],const double C_NS[3]);
void star_W_spin_vector_idealfluid_update(Physics_T *const phys,const char *const region);
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
static double dh_dx0_root_finder_eq(void *params,const double *const x);
static double dh_dx1_root_finder_eq(void *params,const double *const x);
static double dh_dx2_root_finder_eq(void *params,const double *const x);
void star_NS_find_where_denthalpy_is_0(Physics_T *const phys,double xdh0[3]);
int star_NS_keep_center_fixed(Physics_T *const phys);
static void adjust_NS_center_interpolation(Physics_T *const phys);
static void adjust_NS_center_Taylor_expansion(Physics_T *const phys);
void star_populate_psi_alphaPsi_matter_fields_TOV
      (Physics_T *const phys,const char *const region,
      const char *const Psi,const char *const AlphaPsi,
      const char *const Enthalpy,const char *const Rho0,
      const char *const Phi,const char *const W);

  




