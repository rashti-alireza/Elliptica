#include "bbn_headers.h"



/* function prefix DON NOT change this prefix macro some fields 
// using this prefix which then you must change them too.*/
#define KS_func_prefix(name) bbn_ks_##name

/* function names and prototype */
#define KS_func_def_macro(name) double KS_func_prefix(name)

/* function args */
#define KS_func_args_macro  \
  (const double x,const double y, const double z, \
   const double M_BH,const double a_BH,const double phix, \
   const double phiy,const double Bx, \
   const double By,const double Bz,const double B2)

/* pass special argument to each function */   
#define KS_func_pass_args_macro  \
  (x,y,z,M_BH,a_BH,phix,phiy,Bx,By,Bz,B2)

/* pass arguments for the following functions */
#define bbn_ks_k0(x,y,z) bbn_ks_k0 KS_func_pass_args_macro
#define bbn_ks_k1(x,y,z) bbn_ks_k1 KS_func_pass_args_macro
#define bbn_ks_k2(x,y,z) bbn_ks_k2 KS_func_pass_args_macro
#define bbn_ks_kt(x,y,z) bbn_ks_kt KS_func_pass_args_macro
#define bbn_ks_c(x,y,z)  bbn_ks_c  KS_func_pass_args_macro


/* all external functions */
void bbn_free_data_g_gI_analytic(
        Patch_T *const patch,
        double *(*get_v)(const char *const fname,void *params),
        void *params);
void bbn_free_data_dg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params);
void bbn_free_data_ddg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params);
void bbn_free_data_dddg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params);

KS_func_def_macro(K0)KS_func_args_macro;
KS_func_def_macro(dK0_D1)KS_func_args_macro;
KS_func_def_macro(dK0_D0)KS_func_args_macro;
KS_func_def_macro(dK0_D2)KS_func_args_macro;
KS_func_def_macro(ddK0_D1D2)KS_func_args_macro;
KS_func_def_macro(ddK0_D0D0)KS_func_args_macro;
KS_func_def_macro(ddK0_D2D2)KS_func_args_macro;
KS_func_def_macro(ddK0_D1D1)KS_func_args_macro;
KS_func_def_macro(ddK0_D0D2)KS_func_args_macro;
KS_func_def_macro(ddK0_D0D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D2D2D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D2D2D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D0D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D1D2)KS_func_args_macro;
KS_func_def_macro(K1)KS_func_args_macro;
KS_func_def_macro(dK1_D2)KS_func_args_macro;
KS_func_def_macro(dK1_D0)KS_func_args_macro;
KS_func_def_macro(dK1_D1)KS_func_args_macro;
KS_func_def_macro(ddK1_D0D1)KS_func_args_macro;
KS_func_def_macro(ddK1_D1D1)KS_func_args_macro;
KS_func_def_macro(ddK1_D2D2)KS_func_args_macro;
KS_func_def_macro(ddK1_D0D2)KS_func_args_macro;
KS_func_def_macro(ddK1_D1D2)KS_func_args_macro;
KS_func_def_macro(ddK1_D0D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D2D2D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D2D2D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D1D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D0D1)KS_func_args_macro;
KS_func_def_macro(K2)KS_func_args_macro;
KS_func_def_macro(dK2_D0)KS_func_args_macro;
KS_func_def_macro(dK2_D1)KS_func_args_macro;
KS_func_def_macro(dK2_D2)KS_func_args_macro;
KS_func_def_macro(ddK2_D1D2)KS_func_args_macro;
KS_func_def_macro(ddK2_D0D1)KS_func_args_macro;
KS_func_def_macro(ddK2_D2D2)KS_func_args_macro;
KS_func_def_macro(ddK2_D1D1)KS_func_args_macro;
KS_func_def_macro(ddK2_D0D0)KS_func_args_macro;
KS_func_def_macro(ddK2_D0D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D0D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D1D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D2D2D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D2D2D2)KS_func_args_macro;
KS_func_def_macro(H)KS_func_args_macro;
KS_func_def_macro(dH_D1)KS_func_args_macro;
KS_func_def_macro(dH_D0)KS_func_args_macro;
KS_func_def_macro(dH_D2)KS_func_args_macro;
KS_func_def_macro(ddH_D1D1)KS_func_args_macro;
KS_func_def_macro(ddH_D0D2)KS_func_args_macro;
KS_func_def_macro(ddH_D2D2)KS_func_args_macro;
KS_func_def_macro(ddH_D0D1)KS_func_args_macro;
KS_func_def_macro(ddH_D0D0)KS_func_args_macro;
KS_func_def_macro(ddH_D1D2)KS_func_args_macro;
KS_func_def_macro(dddH_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddH_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddH_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddH_D1D1D2)KS_func_args_macro;
KS_func_def_macro(dddH_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddH_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddH_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddH_D2D2D2)KS_func_args_macro;
KS_func_def_macro(dddH_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddH_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddH_D0D0D1)KS_func_args_macro;
KS_func_def_macro(dddH_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddH_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddH_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddH_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddH_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddH_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddH_D2D2D1)KS_func_args_macro;

KS_func_def_macro(k0) KS_func_args_macro;
KS_func_def_macro(dk0_D2) KS_func_args_macro;
KS_func_def_macro(dk0_D1) KS_func_args_macro;
KS_func_def_macro(dk0_D0) KS_func_args_macro;
KS_func_def_macro(ddk0_D1D2) KS_func_args_macro;
KS_func_def_macro(ddk0_D0D0) KS_func_args_macro;
KS_func_def_macro(ddk0_D1D1) KS_func_args_macro;
KS_func_def_macro(ddk0_D0D2) KS_func_args_macro;
KS_func_def_macro(ddk0_D0D1) KS_func_args_macro;
KS_func_def_macro(ddk0_D2D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D2D2) KS_func_args_macro;
KS_func_def_macro(k1) KS_func_args_macro;
KS_func_def_macro(dk1_D2) KS_func_args_macro;
KS_func_def_macro(dk1_D0) KS_func_args_macro;
KS_func_def_macro(dk1_D1) KS_func_args_macro;
KS_func_def_macro(ddk1_D1D1) KS_func_args_macro;
KS_func_def_macro(ddk1_D1D2) KS_func_args_macro;
KS_func_def_macro(ddk1_D2D2) KS_func_args_macro;
KS_func_def_macro(ddk1_D0D2) KS_func_args_macro;
KS_func_def_macro(ddk1_D0D0) KS_func_args_macro;
KS_func_def_macro(ddk1_D0D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D1D2) KS_func_args_macro;
KS_func_def_macro(k2) KS_func_args_macro;
KS_func_def_macro(dk2_D2) KS_func_args_macro;
KS_func_def_macro(dk2_D1) KS_func_args_macro;
KS_func_def_macro(dk2_D0) KS_func_args_macro;
KS_func_def_macro(ddk2_D0D1) KS_func_args_macro;
KS_func_def_macro(ddk2_D0D0) KS_func_args_macro;
KS_func_def_macro(ddk2_D0D2) KS_func_args_macro;
KS_func_def_macro(ddk2_D2D2) KS_func_args_macro;
KS_func_def_macro(ddk2_D1D1) KS_func_args_macro;
KS_func_def_macro(ddk2_D1D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D1D0) KS_func_args_macro;
KS_func_def_macro(kt) KS_func_args_macro;
KS_func_def_macro(dkt_D1) KS_func_args_macro;
KS_func_def_macro(dkt_D2) KS_func_args_macro;
KS_func_def_macro(dkt_D0) KS_func_args_macro;
KS_func_def_macro(ddkt_D2D2) KS_func_args_macro;
KS_func_def_macro(ddkt_D0D0) KS_func_args_macro;
KS_func_def_macro(ddkt_D0D1) KS_func_args_macro;
KS_func_def_macro(ddkt_D1D1) KS_func_args_macro;
KS_func_def_macro(ddkt_D0D2) KS_func_args_macro;
KS_func_def_macro(ddkt_D1D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D1D0) KS_func_args_macro;
KS_func_def_macro(c) KS_func_args_macro;
KS_func_def_macro(dc_D2) KS_func_args_macro;
KS_func_def_macro(dc_D1) KS_func_args_macro;
KS_func_def_macro(dc_D0) KS_func_args_macro;
KS_func_def_macro(ddc_D0D0) KS_func_args_macro;
KS_func_def_macro(ddc_D2D2) KS_func_args_macro;
KS_func_def_macro(ddc_D0D1) KS_func_args_macro;
KS_func_def_macro(ddc_D0D2) KS_func_args_macro;
KS_func_def_macro(ddc_D1D2) KS_func_args_macro;
KS_func_def_macro(ddc_D1D1) KS_func_args_macro;
KS_func_def_macro(dddc_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddc_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddc_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddc_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddc_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddc_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddc_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddc_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddc_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddc_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddc_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddc_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddc_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddc_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddc_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddc_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddc_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddc_D2D2D1) KS_func_args_macro;

KS_func_def_macro(X) KS_func_args_macro;
KS_func_def_macro(dX_D0) KS_func_args_macro;
KS_func_def_macro(dX_D2) KS_func_args_macro;
KS_func_def_macro(dX_D1) KS_func_args_macro;
KS_func_def_macro(ddX_D1D1) KS_func_args_macro;
KS_func_def_macro(ddX_D0D2) KS_func_args_macro;
KS_func_def_macro(ddX_D0D0) KS_func_args_macro;
KS_func_def_macro(ddX_D1D2) KS_func_args_macro;
KS_func_def_macro(ddX_D0D1) KS_func_args_macro;
KS_func_def_macro(ddX_D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddX_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddX_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddX_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddX_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddX_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddX_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddX_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddX_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddX_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddX_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddX_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddX_D1D1D0) KS_func_args_macro;
KS_func_def_macro(Y) KS_func_args_macro;
KS_func_def_macro(dY_D0) KS_func_args_macro;
KS_func_def_macro(dY_D2) KS_func_args_macro;
KS_func_def_macro(dY_D1) KS_func_args_macro;
KS_func_def_macro(ddY_D1D2) KS_func_args_macro;
KS_func_def_macro(ddY_D0D1) KS_func_args_macro;
KS_func_def_macro(ddY_D1D1) KS_func_args_macro;
KS_func_def_macro(ddY_D0D2) KS_func_args_macro;
KS_func_def_macro(ddY_D2D2) KS_func_args_macro;
KS_func_def_macro(ddY_D0D0) KS_func_args_macro;
KS_func_def_macro(dddY_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddY_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddY_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddY_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddY_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddY_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddY_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddY_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddY_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddY_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddY_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddY_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddY_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddY_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D1D1) KS_func_args_macro;
KS_func_def_macro(Z) KS_func_args_macro;
KS_func_def_macro(dZ_D2) KS_func_args_macro;
KS_func_def_macro(dZ_D1) KS_func_args_macro;
KS_func_def_macro(dZ_D0) KS_func_args_macro;
KS_func_def_macro(ddZ_D2D2) KS_func_args_macro;
KS_func_def_macro(ddZ_D0D0) KS_func_args_macro;
KS_func_def_macro(ddZ_D0D1) KS_func_args_macro;
KS_func_def_macro(ddZ_D1D2) KS_func_args_macro;
KS_func_def_macro(ddZ_D1D1) KS_func_args_macro;
KS_func_def_macro(ddZ_D0D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D2D0) KS_func_args_macro;
KS_func_def_macro(R) KS_func_args_macro;
KS_func_def_macro(dR_D2) KS_func_args_macro;
KS_func_def_macro(dR_D1) KS_func_args_macro;
KS_func_def_macro(dR_D0) KS_func_args_macro;
KS_func_def_macro(ddR_D1D2) KS_func_args_macro;
KS_func_def_macro(ddR_D1D1) KS_func_args_macro;
KS_func_def_macro(ddR_D2D2) KS_func_args_macro;
KS_func_def_macro(ddR_D0D2) KS_func_args_macro;
KS_func_def_macro(ddR_D0D0) KS_func_args_macro;
KS_func_def_macro(ddR_D0D1) KS_func_args_macro;
KS_func_def_macro(dddR_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddR_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddR_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddR_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddR_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddR_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddR_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddR_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddR_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddR_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddR_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddR_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddR_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddR_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddR_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddR_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddR_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddR_D1D2D0) KS_func_args_macro;

