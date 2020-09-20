#include "bbn_headers.h"


/* mathematica */
#define Cos(a) cos(a)
#define Sin(a) sin(a)
#define Tan(a) tan(a)
#define Cot(a) (1./tan(a))
#define Sec(a) (1./cos(a))
#define Csc(a) (1./sin(a))
#define Power(a,b) pow(a,b)
#define Sqrt(a)    sqrt(a)

/* nothing for Hold */
#define Hold    

/* nothing for List */
#define List  

/* concatenate */
#define Pattern(a,b) a##_##b


/* function prefix DON NOT change this prefix macro some fields 
// using this prefix which then you must change them too.*/
#define KS_func_prefix(name) bbn_ks_##name

/* function names and prototype */
#define KS_func_def_macro(name) double KS_func_prefix(name)

/* function args */
#define KS_func_args_macro  \
  (const double x,const double y, const double z,/* Cartesian coords */ \
   const double M_BH,const double a_BH,/* mass and spin of BH */ \
   const double phiy, const double phiz,/* rotation angels */\
   const double Bx, const double By,const double Bz,/* B=v/c (boost) */\
   const double B2/* B^i B_i */,\
   const double r0/* roll off radius */,\
   const double Lambda/* flat data => 0, kerr-schild => 1 */)

/* pass special argument to each function */   
#define KS_func_pass_args_macro  \
  (x,y,z,M_BH,a_BH,phiy,phiz,Bx,By,Bz,B2,r0,Lambda)

/* derivative function args */
#define KS_deriv_func_args_macro  \
  (const char *const stem/* ex. bbn_ks_k */,\
   const char *const derivs/* ex."x,y,z" */,\
   const double x,const double y, const double z,/* Cartesian coords */ \
   const double M_BH,const double a_BH,/* mass and spin of BH */ \
   const double phiy, const double phiz,/* rotation angels */\
   const double Bx, const double By,const double Bz,/* B=v/c (boost) */\
   const double B2/* B^i B_i */,\
   const double r0/* roll off radius */,\
   const double Lambda/* flat data => 0, kerr-schild => 1 */)

/* derivative */
#define Derivative(a,...) \
  bbn_ks_derivative \
  (#a,#__VA_ARGS__,x,y,z, \
   M_BH,a_BH,phiy,phiz,Bx,By,Bz,\
   B2,r0,Lambda)
  

/* derivative */
#define D(a,...) \
  bbn_ks_derivative \
  (#a,#__VA_ARGS__,x,y,z, \
   M_BH,a_BH,phiy,phiz,Bx,By,Bz,\
   B2,r0,Lambda)

/* bbn_ks type def function. FIND the derivative data base in the 
// bottom of this file. this is used for various derivative functions. */
typedef double Fbbn_ks_func_t KS_func_args_macro;

/* pass arguments for the following functions */
#define bbn_ks_k0(x,y,z) bbn_ks_k0 KS_func_pass_args_macro
#define bbn_ks_k1(x,y,z) bbn_ks_k1 KS_func_pass_args_macro
#define bbn_ks_k2(x,y,z) bbn_ks_k2 KS_func_pass_args_macro
#define bbn_ks_kt(x,y,z) bbn_ks_kt KS_func_pass_args_macro
#define bbn_ks_c(x,y,z)  bbn_ks_c  KS_func_pass_args_macro

#define bbn_ks_K0(x,y,z) bbn_ks_K0 KS_func_pass_args_macro
#define bbn_ks_K1(x,y,z) bbn_ks_K1 KS_func_pass_args_macro
#define bbn_ks_K2(x,y,z) bbn_ks_K2 KS_func_pass_args_macro


#define bbn_ks_X(x,y,z)  bbn_ks_X  KS_func_pass_args_macro
#define bbn_ks_Y(x,y,z)  bbn_ks_Y  KS_func_pass_args_macro
#define bbn_ks_Z(x,y,z)  bbn_ks_Z  KS_func_pass_args_macro
#define bbn_ks_R(x,y,z)  bbn_ks_R  KS_func_pass_args_macro
#define bbn_ks_H(x,y,z)  bbn_ks_H  KS_func_pass_args_macro

#define bbn_ks_rolloff(x,y,z)  bbn_ks_rolloff KS_func_pass_args_macro

/* replace this functions */
#define XX bbn_ks_X
#define YY bbn_ks_Y
#define ZZ bbn_ks_Z


/* avoid compiler warning since not all of functions need all variables */
#define ks_no_use_macro \
  UNUSED(0*M_BH*a_BH*phiy*phiz*Bx*By*Bz*B2*r0*Lambda)

/* all external functions */
double bbn_ks_derivative KS_deriv_func_args_macro;

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

KS_func_def_macro(rolloff) KS_func_args_macro;
KS_func_def_macro(drolloff_D0) KS_func_args_macro;
KS_func_def_macro(drolloff_D2) KS_func_args_macro;
KS_func_def_macro(drolloff_D1) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D2D2) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D1D2) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D0D1) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D0D2) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D0D0) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D1D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D0D1) KS_func_args_macro;

/* macro needed for db macros */
#define StrinIt(x) #x
#define KS_func_db_macro(x) KS_func_prefix(x)
#define KS_name_db_macro(x) StrinIt(KS_func_prefix(x))

/* derivative function data base. 
// NOTE: derive_func_db MUST match with derive_name_db below. */
static Fbbn_ks_func_t *const derive_func_db[] = 
{

KS_func_db_macro(K0),
KS_func_db_macro(dK0_D1),
KS_func_db_macro(dK0_D0),
KS_func_db_macro(dK0_D2),
KS_func_db_macro(ddK0_D1D2),
KS_func_db_macro(ddK0_D0D0),
KS_func_db_macro(ddK0_D2D2),
KS_func_db_macro(ddK0_D1D1),
KS_func_db_macro(ddK0_D0D2),
KS_func_db_macro(ddK0_D0D1),
KS_func_db_macro(dddK0_D2D2D1),
KS_func_db_macro(dddK0_D2D2D0),
KS_func_db_macro(dddK0_D0D2D1),
KS_func_db_macro(dddK0_D0D1D1),
KS_func_db_macro(dddK0_D2D2D2),
KS_func_db_macro(dddK0_D1D2D2),
KS_func_db_macro(dddK0_D0D2D2),
KS_func_db_macro(dddK0_D1D2D1),
KS_func_db_macro(dddK0_D0D0D0),
KS_func_db_macro(dddK0_D0D0D1),
KS_func_db_macro(dddK0_D1D1D0),
KS_func_db_macro(dddK0_D1D2D0),
KS_func_db_macro(dddK0_D0D2D0),
KS_func_db_macro(dddK0_D0D1D0),
KS_func_db_macro(dddK0_D0D1D2),
KS_func_db_macro(dddK0_D0D0D2),
KS_func_db_macro(dddK0_D1D1D1),
KS_func_db_macro(dddK0_D1D1D2),
KS_func_db_macro(K1),
KS_func_db_macro(dK1_D2),
KS_func_db_macro(dK1_D0),
KS_func_db_macro(dK1_D1),
KS_func_db_macro(ddK1_D0D1),
KS_func_db_macro(ddK1_D1D1),
KS_func_db_macro(ddK1_D2D2),
KS_func_db_macro(ddK1_D0D2),
KS_func_db_macro(ddK1_D1D2),
KS_func_db_macro(ddK1_D0D0),
KS_func_db_macro(dddK1_D0D2D1),
KS_func_db_macro(dddK1_D1D1D0),
KS_func_db_macro(dddK1_D2D2D2),
KS_func_db_macro(dddK1_D1D2D0),
KS_func_db_macro(dddK1_D1D1D1),
KS_func_db_macro(dddK1_D0D1D2),
KS_func_db_macro(dddK1_D0D2D2),
KS_func_db_macro(dddK1_D0D0D0),
KS_func_db_macro(dddK1_D2D2D0),
KS_func_db_macro(dddK1_D2D2D1),
KS_func_db_macro(dddK1_D0D2D0),
KS_func_db_macro(dddK1_D0D0D2),
KS_func_db_macro(dddK1_D1D2D1),
KS_func_db_macro(dddK1_D1D1D2),
KS_func_db_macro(dddK1_D0D1D0),
KS_func_db_macro(dddK1_D1D2D2),
KS_func_db_macro(dddK1_D0D1D1),
KS_func_db_macro(dddK1_D0D0D1),
KS_func_db_macro(K2),
KS_func_db_macro(dK2_D0),
KS_func_db_macro(dK2_D1),
KS_func_db_macro(dK2_D2),
KS_func_db_macro(ddK2_D1D2),
KS_func_db_macro(ddK2_D0D1),
KS_func_db_macro(ddK2_D2D2),
KS_func_db_macro(ddK2_D1D1),
KS_func_db_macro(ddK2_D0D0),
KS_func_db_macro(ddK2_D0D2),
KS_func_db_macro(dddK2_D0D2D0),
KS_func_db_macro(dddK2_D1D2D1),
KS_func_db_macro(dddK2_D0D1D0),
KS_func_db_macro(dddK2_D0D2D2),
KS_func_db_macro(dddK2_D0D0D1),
KS_func_db_macro(dddK2_D0D1D1),
KS_func_db_macro(dddK2_D1D2D2),
KS_func_db_macro(dddK2_D1D1D1),
KS_func_db_macro(dddK2_D0D1D2),
KS_func_db_macro(dddK2_D1D1D2),
KS_func_db_macro(dddK2_D1D2D0),
KS_func_db_macro(dddK2_D1D1D0),
KS_func_db_macro(dddK2_D2D2D0),
KS_func_db_macro(dddK2_D0D0D0),
KS_func_db_macro(dddK2_D0D2D1),
KS_func_db_macro(dddK2_D0D0D2),
KS_func_db_macro(dddK2_D2D2D1),
KS_func_db_macro(dddK2_D2D2D2),
KS_func_db_macro(H),
KS_func_db_macro(dH_D1),
KS_func_db_macro(dH_D0),
KS_func_db_macro(dH_D2),
KS_func_db_macro(ddH_D1D1),
KS_func_db_macro(ddH_D0D2),
KS_func_db_macro(ddH_D2D2),
KS_func_db_macro(ddH_D0D1),
KS_func_db_macro(ddH_D0D0),
KS_func_db_macro(ddH_D1D2),
KS_func_db_macro(dddH_D0D2D2),
KS_func_db_macro(dddH_D1D2D0),
KS_func_db_macro(dddH_D0D1D0),
KS_func_db_macro(dddH_D1D1D2),
KS_func_db_macro(dddH_D0D2D0),
KS_func_db_macro(dddH_D1D2D1),
KS_func_db_macro(dddH_D0D0D0),
KS_func_db_macro(dddH_D2D2D2),
KS_func_db_macro(dddH_D1D2D2),
KS_func_db_macro(dddH_D0D1D2),
KS_func_db_macro(dddH_D0D0D1),
KS_func_db_macro(dddH_D1D1D1),
KS_func_db_macro(dddH_D0D2D1),
KS_func_db_macro(dddH_D1D1D0),
KS_func_db_macro(dddH_D0D0D2),
KS_func_db_macro(dddH_D2D2D0),
KS_func_db_macro(dddH_D0D1D1),
KS_func_db_macro(dddH_D2D2D1),

KS_func_db_macro(k0) ,
KS_func_db_macro(dk0_D2) ,
KS_func_db_macro(dk0_D1) ,
KS_func_db_macro(dk0_D0) ,
KS_func_db_macro(ddk0_D1D2) ,
KS_func_db_macro(ddk0_D0D0) ,
KS_func_db_macro(ddk0_D1D1) ,
KS_func_db_macro(ddk0_D0D2) ,
KS_func_db_macro(ddk0_D0D1) ,
KS_func_db_macro(ddk0_D2D2) ,
KS_func_db_macro(dddk0_D2D2D0) ,
KS_func_db_macro(dddk0_D0D2D2) ,
KS_func_db_macro(dddk0_D0D1D1) ,
KS_func_db_macro(dddk0_D0D0D1) ,
KS_func_db_macro(dddk0_D0D2D0) ,
KS_func_db_macro(dddk0_D0D1D2) ,
KS_func_db_macro(dddk0_D1D1D2) ,
KS_func_db_macro(dddk0_D0D2D1) ,
KS_func_db_macro(dddk0_D1D1D0) ,
KS_func_db_macro(dddk0_D1D2D0) ,
KS_func_db_macro(dddk0_D0D1D0) ,
KS_func_db_macro(dddk0_D1D2D1) ,
KS_func_db_macro(dddk0_D0D0D2) ,
KS_func_db_macro(dddk0_D2D2D1) ,
KS_func_db_macro(dddk0_D1D1D1) ,
KS_func_db_macro(dddk0_D2D2D2) ,
KS_func_db_macro(dddk0_D0D0D0) ,
KS_func_db_macro(dddk0_D1D2D2) ,
KS_func_db_macro(k1) ,
KS_func_db_macro(dk1_D2) ,
KS_func_db_macro(dk1_D0) ,
KS_func_db_macro(dk1_D1) ,
KS_func_db_macro(ddk1_D1D1) ,
KS_func_db_macro(ddk1_D1D2) ,
KS_func_db_macro(ddk1_D2D2) ,
KS_func_db_macro(ddk1_D0D2) ,
KS_func_db_macro(ddk1_D0D0) ,
KS_func_db_macro(ddk1_D0D1) ,
KS_func_db_macro(dddk1_D0D1D2) ,
KS_func_db_macro(dddk1_D0D1D0) ,
KS_func_db_macro(dddk1_D0D2D1) ,
KS_func_db_macro(dddk1_D1D2D1) ,
KS_func_db_macro(dddk1_D2D2D2) ,
KS_func_db_macro(dddk1_D1D2D2) ,
KS_func_db_macro(dddk1_D2D2D1) ,
KS_func_db_macro(dddk1_D0D0D2) ,
KS_func_db_macro(dddk1_D1D2D0) ,
KS_func_db_macro(dddk1_D0D2D0) ,
KS_func_db_macro(dddk1_D1D1D1) ,
KS_func_db_macro(dddk1_D0D2D2) ,
KS_func_db_macro(dddk1_D0D1D1) ,
KS_func_db_macro(dddk1_D0D0D0) ,
KS_func_db_macro(dddk1_D0D0D1) ,
KS_func_db_macro(dddk1_D2D2D0) ,
KS_func_db_macro(dddk1_D1D1D0) ,
KS_func_db_macro(dddk1_D1D1D2) ,
KS_func_db_macro(k2) ,
KS_func_db_macro(dk2_D2) ,
KS_func_db_macro(dk2_D1) ,
KS_func_db_macro(dk2_D0) ,
KS_func_db_macro(ddk2_D0D1) ,
KS_func_db_macro(ddk2_D0D0) ,
KS_func_db_macro(ddk2_D0D2) ,
KS_func_db_macro(ddk2_D2D2) ,
KS_func_db_macro(ddk2_D1D1) ,
KS_func_db_macro(ddk2_D1D2) ,
KS_func_db_macro(dddk2_D0D1D2) ,
KS_func_db_macro(dddk2_D1D2D2) ,
KS_func_db_macro(dddk2_D1D1D2) ,
KS_func_db_macro(dddk2_D2D2D2) ,
KS_func_db_macro(dddk2_D1D1D0) ,
KS_func_db_macro(dddk2_D2D2D0) ,
KS_func_db_macro(dddk2_D0D2D1) ,
KS_func_db_macro(dddk2_D1D2D0) ,
KS_func_db_macro(dddk2_D0D0D1) ,
KS_func_db_macro(dddk2_D0D2D2) ,
KS_func_db_macro(dddk2_D2D2D1) ,
KS_func_db_macro(dddk2_D0D1D1) ,
KS_func_db_macro(dddk2_D0D2D0) ,
KS_func_db_macro(dddk2_D1D1D1) ,
KS_func_db_macro(dddk2_D1D2D1) ,
KS_func_db_macro(dddk2_D0D0D2) ,
KS_func_db_macro(dddk2_D0D0D0) ,
KS_func_db_macro(dddk2_D0D1D0) ,
KS_func_db_macro(kt) ,
KS_func_db_macro(dkt_D1) ,
KS_func_db_macro(dkt_D2) ,
KS_func_db_macro(dkt_D0) ,
KS_func_db_macro(ddkt_D2D2) ,
KS_func_db_macro(ddkt_D0D0) ,
KS_func_db_macro(ddkt_D0D1) ,
KS_func_db_macro(ddkt_D1D1) ,
KS_func_db_macro(ddkt_D0D2) ,
KS_func_db_macro(ddkt_D1D2) ,
KS_func_db_macro(dddkt_D0D1D1) ,
KS_func_db_macro(dddkt_D1D2D2) ,
KS_func_db_macro(dddkt_D0D1D2) ,
KS_func_db_macro(dddkt_D0D0D0) ,
KS_func_db_macro(dddkt_D0D2D0) ,
KS_func_db_macro(dddkt_D0D2D1) ,
KS_func_db_macro(dddkt_D0D2D2) ,
KS_func_db_macro(dddkt_D2D2D2) ,
KS_func_db_macro(dddkt_D1D2D0) ,
KS_func_db_macro(dddkt_D0D0D2) ,
KS_func_db_macro(dddkt_D1D1D0) ,
KS_func_db_macro(dddkt_D1D1D2) ,
KS_func_db_macro(dddkt_D2D2D1) ,
KS_func_db_macro(dddkt_D1D1D1) ,
KS_func_db_macro(dddkt_D2D2D0) ,
KS_func_db_macro(dddkt_D1D2D1) ,
KS_func_db_macro(dddkt_D0D0D1) ,
KS_func_db_macro(dddkt_D0D1D0) ,
KS_func_db_macro(c) ,
KS_func_db_macro(dc_D2) ,
KS_func_db_macro(dc_D1) ,
KS_func_db_macro(dc_D0) ,
KS_func_db_macro(ddc_D0D0) ,
KS_func_db_macro(ddc_D2D2) ,
KS_func_db_macro(ddc_D0D1) ,
KS_func_db_macro(ddc_D0D2) ,
KS_func_db_macro(ddc_D1D2) ,
KS_func_db_macro(ddc_D1D1) ,
KS_func_db_macro(dddc_D0D1D1) ,
KS_func_db_macro(dddc_D1D2D2) ,
KS_func_db_macro(dddc_D2D2D0) ,
KS_func_db_macro(dddc_D0D2D1) ,
KS_func_db_macro(dddc_D2D2D2) ,
KS_func_db_macro(dddc_D1D1D2) ,
KS_func_db_macro(dddc_D0D1D2) ,
KS_func_db_macro(dddc_D1D1D0) ,
KS_func_db_macro(dddc_D1D2D1) ,
KS_func_db_macro(dddc_D0D0D2) ,
KS_func_db_macro(dddc_D0D2D0) ,
KS_func_db_macro(dddc_D1D1D1) ,
KS_func_db_macro(dddc_D1D2D0) ,
KS_func_db_macro(dddc_D0D1D0) ,
KS_func_db_macro(dddc_D0D0D0) ,
KS_func_db_macro(dddc_D0D0D1) ,
KS_func_db_macro(dddc_D0D2D2) ,
KS_func_db_macro(dddc_D2D2D1) ,

KS_func_db_macro(X) ,
KS_func_db_macro(dX_D0) ,
KS_func_db_macro(dX_D2) ,
KS_func_db_macro(dX_D1) ,
KS_func_db_macro(ddX_D1D1) ,
KS_func_db_macro(ddX_D0D2) ,
KS_func_db_macro(ddX_D0D0) ,
KS_func_db_macro(ddX_D1D2) ,
KS_func_db_macro(ddX_D0D1) ,
KS_func_db_macro(ddX_D2D2) ,
KS_func_db_macro(dddX_D2D2D0) ,
KS_func_db_macro(dddX_D0D2D0) ,
KS_func_db_macro(dddX_D0D2D2) ,
KS_func_db_macro(dddX_D0D2D1) ,
KS_func_db_macro(dddX_D1D2D0) ,
KS_func_db_macro(dddX_D1D1D2) ,
KS_func_db_macro(dddX_D1D2D2) ,
KS_func_db_macro(dddX_D2D2D2) ,
KS_func_db_macro(dddX_D0D1D1) ,
KS_func_db_macro(dddX_D0D0D0) ,
KS_func_db_macro(dddX_D0D0D1) ,
KS_func_db_macro(dddX_D1D1D1) ,
KS_func_db_macro(dddX_D2D2D1) ,
KS_func_db_macro(dddX_D1D2D1) ,
KS_func_db_macro(dddX_D0D1D2) ,
KS_func_db_macro(dddX_D0D0D2) ,
KS_func_db_macro(dddX_D0D1D0) ,
KS_func_db_macro(dddX_D1D1D0) ,
KS_func_db_macro(Y) ,
KS_func_db_macro(dY_D0) ,
KS_func_db_macro(dY_D2) ,
KS_func_db_macro(dY_D1) ,
KS_func_db_macro(ddY_D1D2) ,
KS_func_db_macro(ddY_D0D1) ,
KS_func_db_macro(ddY_D1D1) ,
KS_func_db_macro(ddY_D0D2) ,
KS_func_db_macro(ddY_D2D2) ,
KS_func_db_macro(ddY_D0D0) ,
KS_func_db_macro(dddY_D1D1D1) ,
KS_func_db_macro(dddY_D1D1D2) ,
KS_func_db_macro(dddY_D2D2D2) ,
KS_func_db_macro(dddY_D0D1D2) ,
KS_func_db_macro(dddY_D1D2D2) ,
KS_func_db_macro(dddY_D0D0D0) ,
KS_func_db_macro(dddY_D0D1D0) ,
KS_func_db_macro(dddY_D0D2D2) ,
KS_func_db_macro(dddY_D1D2D1) ,
KS_func_db_macro(dddY_D0D2D1) ,
KS_func_db_macro(dddY_D0D2D0) ,
KS_func_db_macro(dddY_D0D0D1) ,
KS_func_db_macro(dddY_D2D2D0) ,
KS_func_db_macro(dddY_D2D2D1) ,
KS_func_db_macro(dddY_D0D0D2) ,
KS_func_db_macro(dddY_D1D2D0) ,
KS_func_db_macro(dddY_D1D1D0) ,
KS_func_db_macro(dddY_D0D1D1) ,
KS_func_db_macro(Z) ,
KS_func_db_macro(dZ_D2) ,
KS_func_db_macro(dZ_D1) ,
KS_func_db_macro(dZ_D0) ,
KS_func_db_macro(ddZ_D2D2) ,
KS_func_db_macro(ddZ_D0D0) ,
KS_func_db_macro(ddZ_D0D1) ,
KS_func_db_macro(ddZ_D1D2) ,
KS_func_db_macro(ddZ_D1D1) ,
KS_func_db_macro(ddZ_D0D2) ,
KS_func_db_macro(dddZ_D0D0D0) ,
KS_func_db_macro(dddZ_D1D2D2) ,
KS_func_db_macro(dddZ_D0D0D2) ,
KS_func_db_macro(dddZ_D0D1D0) ,
KS_func_db_macro(dddZ_D1D1D1) ,
KS_func_db_macro(dddZ_D0D2D2) ,
KS_func_db_macro(dddZ_D1D2D1) ,
KS_func_db_macro(dddZ_D0D2D0) ,
KS_func_db_macro(dddZ_D0D1D2) ,
KS_func_db_macro(dddZ_D0D0D1) ,
KS_func_db_macro(dddZ_D1D1D0) ,
KS_func_db_macro(dddZ_D2D2D0) ,
KS_func_db_macro(dddZ_D1D1D2) ,
KS_func_db_macro(dddZ_D0D2D1) ,
KS_func_db_macro(dddZ_D0D1D1) ,
KS_func_db_macro(dddZ_D2D2D1) ,
KS_func_db_macro(dddZ_D2D2D2) ,
KS_func_db_macro(dddZ_D1D2D0) ,
KS_func_db_macro(R) ,
KS_func_db_macro(dR_D2) ,
KS_func_db_macro(dR_D1) ,
KS_func_db_macro(dR_D0) ,
KS_func_db_macro(ddR_D1D2) ,
KS_func_db_macro(ddR_D1D1) ,
KS_func_db_macro(ddR_D2D2) ,
KS_func_db_macro(ddR_D0D2) ,
KS_func_db_macro(ddR_D0D0) ,
KS_func_db_macro(ddR_D0D1) ,
KS_func_db_macro(dddR_D0D0D2) ,
KS_func_db_macro(dddR_D0D2D2) ,
KS_func_db_macro(dddR_D2D2D2) ,
KS_func_db_macro(dddR_D2D2D0) ,
KS_func_db_macro(dddR_D1D1D0) ,
KS_func_db_macro(dddR_D2D2D1) ,
KS_func_db_macro(dddR_D0D1D1) ,
KS_func_db_macro(dddR_D1D2D2) ,
KS_func_db_macro(dddR_D0D0D0) ,
KS_func_db_macro(dddR_D0D1D0) ,
KS_func_db_macro(dddR_D0D2D0) ,
KS_func_db_macro(dddR_D0D0D1) ,
KS_func_db_macro(dddR_D0D2D1) ,
KS_func_db_macro(dddR_D0D1D2) ,
KS_func_db_macro(dddR_D1D1D2) ,
KS_func_db_macro(dddR_D1D2D1) ,
KS_func_db_macro(dddR_D1D1D1) ,
KS_func_db_macro(dddR_D1D2D0) ,

KS_func_db_macro(rolloff) ,
KS_func_db_macro(drolloff_D0) ,
KS_func_db_macro(drolloff_D2) ,
KS_func_db_macro(drolloff_D1) ,
KS_func_db_macro(ddrolloff_D2D2) ,
KS_func_db_macro(ddrolloff_D1D2) ,
KS_func_db_macro(ddrolloff_D0D1) ,
KS_func_db_macro(ddrolloff_D0D2) ,
KS_func_db_macro(ddrolloff_D0D0) ,
KS_func_db_macro(ddrolloff_D1D1) ,
KS_func_db_macro(dddrolloff_D0D1D0) ,
KS_func_db_macro(dddrolloff_D1D1D2) ,
KS_func_db_macro(dddrolloff_D0D2D0) ,
KS_func_db_macro(dddrolloff_D2D2D1) ,
KS_func_db_macro(dddrolloff_D0D0D0) ,
KS_func_db_macro(dddrolloff_D0D1D2) ,
KS_func_db_macro(dddrolloff_D0D2D2) ,
KS_func_db_macro(dddrolloff_D1D1D0) ,
KS_func_db_macro(dddrolloff_D1D2D1) ,
KS_func_db_macro(dddrolloff_D2D2D2) ,
KS_func_db_macro(dddrolloff_D1D2D0) ,
KS_func_db_macro(dddrolloff_D1D2D2) ,
KS_func_db_macro(dddrolloff_D2D2D0) ,
KS_func_db_macro(dddrolloff_D0D2D1) ,
KS_func_db_macro(dddrolloff_D0D1D1) ,
KS_func_db_macro(dddrolloff_D1D1D1) ,
KS_func_db_macro(dddrolloff_D0D0D2) ,
KS_func_db_macro(dddrolloff_D0D0D1) ,
0
};

/* derivative name data base */
static char *const derive_name_db[] = 
{
KS_name_db_macro(K0),
KS_name_db_macro(dK0_D1),
KS_name_db_macro(dK0_D0),
KS_name_db_macro(dK0_D2),
KS_name_db_macro(ddK0_D1D2),
KS_name_db_macro(ddK0_D0D0),
KS_name_db_macro(ddK0_D2D2),
KS_name_db_macro(ddK0_D1D1),
KS_name_db_macro(ddK0_D0D2),
KS_name_db_macro(ddK0_D0D1),
KS_name_db_macro(dddK0_D2D2D1),
KS_name_db_macro(dddK0_D2D2D0),
KS_name_db_macro(dddK0_D0D2D1),
KS_name_db_macro(dddK0_D0D1D1),
KS_name_db_macro(dddK0_D2D2D2),
KS_name_db_macro(dddK0_D1D2D2),
KS_name_db_macro(dddK0_D0D2D2),
KS_name_db_macro(dddK0_D1D2D1),
KS_name_db_macro(dddK0_D0D0D0),
KS_name_db_macro(dddK0_D0D0D1),
KS_name_db_macro(dddK0_D1D1D0),
KS_name_db_macro(dddK0_D1D2D0),
KS_name_db_macro(dddK0_D0D2D0),
KS_name_db_macro(dddK0_D0D1D0),
KS_name_db_macro(dddK0_D0D1D2),
KS_name_db_macro(dddK0_D0D0D2),
KS_name_db_macro(dddK0_D1D1D1),
KS_name_db_macro(dddK0_D1D1D2),
KS_name_db_macro(K1),
KS_name_db_macro(dK1_D2),
KS_name_db_macro(dK1_D0),
KS_name_db_macro(dK1_D1),
KS_name_db_macro(ddK1_D0D1),
KS_name_db_macro(ddK1_D1D1),
KS_name_db_macro(ddK1_D2D2),
KS_name_db_macro(ddK1_D0D2),
KS_name_db_macro(ddK1_D1D2),
KS_name_db_macro(ddK1_D0D0),
KS_name_db_macro(dddK1_D0D2D1),
KS_name_db_macro(dddK1_D1D1D0),
KS_name_db_macro(dddK1_D2D2D2),
KS_name_db_macro(dddK1_D1D2D0),
KS_name_db_macro(dddK1_D1D1D1),
KS_name_db_macro(dddK1_D0D1D2),
KS_name_db_macro(dddK1_D0D2D2),
KS_name_db_macro(dddK1_D0D0D0),
KS_name_db_macro(dddK1_D2D2D0),
KS_name_db_macro(dddK1_D2D2D1),
KS_name_db_macro(dddK1_D0D2D0),
KS_name_db_macro(dddK1_D0D0D2),
KS_name_db_macro(dddK1_D1D2D1),
KS_name_db_macro(dddK1_D1D1D2),
KS_name_db_macro(dddK1_D0D1D0),
KS_name_db_macro(dddK1_D1D2D2),
KS_name_db_macro(dddK1_D0D1D1),
KS_name_db_macro(dddK1_D0D0D1),
KS_name_db_macro(K2),
KS_name_db_macro(dK2_D0),
KS_name_db_macro(dK2_D1),
KS_name_db_macro(dK2_D2),
KS_name_db_macro(ddK2_D1D2),
KS_name_db_macro(ddK2_D0D1),
KS_name_db_macro(ddK2_D2D2),
KS_name_db_macro(ddK2_D1D1),
KS_name_db_macro(ddK2_D0D0),
KS_name_db_macro(ddK2_D0D2),
KS_name_db_macro(dddK2_D0D2D0),
KS_name_db_macro(dddK2_D1D2D1),
KS_name_db_macro(dddK2_D0D1D0),
KS_name_db_macro(dddK2_D0D2D2),
KS_name_db_macro(dddK2_D0D0D1),
KS_name_db_macro(dddK2_D0D1D1),
KS_name_db_macro(dddK2_D1D2D2),
KS_name_db_macro(dddK2_D1D1D1),
KS_name_db_macro(dddK2_D0D1D2),
KS_name_db_macro(dddK2_D1D1D2),
KS_name_db_macro(dddK2_D1D2D0),
KS_name_db_macro(dddK2_D1D1D0),
KS_name_db_macro(dddK2_D2D2D0),
KS_name_db_macro(dddK2_D0D0D0),
KS_name_db_macro(dddK2_D0D2D1),
KS_name_db_macro(dddK2_D0D0D2),
KS_name_db_macro(dddK2_D2D2D1),
KS_name_db_macro(dddK2_D2D2D2),
KS_name_db_macro(H),
KS_name_db_macro(dH_D1),
KS_name_db_macro(dH_D0),
KS_name_db_macro(dH_D2),
KS_name_db_macro(ddH_D1D1),
KS_name_db_macro(ddH_D0D2),
KS_name_db_macro(ddH_D2D2),
KS_name_db_macro(ddH_D0D1),
KS_name_db_macro(ddH_D0D0),
KS_name_db_macro(ddH_D1D2),
KS_name_db_macro(dddH_D0D2D2),
KS_name_db_macro(dddH_D1D2D0),
KS_name_db_macro(dddH_D0D1D0),
KS_name_db_macro(dddH_D1D1D2),
KS_name_db_macro(dddH_D0D2D0),
KS_name_db_macro(dddH_D1D2D1),
KS_name_db_macro(dddH_D0D0D0),
KS_name_db_macro(dddH_D2D2D2),
KS_name_db_macro(dddH_D1D2D2),
KS_name_db_macro(dddH_D0D1D2),
KS_name_db_macro(dddH_D0D0D1),
KS_name_db_macro(dddH_D1D1D1),
KS_name_db_macro(dddH_D0D2D1),
KS_name_db_macro(dddH_D1D1D0),
KS_name_db_macro(dddH_D0D0D2),
KS_name_db_macro(dddH_D2D2D0),
KS_name_db_macro(dddH_D0D1D1),
KS_name_db_macro(dddH_D2D2D1),

KS_name_db_macro(k0) ,
KS_name_db_macro(dk0_D2) ,
KS_name_db_macro(dk0_D1) ,
KS_name_db_macro(dk0_D0) ,
KS_name_db_macro(ddk0_D1D2) ,
KS_name_db_macro(ddk0_D0D0) ,
KS_name_db_macro(ddk0_D1D1) ,
KS_name_db_macro(ddk0_D0D2) ,
KS_name_db_macro(ddk0_D0D1) ,
KS_name_db_macro(ddk0_D2D2) ,
KS_name_db_macro(dddk0_D2D2D0) ,
KS_name_db_macro(dddk0_D0D2D2) ,
KS_name_db_macro(dddk0_D0D1D1) ,
KS_name_db_macro(dddk0_D0D0D1) ,
KS_name_db_macro(dddk0_D0D2D0) ,
KS_name_db_macro(dddk0_D0D1D2) ,
KS_name_db_macro(dddk0_D1D1D2) ,
KS_name_db_macro(dddk0_D0D2D1) ,
KS_name_db_macro(dddk0_D1D1D0) ,
KS_name_db_macro(dddk0_D1D2D0) ,
KS_name_db_macro(dddk0_D0D1D0) ,
KS_name_db_macro(dddk0_D1D2D1) ,
KS_name_db_macro(dddk0_D0D0D2) ,
KS_name_db_macro(dddk0_D2D2D1) ,
KS_name_db_macro(dddk0_D1D1D1) ,
KS_name_db_macro(dddk0_D2D2D2) ,
KS_name_db_macro(dddk0_D0D0D0) ,
KS_name_db_macro(dddk0_D1D2D2) ,
KS_name_db_macro(k1) ,
KS_name_db_macro(dk1_D2) ,
KS_name_db_macro(dk1_D0) ,
KS_name_db_macro(dk1_D1) ,
KS_name_db_macro(ddk1_D1D1) ,
KS_name_db_macro(ddk1_D1D2) ,
KS_name_db_macro(ddk1_D2D2) ,
KS_name_db_macro(ddk1_D0D2) ,
KS_name_db_macro(ddk1_D0D0) ,
KS_name_db_macro(ddk1_D0D1) ,
KS_name_db_macro(dddk1_D0D1D2) ,
KS_name_db_macro(dddk1_D0D1D0) ,
KS_name_db_macro(dddk1_D0D2D1) ,
KS_name_db_macro(dddk1_D1D2D1) ,
KS_name_db_macro(dddk1_D2D2D2) ,
KS_name_db_macro(dddk1_D1D2D2) ,
KS_name_db_macro(dddk1_D2D2D1) ,
KS_name_db_macro(dddk1_D0D0D2) ,
KS_name_db_macro(dddk1_D1D2D0) ,
KS_name_db_macro(dddk1_D0D2D0) ,
KS_name_db_macro(dddk1_D1D1D1) ,
KS_name_db_macro(dddk1_D0D2D2) ,
KS_name_db_macro(dddk1_D0D1D1) ,
KS_name_db_macro(dddk1_D0D0D0) ,
KS_name_db_macro(dddk1_D0D0D1) ,
KS_name_db_macro(dddk1_D2D2D0) ,
KS_name_db_macro(dddk1_D1D1D0) ,
KS_name_db_macro(dddk1_D1D1D2) ,
KS_name_db_macro(k2) ,
KS_name_db_macro(dk2_D2) ,
KS_name_db_macro(dk2_D1) ,
KS_name_db_macro(dk2_D0) ,
KS_name_db_macro(ddk2_D0D1) ,
KS_name_db_macro(ddk2_D0D0) ,
KS_name_db_macro(ddk2_D0D2) ,
KS_name_db_macro(ddk2_D2D2) ,
KS_name_db_macro(ddk2_D1D1) ,
KS_name_db_macro(ddk2_D1D2) ,
KS_name_db_macro(dddk2_D0D1D2) ,
KS_name_db_macro(dddk2_D1D2D2) ,
KS_name_db_macro(dddk2_D1D1D2) ,
KS_name_db_macro(dddk2_D2D2D2) ,
KS_name_db_macro(dddk2_D1D1D0) ,
KS_name_db_macro(dddk2_D2D2D0) ,
KS_name_db_macro(dddk2_D0D2D1) ,
KS_name_db_macro(dddk2_D1D2D0) ,
KS_name_db_macro(dddk2_D0D0D1) ,
KS_name_db_macro(dddk2_D0D2D2) ,
KS_name_db_macro(dddk2_D2D2D1) ,
KS_name_db_macro(dddk2_D0D1D1) ,
KS_name_db_macro(dddk2_D0D2D0) ,
KS_name_db_macro(dddk2_D1D1D1) ,
KS_name_db_macro(dddk2_D1D2D1) ,
KS_name_db_macro(dddk2_D0D0D2) ,
KS_name_db_macro(dddk2_D0D0D0) ,
KS_name_db_macro(dddk2_D0D1D0) ,
KS_name_db_macro(kt) ,
KS_name_db_macro(dkt_D1) ,
KS_name_db_macro(dkt_D2) ,
KS_name_db_macro(dkt_D0) ,
KS_name_db_macro(ddkt_D2D2) ,
KS_name_db_macro(ddkt_D0D0) ,
KS_name_db_macro(ddkt_D0D1) ,
KS_name_db_macro(ddkt_D1D1) ,
KS_name_db_macro(ddkt_D0D2) ,
KS_name_db_macro(ddkt_D1D2) ,
KS_name_db_macro(dddkt_D0D1D1) ,
KS_name_db_macro(dddkt_D1D2D2) ,
KS_name_db_macro(dddkt_D0D1D2) ,
KS_name_db_macro(dddkt_D0D0D0) ,
KS_name_db_macro(dddkt_D0D2D0) ,
KS_name_db_macro(dddkt_D0D2D1) ,
KS_name_db_macro(dddkt_D0D2D2) ,
KS_name_db_macro(dddkt_D2D2D2) ,
KS_name_db_macro(dddkt_D1D2D0) ,
KS_name_db_macro(dddkt_D0D0D2) ,
KS_name_db_macro(dddkt_D1D1D0) ,
KS_name_db_macro(dddkt_D1D1D2) ,
KS_name_db_macro(dddkt_D2D2D1) ,
KS_name_db_macro(dddkt_D1D1D1) ,
KS_name_db_macro(dddkt_D2D2D0) ,
KS_name_db_macro(dddkt_D1D2D1) ,
KS_name_db_macro(dddkt_D0D0D1) ,
KS_name_db_macro(dddkt_D0D1D0) ,
KS_name_db_macro(c) ,
KS_name_db_macro(dc_D2) ,
KS_name_db_macro(dc_D1) ,
KS_name_db_macro(dc_D0) ,
KS_name_db_macro(ddc_D0D0) ,
KS_name_db_macro(ddc_D2D2) ,
KS_name_db_macro(ddc_D0D1) ,
KS_name_db_macro(ddc_D0D2) ,
KS_name_db_macro(ddc_D1D2) ,
KS_name_db_macro(ddc_D1D1) ,
KS_name_db_macro(dddc_D0D1D1) ,
KS_name_db_macro(dddc_D1D2D2) ,
KS_name_db_macro(dddc_D2D2D0) ,
KS_name_db_macro(dddc_D0D2D1) ,
KS_name_db_macro(dddc_D2D2D2) ,
KS_name_db_macro(dddc_D1D1D2) ,
KS_name_db_macro(dddc_D0D1D2) ,
KS_name_db_macro(dddc_D1D1D0) ,
KS_name_db_macro(dddc_D1D2D1) ,
KS_name_db_macro(dddc_D0D0D2) ,
KS_name_db_macro(dddc_D0D2D0) ,
KS_name_db_macro(dddc_D1D1D1) ,
KS_name_db_macro(dddc_D1D2D0) ,
KS_name_db_macro(dddc_D0D1D0) ,
KS_name_db_macro(dddc_D0D0D0) ,
KS_name_db_macro(dddc_D0D0D1) ,
KS_name_db_macro(dddc_D0D2D2) ,
KS_name_db_macro(dddc_D2D2D1) ,

KS_name_db_macro(X) ,
KS_name_db_macro(dX_D0) ,
KS_name_db_macro(dX_D2) ,
KS_name_db_macro(dX_D1) ,
KS_name_db_macro(ddX_D1D1) ,
KS_name_db_macro(ddX_D0D2) ,
KS_name_db_macro(ddX_D0D0) ,
KS_name_db_macro(ddX_D1D2) ,
KS_name_db_macro(ddX_D0D1) ,
KS_name_db_macro(ddX_D2D2) ,
KS_name_db_macro(dddX_D2D2D0) ,
KS_name_db_macro(dddX_D0D2D0) ,
KS_name_db_macro(dddX_D0D2D2) ,
KS_name_db_macro(dddX_D0D2D1) ,
KS_name_db_macro(dddX_D1D2D0) ,
KS_name_db_macro(dddX_D1D1D2) ,
KS_name_db_macro(dddX_D1D2D2) ,
KS_name_db_macro(dddX_D2D2D2) ,
KS_name_db_macro(dddX_D0D1D1) ,
KS_name_db_macro(dddX_D0D0D0) ,
KS_name_db_macro(dddX_D0D0D1) ,
KS_name_db_macro(dddX_D1D1D1) ,
KS_name_db_macro(dddX_D2D2D1) ,
KS_name_db_macro(dddX_D1D2D1) ,
KS_name_db_macro(dddX_D0D1D2) ,
KS_name_db_macro(dddX_D0D0D2) ,
KS_name_db_macro(dddX_D0D1D0) ,
KS_name_db_macro(dddX_D1D1D0) ,
KS_name_db_macro(Y) ,
KS_name_db_macro(dY_D0) ,
KS_name_db_macro(dY_D2) ,
KS_name_db_macro(dY_D1) ,
KS_name_db_macro(ddY_D1D2) ,
KS_name_db_macro(ddY_D0D1) ,
KS_name_db_macro(ddY_D1D1) ,
KS_name_db_macro(ddY_D0D2) ,
KS_name_db_macro(ddY_D2D2) ,
KS_name_db_macro(ddY_D0D0) ,
KS_name_db_macro(dddY_D1D1D1) ,
KS_name_db_macro(dddY_D1D1D2) ,
KS_name_db_macro(dddY_D2D2D2) ,
KS_name_db_macro(dddY_D0D1D2) ,
KS_name_db_macro(dddY_D1D2D2) ,
KS_name_db_macro(dddY_D0D0D0) ,
KS_name_db_macro(dddY_D0D1D0) ,
KS_name_db_macro(dddY_D0D2D2) ,
KS_name_db_macro(dddY_D1D2D1) ,
KS_name_db_macro(dddY_D0D2D1) ,
KS_name_db_macro(dddY_D0D2D0) ,
KS_name_db_macro(dddY_D0D0D1) ,
KS_name_db_macro(dddY_D2D2D0) ,
KS_name_db_macro(dddY_D2D2D1) ,
KS_name_db_macro(dddY_D0D0D2) ,
KS_name_db_macro(dddY_D1D2D0) ,
KS_name_db_macro(dddY_D1D1D0) ,
KS_name_db_macro(dddY_D0D1D1) ,
KS_name_db_macro(Z) ,
KS_name_db_macro(dZ_D2) ,
KS_name_db_macro(dZ_D1) ,
KS_name_db_macro(dZ_D0) ,
KS_name_db_macro(ddZ_D2D2) ,
KS_name_db_macro(ddZ_D0D0) ,
KS_name_db_macro(ddZ_D0D1) ,
KS_name_db_macro(ddZ_D1D2) ,
KS_name_db_macro(ddZ_D1D1) ,
KS_name_db_macro(ddZ_D0D2) ,
KS_name_db_macro(dddZ_D0D0D0) ,
KS_name_db_macro(dddZ_D1D2D2) ,
KS_name_db_macro(dddZ_D0D0D2) ,
KS_name_db_macro(dddZ_D0D1D0) ,
KS_name_db_macro(dddZ_D1D1D1) ,
KS_name_db_macro(dddZ_D0D2D2) ,
KS_name_db_macro(dddZ_D1D2D1) ,
KS_name_db_macro(dddZ_D0D2D0) ,
KS_name_db_macro(dddZ_D0D1D2) ,
KS_name_db_macro(dddZ_D0D0D1) ,
KS_name_db_macro(dddZ_D1D1D0) ,
KS_name_db_macro(dddZ_D2D2D0) ,
KS_name_db_macro(dddZ_D1D1D2) ,
KS_name_db_macro(dddZ_D0D2D1) ,
KS_name_db_macro(dddZ_D0D1D1) ,
KS_name_db_macro(dddZ_D2D2D1) ,
KS_name_db_macro(dddZ_D2D2D2) ,
KS_name_db_macro(dddZ_D1D2D0) ,
KS_name_db_macro(R) ,
KS_name_db_macro(dR_D2) ,
KS_name_db_macro(dR_D1) ,
KS_name_db_macro(dR_D0) ,
KS_name_db_macro(ddR_D1D2) ,
KS_name_db_macro(ddR_D1D1) ,
KS_name_db_macro(ddR_D2D2) ,
KS_name_db_macro(ddR_D0D2) ,
KS_name_db_macro(ddR_D0D0) ,
KS_name_db_macro(ddR_D0D1) ,
KS_name_db_macro(dddR_D0D0D2) ,
KS_name_db_macro(dddR_D0D2D2) ,
KS_name_db_macro(dddR_D2D2D2) ,
KS_name_db_macro(dddR_D2D2D0) ,
KS_name_db_macro(dddR_D1D1D0) ,
KS_name_db_macro(dddR_D2D2D1) ,
KS_name_db_macro(dddR_D0D1D1) ,
KS_name_db_macro(dddR_D1D2D2) ,
KS_name_db_macro(dddR_D0D0D0) ,
KS_name_db_macro(dddR_D0D1D0) ,
KS_name_db_macro(dddR_D0D2D0) ,
KS_name_db_macro(dddR_D0D0D1) ,
KS_name_db_macro(dddR_D0D2D1) ,
KS_name_db_macro(dddR_D0D1D2) ,
KS_name_db_macro(dddR_D1D1D2) ,
KS_name_db_macro(dddR_D1D2D1) ,
KS_name_db_macro(dddR_D1D1D1) ,
KS_name_db_macro(dddR_D1D2D0) ,

KS_name_db_macro(rolloff) ,
KS_name_db_macro(drolloff_D0) ,
KS_name_db_macro(drolloff_D2) ,
KS_name_db_macro(drolloff_D1) ,
KS_name_db_macro(ddrolloff_D2D2) ,
KS_name_db_macro(ddrolloff_D1D2) ,
KS_name_db_macro(ddrolloff_D0D1) ,
KS_name_db_macro(ddrolloff_D0D2) ,
KS_name_db_macro(ddrolloff_D0D0) ,
KS_name_db_macro(ddrolloff_D1D1) ,
KS_name_db_macro(dddrolloff_D0D1D0) ,
KS_name_db_macro(dddrolloff_D1D1D2) ,
KS_name_db_macro(dddrolloff_D0D2D0) ,
KS_name_db_macro(dddrolloff_D2D2D1) ,
KS_name_db_macro(dddrolloff_D0D0D0) ,
KS_name_db_macro(dddrolloff_D0D1D2) ,
KS_name_db_macro(dddrolloff_D0D2D2) ,
KS_name_db_macro(dddrolloff_D1D1D0) ,
KS_name_db_macro(dddrolloff_D1D2D1) ,
KS_name_db_macro(dddrolloff_D2D2D2) ,
KS_name_db_macro(dddrolloff_D1D2D0) ,
KS_name_db_macro(dddrolloff_D1D2D2) ,
KS_name_db_macro(dddrolloff_D2D2D0) ,
KS_name_db_macro(dddrolloff_D0D2D1) ,
KS_name_db_macro(dddrolloff_D0D1D1) ,
KS_name_db_macro(dddrolloff_D1D1D1) ,
KS_name_db_macro(dddrolloff_D0D0D2) ,
KS_name_db_macro(dddrolloff_D0D0D1) ,
0
};

