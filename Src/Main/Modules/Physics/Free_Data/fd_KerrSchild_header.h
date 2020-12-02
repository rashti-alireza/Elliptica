#include "fd_header.h"

/* for external variables DON'T CHANGE THIS */
#define KS_glob_var(x) (fd_ks_glob##x)

#define M_BH    KS_glob_var(M_BH) /* BH mass */
#define a_BH    KS_glob_var(a_BH) /* BH spin */
#define phiy    KS_glob_var(phiy) /* rotation angel */
#define phiz    KS_glob_var(phiz) /* rotation angel */
#define Bx      KS_glob_var(Bx) /* B=v/c (boost) */
#define By      KS_glob_var(By) /* B=v/c (boost) */
#define Bz      KS_glob_var(Bz) /* B=v/c (boost) */
#define B2      KS_glob_var(B2) /* B^i B_i */
#define Lambda  KS_glob_var(Lambda) /* flat data => 0, kerr-schild => 1 */

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

/* concatenate */
#define Pattern(a,b) a##_##b

/* function prefix DON NOT change this prefix macro some fields 
// using this prefix which then you must change them too.*/
#define KS_prefix           "fd_ks_"
#define KS_func_prefix(name) fd_ks_##name

/* function names and prototype */
#define KS_func_def_macro(name) double KS_func_prefix(name)

/* function args */
#define KS_func_args_macro  (const struct Analytic_Func_Arg_S *const farg __attribute__((unused)))

/* pass special argument to each function */   
#define KS_func_pass_args_macro  (farg)

/* pass arguments for the following functions */
#define fd_ks_rolloff(x,y,z)  (fd_ks_rolloff KS_func_pass_args_macro)
#define fd_ks_H(x,y,z)  (fd_ks_H  KS_func_pass_args_macro)
#define fd_ks_k0(x,y,z) (fd_ks_k0 KS_func_pass_args_macro)
#define fd_ks_k1(x,y,z) (fd_ks_k1 KS_func_pass_args_macro)
#define fd_ks_k2(x,y,z) (fd_ks_k2 KS_func_pass_args_macro)
#define fd_ks_kt(x,y,z) (fd_ks_kt KS_func_pass_args_macro)
#define fd_ks_c(x,y,z)  (fd_ks_c  KS_func_pass_args_macro)

#define fd_ks_K0(x,y,z) (fd_ks_K0 KS_func_pass_args_macro)
#define fd_ks_K1(x,y,z) (fd_ks_K1 KS_func_pass_args_macro)
#define fd_ks_K2(x,y,z) (fd_ks_K2 KS_func_pass_args_macro)

#define fd_ks_X(x,y,z)  (farg->X)
#define fd_ks_Y(x,y,z)  (farg->Y)
#define fd_ks_Z(x,y,z)  (farg->Z)
#define fd_ks_R(x,y,z)  (farg->R)

#define fd_ks_dX_D0_  (farg->dX_D0)
#define fd_ks_dX_D1_  (farg->dX_D1)
#define fd_ks_dX_D2_  (farg->dX_D2)

#define fd_ks_dY_D0_  (farg->dY_D0)
#define fd_ks_dY_D1_  (farg->dY_D1)
#define fd_ks_dY_D2_  (farg->dY_D2)

#define fd_ks_dZ_D0_  (farg->dZ_D0)
#define fd_ks_dZ_D1_  (farg->dZ_D1)
#define fd_ks_dZ_D2_  (farg->dZ_D2)

#define x  (farg->x)
#define y  (farg->y)
#define z  (farg->z)


/* global variables */
extern double M_BH,a_BH;/* mass and spin of BH */
extern double phiy,phiz;/* rotation angels */
extern double Bx,By,Bz;/* B=v/c (boost) */
extern double B2;/* B^i B_i */
extern double Lambda;/* flat data => 0, kerr-schild => 1 */

void fd_kerr_schild_g_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);
        
void fd_kerr_schild_dg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);
	
void fd_kerr_schild_ddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);

void fd_kerr_schild_dddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);

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

