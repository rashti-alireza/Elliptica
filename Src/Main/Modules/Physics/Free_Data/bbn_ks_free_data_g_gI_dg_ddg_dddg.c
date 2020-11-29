#include "frda_header.h"
#include "bbn_ks_free_data_analytic.h"
#undef x
#undef y
#undef z

#define STR_LEN (99)

#define get_v(stem,indices) \
 (sprintf(f_name,"%s_%s",stem,indices)?patch->fields[Ind(f_name)]->v:0)
 
#define KS_func_pass_args_sed KS_func_pass_args_macro

#define KS_set_args \
  struct KS_Arg_S ksa[1];\
  ksa->x = x;\
  ksa->y = y;\
  ksa->z = z;\
  ksa->X = bbn_ks_X KS_func_pass_args_macro;\
  ksa->Y = bbn_ks_Y KS_func_pass_args_macro;\
  ksa->Z = bbn_ks_Z KS_func_pass_args_macro;\
  ksa->R = bbn_ks_R KS_func_pass_args_macro;\
  ksa->dX_D0 = bbn_ks_dX_D0 KS_func_pass_args_sed;\
  ksa->dX_D1 = bbn_ks_dX_D1 KS_func_pass_args_sed;\
  ksa->dX_D2 = bbn_ks_dX_D2 KS_func_pass_args_sed;\
  ksa->dY_D0 = bbn_ks_dY_D0 KS_func_pass_args_sed;\
  ksa->dY_D1 = bbn_ks_dY_D1 KS_func_pass_args_sed;\
  ksa->dY_D2 = bbn_ks_dY_D2 KS_func_pass_args_sed;\
  ksa->dZ_D0 = bbn_ks_dZ_D0 KS_func_pass_args_sed;\
  ksa->dZ_D1 = bbn_ks_dZ_D1 KS_func_pass_args_sed;\
  ksa->dZ_D2 = bbn_ks_dZ_D2 KS_func_pass_args_sed;

void bbn_free_data_g_gI_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);
        
void bbn_free_data_dg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);
	
void bbn_free_data_ddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);

void bbn_free_data_dddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);

void bbn_free_data_g_gI_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem)

{
  char f_name[STR_LEN];
  const unsigned nn = patch->nn;
  unsigned ijk;
    double *const _gamma_D2D2 = get_v(gConf_stem,"D2D2");
    double *const _gamma_D0D2 = get_v(gConf_stem,"D0D2");
    double *const _gamma_D0D0 = get_v(gConf_stem,"D0D0");
    double *const _gamma_D0D1 = get_v(gConf_stem,"D0D1");
    double *const _gamma_D1D2 = get_v(gConf_stem,"D1D2");
    double *const _gamma_D1D1 = get_v(gConf_stem,"D1D1");
    double *const _gammaI_U0U2 = get_v(igConf_stem,"U0U2");
    double *const _gammaI_U0U0 = get_v(igConf_stem,"U0U0");
    double *const _gammaI_U0U1 = get_v(igConf_stem,"U0U1");
    double *const _gammaI_U1U2 = get_v(igConf_stem,"U1U2");
    double *const _gammaI_U1U1 = get_v(igConf_stem,"U1U1");
    double *const _gammaI_U2U2 = get_v(igConf_stem,"U2U2");
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x,y,z;
      x = patch->node[ijk]->x[0]-BH_center_x;
      y = patch->node[ijk]->x[1]-BH_center_y;
      z = patch->node[ijk]->x[2]-BH_center_z;
      KS_set_args
      _gamma_D0D0[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
bbn_ks_c(x, y, z)*pow(bbn_ks_k0(x, y, z), 2) + 1.0;
      _gamma_D0D1[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z);
      _gamma_D0D2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z);
      _gamma_D1D1[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k1
bbn_ks_c(x, y, z)*pow(bbn_ks_k1(x, y, z), 2) + 1.0;
      _gamma_D1D2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z);
      _gamma_D2D2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k2
bbn_ks_c(x, y, z)*pow(bbn_ks_k2(x, y, z), 2) + 1.0;
      _gammaI_U0U0[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
((pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*bbn_ks_c(x, y, z) + 1)/
((pow(bbn_ks_k0(x, y, z), 2) + pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*
bbn_ks_c(x, y, z) + 1);
      _gammaI_U0U1[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
-bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)/((pow(bbn_ks_k0(x, y, z), 2) + 
pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*bbn_ks_c(x, y, z) + 
1);
      _gammaI_U0U2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
-bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)/((pow(bbn_ks_k0(x, y, z), 2) + 
pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*bbn_ks_c(x, y, z) + 
1);
      _gammaI_U1U1[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
((pow(bbn_ks_k0(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*bbn_ks_c(x, y, z) + 1)/
((pow(bbn_ks_k0(x, y, z), 2) + pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*
bbn_ks_c(x, y, z) + 1);
      _gammaI_U1U2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
-bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)/((pow(bbn_ks_k0(x, y, z), 2) + 
pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*bbn_ks_c(x, y, z) + 
1);
      _gammaI_U2U2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
// k2
((pow(bbn_ks_k0(x, y, z), 2) + pow(bbn_ks_k1(x, y, z), 2))*bbn_ks_c(x, y, z) + 1)/
((pow(bbn_ks_k0(x, y, z), 2) + pow(bbn_ks_k1(x, y, z), 2) + pow(bbn_ks_k2(x, y, z), 2))*
bbn_ks_c(x, y, z) + 1);
      /* quick test check _gamma * _gammaI = delta */
      if (0)
      {
          double delta_U0D0 = 
        _gammaI_U0U0[ijk]*_gamma_D0D0[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D0D2[ijk];
          double delta_U0D1 = 
        _gammaI_U0U0[ijk]*_gamma_D0D1[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D1D2[ijk];
          double delta_U0D2 = 
        _gammaI_U0U0[ijk]*_gamma_D0D2[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U0U2[ijk]*_gamma_D2D2[ijk];
          double delta_U1D2 = 
        _gammaI_U0U1[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U1U2[ijk]*_gamma_D2D2[ijk];
          double delta_U1D0 = 
        _gammaI_U0U1[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D0D2[ijk];
         double delta_U1D1 = 
        _gammaI_U0U1[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D1D2[ijk];
        
          double delta_U2D2 = 
        _gammaI_U0U2[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U2U2[ijk]*_gamma_D2D2[ijk];
          double delta_U2D0 = 
        _gammaI_U0U2[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D0D2[ijk];
          double delta_U2D1 = 
        _gammaI_U0U2[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D1D2[ijk];
        if(!EQL(delta_U1D1,1)||!isfinite(delta_U1D1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D1,0)||!isfinite(delta_U0D1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D2,0)||!isfinite(delta_U0D2))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U1D2,0)||!isfinite(delta_U1D2))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D0,1)||!isfinite(delta_U0D0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D1,0)||!isfinite(delta_U2D1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D2,1)||!isfinite(delta_U2D2))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D0,0)||!isfinite(delta_U2D0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U1D0,0)||!isfinite(delta_U1D0))  Error0("_gammaI is not correct!\n");
      }
    }
}
void bbn_free_data_dg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem)

{

  char f_name[STR_LEN];
  const unsigned nn = patch->nn;
  unsigned ijk;
  
double *const _dgamma_D1D2D2 = get_v(dgConf_stem,"D1D2D2");
double *const _dgamma_D0D0D1 = get_v(dgConf_stem,"D0D0D1");
double *const _dgamma_D0D2D1 = get_v(dgConf_stem,"D0D2D1");
double *const _dgamma_D0D1D0 = get_v(dgConf_stem,"D0D1D0");
double *const _dgamma_D1D2D1 = get_v(dgConf_stem,"D1D2D1");
double *const _dgamma_D2D2D0 = get_v(dgConf_stem,"D2D2D0");
double *const _dgamma_D0D0D0 = get_v(dgConf_stem,"D0D0D0");
double *const _dgamma_D0D0D2 = get_v(dgConf_stem,"D0D0D2");
double *const _dgamma_D0D2D2 = get_v(dgConf_stem,"D0D2D2");
double *const _dgamma_D2D2D1 = get_v(dgConf_stem,"D2D2D1");
double *const _dgamma_D0D1D1 = get_v(dgConf_stem,"D0D1D1");
double *const _dgamma_D0D2D0 = get_v(dgConf_stem,"D0D2D0");
double *const _dgamma_D1D2D0 = get_v(dgConf_stem,"D1D2D0");
double *const _dgamma_D1D1D1 = get_v(dgConf_stem,"D1D1D1");
double *const _dgamma_D0D1D2 = get_v(dgConf_stem,"D0D1D2");
double *const _dgamma_D1D1D0 = get_v(dgConf_stem,"D1D1D0");
double *const _dgamma_D1D1D2 = get_v(dgConf_stem,"D1D1D2");
double *const _dgamma_D2D2D2 = get_v(dgConf_stem,"D2D2D2");

  for (ijk = 0; ijk < nn; ++ijk) 
  {
    double x,y,z;
    x = patch->node[ijk]->x[0]-BH_center_x;
    y = patch->node[ijk]->x[1]-BH_center_y;
    z = patch->node[ijk]->x[2]-BH_center_z;
    KS_set_args

    _dgamma_D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro );

    _dgamma_D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k0
(2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro ))*bbn_ks_k0(x, y, z);

    _dgamma_D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro );

    _dgamma_D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro );

    _dgamma_D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro );

    _dgamma_D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k2
(2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro ))*bbn_ks_k2(x, y, z);

    _dgamma_D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k0
(2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro ))*bbn_ks_k0(x, y, z);

    _dgamma_D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k0
(2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro ))*bbn_ks_k0(x, y, z);

    _dgamma_D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro );

    _dgamma_D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k2
(2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro ))*bbn_ks_k2(x, y, z);

    _dgamma_D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro );

    _dgamma_D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro );

    _dgamma_D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro );

    _dgamma_D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k1
(2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro ))*bbn_ks_k1(x, y, z);

    _dgamma_D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro );

    _dgamma_D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k1
(2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro ))*bbn_ks_k1(x, y, z);

    _dgamma_D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k1
(2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro ))*bbn_ks_k1(x, y, z);

    _dgamma_D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k2
(2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro ))*bbn_ks_k2(x, y, z);

  }
}

void bbn_free_data_ddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem)

{

  char f_name[STR_LEN];
  const unsigned nn = patch->nn;
  unsigned ijk;

double *const _ddgamma_D2D2D0D2 = get_v(ddgConf_stem,"D2D2D0D2");
double *const _ddgamma_D0D1D0D1 = get_v(ddgConf_stem,"D0D1D0D1");
double *const _ddgamma_D1D1D0D1 = get_v(ddgConf_stem,"D1D1D0D1");
double *const _ddgamma_D1D2D0D1 = get_v(ddgConf_stem,"D1D2D0D1");
double *const _ddgamma_D0D0D1D0 = get_v(ddgConf_stem,"D0D0D1D0");
double *const _ddgamma_D1D1D2D0 = get_v(ddgConf_stem,"D1D1D2D0");
double *const _ddgamma_D0D0D0D1 = get_v(ddgConf_stem,"D0D0D0D1");
double *const _ddgamma_D1D1D0D0 = get_v(ddgConf_stem,"D1D1D0D0");
double *const _ddgamma_D0D0D2D2 = get_v(ddgConf_stem,"D0D0D2D2");
double *const _ddgamma_D2D2D1D1 = get_v(ddgConf_stem,"D2D2D1D1");
double *const _ddgamma_D2D2D2D2 = get_v(ddgConf_stem,"D2D2D2D2");
double *const _ddgamma_D0D0D2D0 = get_v(ddgConf_stem,"D0D0D2D0");
double *const _ddgamma_D0D2D2D1 = get_v(ddgConf_stem,"D0D2D2D1");
double *const _ddgamma_D1D2D0D2 = get_v(ddgConf_stem,"D1D2D0D2");
double *const _ddgamma_D1D2D1D0 = get_v(ddgConf_stem,"D1D2D1D0");
double *const _ddgamma_D2D2D2D0 = get_v(ddgConf_stem,"D2D2D2D0");
double *const _ddgamma_D0D0D0D0 = get_v(ddgConf_stem,"D0D0D0D0");
double *const _ddgamma_D0D0D1D1 = get_v(ddgConf_stem,"D0D0D1D1");
double *const _ddgamma_D0D1D0D2 = get_v(ddgConf_stem,"D0D1D0D2");
double *const _ddgamma_D1D2D1D1 = get_v(ddgConf_stem,"D1D2D1D1");
double *const _ddgamma_D0D2D1D1 = get_v(ddgConf_stem,"D0D2D1D1");
double *const _ddgamma_D1D2D2D2 = get_v(ddgConf_stem,"D1D2D2D2");
double *const _ddgamma_D2D2D1D2 = get_v(ddgConf_stem,"D2D2D1D2");
double *const _ddgamma_D0D1D2D0 = get_v(ddgConf_stem,"D0D1D2D0");
double *const _ddgamma_D0D2D2D0 = get_v(ddgConf_stem,"D0D2D2D0");
double *const _ddgamma_D0D1D1D0 = get_v(ddgConf_stem,"D0D1D1D0");
double *const _ddgamma_D1D2D2D0 = get_v(ddgConf_stem,"D1D2D2D0");
double *const _ddgamma_D1D1D1D0 = get_v(ddgConf_stem,"D1D1D1D0");
double *const _ddgamma_D1D2D1D2 = get_v(ddgConf_stem,"D1D2D1D2");
double *const _ddgamma_D1D2D0D0 = get_v(ddgConf_stem,"D1D2D0D0");
double *const _ddgamma_D0D0D0D2 = get_v(ddgConf_stem,"D0D0D0D2");
double *const _ddgamma_D2D2D2D1 = get_v(ddgConf_stem,"D2D2D2D1");
double *const _ddgamma_D0D2D0D0 = get_v(ddgConf_stem,"D0D2D0D0");
double *const _ddgamma_D0D2D0D1 = get_v(ddgConf_stem,"D0D2D0D1");
double *const _ddgamma_D1D1D1D1 = get_v(ddgConf_stem,"D1D1D1D1");
double *const _ddgamma_D0D1D0D0 = get_v(ddgConf_stem,"D0D1D0D0");
double *const _ddgamma_D1D1D2D2 = get_v(ddgConf_stem,"D1D1D2D2");
double *const _ddgamma_D0D2D1D2 = get_v(ddgConf_stem,"D0D2D1D2");
double *const _ddgamma_D0D1D2D1 = get_v(ddgConf_stem,"D0D1D2D1");
double *const _ddgamma_D0D1D1D2 = get_v(ddgConf_stem,"D0D1D1D2");
double *const _ddgamma_D1D1D2D1 = get_v(ddgConf_stem,"D1D1D2D1");
double *const _ddgamma_D0D1D1D1 = get_v(ddgConf_stem,"D0D1D1D1");
double *const _ddgamma_D0D2D1D0 = get_v(ddgConf_stem,"D0D2D1D0");
double *const _ddgamma_D2D2D1D0 = get_v(ddgConf_stem,"D2D2D1D0");
double *const _ddgamma_D0D1D2D2 = get_v(ddgConf_stem,"D0D1D2D2");
double *const _ddgamma_D2D2D0D0 = get_v(ddgConf_stem,"D2D2D0D0");
double *const _ddgamma_D0D0D2D1 = get_v(ddgConf_stem,"D0D0D2D1");
double *const _ddgamma_D0D2D2D2 = get_v(ddgConf_stem,"D0D2D2D2");
double *const _ddgamma_D1D1D1D2 = get_v(ddgConf_stem,"D1D1D1D2");
double *const _ddgamma_D1D2D2D1 = get_v(ddgConf_stem,"D1D2D2D1");
double *const _ddgamma_D0D0D1D2 = get_v(ddgConf_stem,"D0D0D1D2");
double *const _ddgamma_D0D2D0D2 = get_v(ddgConf_stem,"D0D2D0D2");
double *const _ddgamma_D1D1D0D2 = get_v(ddgConf_stem,"D1D1D0D2");
double *const _ddgamma_D2D2D0D1 = get_v(ddgConf_stem,"D2D2D0D1");

  for (ijk = 0; ijk < nn; ++ijk) 
  {
    double x,y,z;
    x = patch->node[ijk]->x[0]-BH_center_x;
    y = patch->node[ijk]->x[1]-BH_center_y;
    z = patch->node[ijk]->x[2]-BH_center_z;
    KS_set_args

    _ddgamma_D2D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _ddgamma_D0D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D1D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D0D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D0D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*pow((bbn_ks_dk1_D0 KS_func_pass_args_macro ), 2) + pow(bbn_ks_k1(x, y, z), 2)*
(bbn_ks_ddc_D0D0 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D0D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*pow((bbn_ks_dk0_D2 KS_func_pass_args_macro ), 2) + pow(bbn_ks_k0(x, y, z), 2)*
(bbn_ks_ddc_D2D2 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro );

    _ddgamma_D2D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*pow((bbn_ks_dk2_D1 KS_func_pass_args_macro ), 2) + pow(bbn_ks_k2(x, y, z), 2)*
(bbn_ks_ddc_D1D1 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _ddgamma_D2D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*pow((bbn_ks_dk2_D2 KS_func_pass_args_macro ), 2) + pow(bbn_ks_k2(x, y, z), 2)*
(bbn_ks_ddc_D2D2 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro );

    _ddgamma_D0D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D0D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _ddgamma_D1D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D1D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D2D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _ddgamma_D0D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*pow((bbn_ks_dk0_D0 KS_func_pass_args_macro ), 2) + pow(bbn_ks_k0(x, y, z), 2)*
(bbn_ks_ddc_D0D0 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D0D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*pow((bbn_ks_dk0_D1 KS_func_pass_args_macro ), 2) + pow(bbn_ks_k0(x, y, z), 2)*
(bbn_ks_ddc_D1D1 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _ddgamma_D0D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _ddgamma_D0D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _ddgamma_D1D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro );

    _ddgamma_D2D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _ddgamma_D0D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D0D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D0D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D1D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D1D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _ddgamma_D1D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D0D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D2D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _ddgamma_D0D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D0D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*pow((bbn_ks_dk1_D1 KS_func_pass_args_macro ), 2) + pow(bbn_ks_k1(x, y, z), 2)*
(bbn_ks_ddc_D1D1 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _ddgamma_D0D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*pow((bbn_ks_dk1_D2 KS_func_pass_args_macro ), 2) + pow(bbn_ks_k1(x, y, z), 2)*
(bbn_ks_ddc_D2D2 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro );

    _ddgamma_D0D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _ddgamma_D0D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _ddgamma_D0D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _ddgamma_D1D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _ddgamma_D0D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _ddgamma_D0D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D2D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _ddgamma_D0D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro );

    _ddgamma_D2D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*pow((bbn_ks_dk2_D0 KS_func_pass_args_macro ), 2) + pow(bbn_ks_k2(x, y, z), 2)*
(bbn_ks_ddc_D0D0 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _ddgamma_D0D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _ddgamma_D0D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro );

    _ddgamma_D1D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _ddgamma_D1D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _ddgamma_D0D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _ddgamma_D0D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D2D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );
  }	
}

void bbn_free_data_dddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem)

{

  char f_name[STR_LEN];
  const unsigned nn = patch->nn;
  unsigned ijk;

double *const _dddgamma_D0D1D2D2D1 = get_v(dddgConf_stem,"D0D1D2D2D1");
double *const _dddgamma_D0D1D2D0D1 = get_v(dddgConf_stem,"D0D1D2D0D1");
double *const _dddgamma_D1D2D0D2D0 = get_v(dddgConf_stem,"D1D2D0D2D0");
double *const _dddgamma_D1D1D1D2D1 = get_v(dddgConf_stem,"D1D1D1D2D1");
double *const _dddgamma_D0D2D0D0D1 = get_v(dddgConf_stem,"D0D2D0D0D1");
double *const _dddgamma_D0D2D0D2D2 = get_v(dddgConf_stem,"D0D2D0D2D2");
double *const _dddgamma_D2D2D2D1D2 = get_v(dddgConf_stem,"D2D2D2D1D2");
double *const _dddgamma_D0D1D0D2D0 = get_v(dddgConf_stem,"D0D1D0D2D0");
double *const _dddgamma_D1D1D2D2D0 = get_v(dddgConf_stem,"D1D1D2D2D0");
double *const _dddgamma_D0D2D0D2D1 = get_v(dddgConf_stem,"D0D2D0D2D1");
double *const _dddgamma_D1D1D2D1D0 = get_v(dddgConf_stem,"D1D1D2D1D0");
double *const _dddgamma_D1D1D0D1D2 = get_v(dddgConf_stem,"D1D1D0D1D2");
double *const _dddgamma_D0D2D1D1D0 = get_v(dddgConf_stem,"D0D2D1D1D0");
double *const _dddgamma_D1D1D1D0D1 = get_v(dddgConf_stem,"D1D1D1D0D1");
double *const _dddgamma_D0D0D2D1D1 = get_v(dddgConf_stem,"D0D0D2D1D1");
double *const _dddgamma_D0D1D2D1D2 = get_v(dddgConf_stem,"D0D1D2D1D2");
double *const _dddgamma_D0D0D0D1D1 = get_v(dddgConf_stem,"D0D0D0D1D1");
double *const _dddgamma_D2D2D1D2D1 = get_v(dddgConf_stem,"D2D2D1D2D1");
double *const _dddgamma_D1D1D1D2D2 = get_v(dddgConf_stem,"D1D1D1D2D2");
double *const _dddgamma_D2D2D2D2D1 = get_v(dddgConf_stem,"D2D2D2D2D1");
double *const _dddgamma_D2D2D1D0D0 = get_v(dddgConf_stem,"D2D2D1D0D0");
double *const _dddgamma_D0D1D2D1D1 = get_v(dddgConf_stem,"D0D1D2D1D1");
double *const _dddgamma_D1D1D2D2D2 = get_v(dddgConf_stem,"D1D1D2D2D2");
double *const _dddgamma_D1D2D1D0D1 = get_v(dddgConf_stem,"D1D2D1D0D1");
double *const _dddgamma_D0D1D1D1D0 = get_v(dddgConf_stem,"D0D1D1D1D0");
double *const _dddgamma_D2D2D2D0D1 = get_v(dddgConf_stem,"D2D2D2D0D1");
double *const _dddgamma_D1D2D2D2D2 = get_v(dddgConf_stem,"D1D2D2D2D2");
double *const _dddgamma_D1D2D1D2D1 = get_v(dddgConf_stem,"D1D2D1D2D1");
double *const _dddgamma_D2D2D0D2D1 = get_v(dddgConf_stem,"D2D2D0D2D1");
double *const _dddgamma_D1D1D2D1D2 = get_v(dddgConf_stem,"D1D1D2D1D2");
double *const _dddgamma_D1D2D2D1D0 = get_v(dddgConf_stem,"D1D2D2D1D0");
double *const _dddgamma_D1D2D0D0D1 = get_v(dddgConf_stem,"D1D2D0D0D1");
double *const _dddgamma_D0D2D0D1D1 = get_v(dddgConf_stem,"D0D2D0D1D1");
double *const _dddgamma_D0D2D1D0D1 = get_v(dddgConf_stem,"D0D2D1D0D1");
double *const _dddgamma_D0D0D0D1D2 = get_v(dddgConf_stem,"D0D0D0D1D2");
double *const _dddgamma_D0D0D2D1D2 = get_v(dddgConf_stem,"D0D0D2D1D2");
double *const _dddgamma_D0D1D0D1D0 = get_v(dddgConf_stem,"D0D1D0D1D0");
double *const _dddgamma_D1D2D1D0D0 = get_v(dddgConf_stem,"D1D2D1D0D0");
double *const _dddgamma_D0D2D1D0D0 = get_v(dddgConf_stem,"D0D2D1D0D0");
double *const _dddgamma_D1D2D0D1D1 = get_v(dddgConf_stem,"D1D2D0D1D1");
double *const _dddgamma_D0D1D2D0D2 = get_v(dddgConf_stem,"D0D1D2D0D2");
double *const _dddgamma_D1D1D2D0D0 = get_v(dddgConf_stem,"D1D1D2D0D0");
double *const _dddgamma_D0D1D2D2D2 = get_v(dddgConf_stem,"D0D1D2D2D2");
double *const _dddgamma_D0D0D2D0D2 = get_v(dddgConf_stem,"D0D0D2D0D2");
double *const _dddgamma_D0D0D0D0D2 = get_v(dddgConf_stem,"D0D0D0D0D2");
double *const _dddgamma_D1D1D0D1D1 = get_v(dddgConf_stem,"D1D1D0D1D1");
double *const _dddgamma_D2D2D0D0D1 = get_v(dddgConf_stem,"D2D2D0D0D1");
double *const _dddgamma_D0D0D1D2D2 = get_v(dddgConf_stem,"D0D0D1D2D2");
double *const _dddgamma_D1D1D0D2D0 = get_v(dddgConf_stem,"D1D1D0D2D0");
double *const _dddgamma_D0D1D1D2D0 = get_v(dddgConf_stem,"D0D1D1D2D0");
double *const _dddgamma_D1D1D0D0D1 = get_v(dddgConf_stem,"D1D1D0D0D1");
double *const _dddgamma_D1D1D0D0D2 = get_v(dddgConf_stem,"D1D1D0D0D2");
double *const _dddgamma_D1D2D0D1D2 = get_v(dddgConf_stem,"D1D2D0D1D2");
double *const _dddgamma_D0D0D1D0D2 = get_v(dddgConf_stem,"D0D0D1D0D2");
double *const _dddgamma_D0D0D0D2D0 = get_v(dddgConf_stem,"D0D0D0D2D0");
double *const _dddgamma_D2D2D2D1D1 = get_v(dddgConf_stem,"D2D2D2D1D1");
double *const _dddgamma_D0D2D0D0D2 = get_v(dddgConf_stem,"D0D2D0D0D2");
double *const _dddgamma_D2D2D2D2D0 = get_v(dddgConf_stem,"D2D2D2D2D0");
double *const _dddgamma_D0D2D0D0D0 = get_v(dddgConf_stem,"D0D2D0D0D0");
double *const _dddgamma_D2D2D1D2D0 = get_v(dddgConf_stem,"D2D2D1D2D0");
double *const _dddgamma_D1D1D0D1D0 = get_v(dddgConf_stem,"D1D1D0D1D0");
double *const _dddgamma_D1D2D2D0D2 = get_v(dddgConf_stem,"D1D2D2D0D2");
double *const _dddgamma_D0D1D0D1D1 = get_v(dddgConf_stem,"D0D1D0D1D1");
double *const _dddgamma_D2D2D0D2D2 = get_v(dddgConf_stem,"D2D2D0D2D2");
double *const _dddgamma_D1D2D2D1D2 = get_v(dddgConf_stem,"D1D2D2D1D2");
double *const _dddgamma_D2D2D1D1D1 = get_v(dddgConf_stem,"D2D2D1D1D1");
double *const _dddgamma_D0D0D0D0D0 = get_v(dddgConf_stem,"D0D0D0D0D0");
double *const _dddgamma_D0D0D1D2D1 = get_v(dddgConf_stem,"D0D0D1D2D1");
double *const _dddgamma_D1D1D1D1D0 = get_v(dddgConf_stem,"D1D1D1D1D0");
double *const _dddgamma_D1D1D1D1D1 = get_v(dddgConf_stem,"D1D1D1D1D1");
double *const _dddgamma_D1D1D1D2D0 = get_v(dddgConf_stem,"D1D1D1D2D0");
double *const _dddgamma_D0D1D0D2D1 = get_v(dddgConf_stem,"D0D1D0D2D1");
double *const _dddgamma_D2D2D1D1D0 = get_v(dddgConf_stem,"D2D2D1D1D0");
double *const _dddgamma_D1D2D2D2D0 = get_v(dddgConf_stem,"D1D2D2D2D0");
double *const _dddgamma_D0D0D1D1D2 = get_v(dddgConf_stem,"D0D0D1D1D2");
double *const _dddgamma_D0D2D0D2D0 = get_v(dddgConf_stem,"D0D2D0D2D0");
double *const _dddgamma_D0D0D1D0D0 = get_v(dddgConf_stem,"D0D0D1D0D0");
double *const _dddgamma_D0D2D2D1D0 = get_v(dddgConf_stem,"D0D2D2D1D0");
double *const _dddgamma_D0D2D2D0D0 = get_v(dddgConf_stem,"D0D2D2D0D0");
double *const _dddgamma_D1D1D1D1D2 = get_v(dddgConf_stem,"D1D1D1D1D2");
double *const _dddgamma_D0D0D2D2D1 = get_v(dddgConf_stem,"D0D0D2D2D1");
double *const _dddgamma_D2D2D0D1D0 = get_v(dddgConf_stem,"D2D2D0D1D0");
double *const _dddgamma_D1D1D2D0D2 = get_v(dddgConf_stem,"D1D1D2D0D2");
double *const _dddgamma_D0D1D1D1D1 = get_v(dddgConf_stem,"D0D1D1D1D1");
double *const _dddgamma_D1D2D0D0D0 = get_v(dddgConf_stem,"D1D2D0D0D0");
double *const _dddgamma_D1D2D0D0D2 = get_v(dddgConf_stem,"D1D2D0D0D2");
double *const _dddgamma_D2D2D0D1D2 = get_v(dddgConf_stem,"D2D2D0D1D2");
double *const _dddgamma_D0D0D1D1D1 = get_v(dddgConf_stem,"D0D0D1D1D1");
double *const _dddgamma_D0D2D0D1D0 = get_v(dddgConf_stem,"D0D2D0D1D0");
double *const _dddgamma_D0D2D1D2D1 = get_v(dddgConf_stem,"D0D2D1D2D1");
double *const _dddgamma_D1D2D1D2D2 = get_v(dddgConf_stem,"D1D2D1D2D2");
double *const _dddgamma_D0D1D1D0D0 = get_v(dddgConf_stem,"D0D1D1D0D0");
double *const _dddgamma_D0D0D1D2D0 = get_v(dddgConf_stem,"D0D0D1D2D0");
double *const _dddgamma_D1D2D1D1D0 = get_v(dddgConf_stem,"D1D2D1D1D0");
double *const _dddgamma_D1D2D2D2D1 = get_v(dddgConf_stem,"D1D2D2D2D1");
double *const _dddgamma_D0D1D0D0D0 = get_v(dddgConf_stem,"D0D1D0D0D0");
double *const _dddgamma_D1D2D2D1D1 = get_v(dddgConf_stem,"D1D2D2D1D1");
double *const _dddgamma_D0D0D0D1D0 = get_v(dddgConf_stem,"D0D0D0D1D0");
double *const _dddgamma_D0D0D0D0D1 = get_v(dddgConf_stem,"D0D0D0D0D1");
double *const _dddgamma_D0D0D0D2D2 = get_v(dddgConf_stem,"D0D0D0D2D2");
double *const _dddgamma_D1D2D0D2D1 = get_v(dddgConf_stem,"D1D2D0D2D1");
double *const _dddgamma_D2D2D1D1D2 = get_v(dddgConf_stem,"D2D2D1D1D2");
double *const _dddgamma_D1D2D1D1D1 = get_v(dddgConf_stem,"D1D2D1D1D1");
double *const _dddgamma_D0D0D1D0D1 = get_v(dddgConf_stem,"D0D0D1D0D1");
double *const _dddgamma_D0D1D0D1D2 = get_v(dddgConf_stem,"D0D1D0D1D2");
double *const _dddgamma_D0D2D1D2D2 = get_v(dddgConf_stem,"D0D2D1D2D2");
double *const _dddgamma_D1D1D0D2D2 = get_v(dddgConf_stem,"D1D1D0D2D2");
double *const _dddgamma_D0D1D0D0D2 = get_v(dddgConf_stem,"D0D1D0D0D2");
double *const _dddgamma_D0D1D0D2D2 = get_v(dddgConf_stem,"D0D1D0D2D2");
double *const _dddgamma_D0D2D1D2D0 = get_v(dddgConf_stem,"D0D2D1D2D0");
double *const _dddgamma_D0D1D0D0D1 = get_v(dddgConf_stem,"D0D1D0D0D1");
double *const _dddgamma_D1D2D1D1D2 = get_v(dddgConf_stem,"D1D2D1D1D2");
double *const _dddgamma_D2D2D0D1D1 = get_v(dddgConf_stem,"D2D2D0D1D1");
double *const _dddgamma_D1D1D2D2D1 = get_v(dddgConf_stem,"D1D1D2D2D1");
double *const _dddgamma_D0D1D2D1D0 = get_v(dddgConf_stem,"D0D1D2D1D0");
double *const _dddgamma_D0D2D1D0D2 = get_v(dddgConf_stem,"D0D2D1D0D2");
double *const _dddgamma_D2D2D2D1D0 = get_v(dddgConf_stem,"D2D2D2D1D0");
double *const _dddgamma_D0D1D1D1D2 = get_v(dddgConf_stem,"D0D1D1D1D2");
double *const _dddgamma_D1D2D1D0D2 = get_v(dddgConf_stem,"D1D2D1D0D2");
double *const _dddgamma_D0D0D2D2D0 = get_v(dddgConf_stem,"D0D0D2D2D0");
double *const _dddgamma_D0D1D1D0D1 = get_v(dddgConf_stem,"D0D1D1D0D1");
double *const _dddgamma_D0D1D2D2D0 = get_v(dddgConf_stem,"D0D1D2D2D0");
double *const _dddgamma_D0D2D0D1D2 = get_v(dddgConf_stem,"D0D2D0D1D2");
double *const _dddgamma_D0D2D2D2D1 = get_v(dddgConf_stem,"D0D2D2D2D1");
double *const _dddgamma_D1D1D0D2D1 = get_v(dddgConf_stem,"D1D1D0D2D1");
double *const _dddgamma_D0D0D2D1D0 = get_v(dddgConf_stem,"D0D0D2D1D0");
double *const _dddgamma_D1D2D0D1D0 = get_v(dddgConf_stem,"D1D2D0D1D0");
double *const _dddgamma_D0D2D2D1D2 = get_v(dddgConf_stem,"D0D2D2D1D2");
double *const _dddgamma_D0D1D1D2D1 = get_v(dddgConf_stem,"D0D1D1D2D1");
double *const _dddgamma_D1D2D2D0D0 = get_v(dddgConf_stem,"D1D2D2D0D0");
double *const _dddgamma_D0D1D2D0D0 = get_v(dddgConf_stem,"D0D1D2D0D0");
double *const _dddgamma_D1D2D1D2D0 = get_v(dddgConf_stem,"D1D2D1D2D0");
double *const _dddgamma_D2D2D1D0D1 = get_v(dddgConf_stem,"D2D2D1D0D1");
double *const _dddgamma_D2D2D1D2D2 = get_v(dddgConf_stem,"D2D2D1D2D2");
double *const _dddgamma_D1D1D1D0D0 = get_v(dddgConf_stem,"D1D1D1D0D0");
double *const _dddgamma_D0D1D1D0D2 = get_v(dddgConf_stem,"D0D1D1D0D2");
double *const _dddgamma_D0D2D2D2D0 = get_v(dddgConf_stem,"D0D2D2D2D0");
double *const _dddgamma_D2D2D2D2D2 = get_v(dddgConf_stem,"D2D2D2D2D2");
double *const _dddgamma_D0D2D2D0D1 = get_v(dddgConf_stem,"D0D2D2D0D1");
double *const _dddgamma_D0D2D1D1D2 = get_v(dddgConf_stem,"D0D2D1D1D2");
double *const _dddgamma_D2D2D0D0D2 = get_v(dddgConf_stem,"D2D2D0D0D2");
double *const _dddgamma_D2D2D2D0D0 = get_v(dddgConf_stem,"D2D2D2D0D0");
double *const _dddgamma_D2D2D0D2D0 = get_v(dddgConf_stem,"D2D2D0D2D0");
double *const _dddgamma_D1D1D2D0D1 = get_v(dddgConf_stem,"D1D1D2D0D1");
double *const _dddgamma_D2D2D0D0D0 = get_v(dddgConf_stem,"D2D2D0D0D0");
double *const _dddgamma_D2D2D2D0D2 = get_v(dddgConf_stem,"D2D2D2D0D2");
double *const _dddgamma_D0D2D2D0D2 = get_v(dddgConf_stem,"D0D2D2D0D2");
double *const _dddgamma_D0D2D2D2D2 = get_v(dddgConf_stem,"D0D2D2D2D2");
double *const _dddgamma_D1D1D0D0D0 = get_v(dddgConf_stem,"D1D1D0D0D0");
double *const _dddgamma_D1D2D2D0D1 = get_v(dddgConf_stem,"D1D2D2D0D1");
double *const _dddgamma_D0D0D1D1D0 = get_v(dddgConf_stem,"D0D0D1D1D0");
double *const _dddgamma_D1D1D1D0D2 = get_v(dddgConf_stem,"D1D1D1D0D2");
double *const _dddgamma_D2D2D1D0D2 = get_v(dddgConf_stem,"D2D2D1D0D2");
double *const _dddgamma_D0D2D1D1D1 = get_v(dddgConf_stem,"D0D2D1D1D1");
double *const _dddgamma_D0D0D2D0D0 = get_v(dddgConf_stem,"D0D0D2D0D0");
double *const _dddgamma_D0D1D1D2D2 = get_v(dddgConf_stem,"D0D1D1D2D2");
double *const _dddgamma_D0D0D2D0D1 = get_v(dddgConf_stem,"D0D0D2D0D1");
double *const _dddgamma_D0D0D0D2D1 = get_v(dddgConf_stem,"D0D0D0D2D1");
double *const _dddgamma_D0D2D2D1D1 = get_v(dddgConf_stem,"D0D2D2D1D1");
double *const _dddgamma_D0D0D2D2D2 = get_v(dddgConf_stem,"D0D0D2D2D2");
double *const _dddgamma_D1D1D2D1D1 = get_v(dddgConf_stem,"D1D1D2D1D1");
double *const _dddgamma_D1D2D0D2D2 = get_v(dddgConf_stem,"D1D2D0D2D2");


  for (ijk = 0; ijk < nn; ++ijk) 
  {
    double x,y,z;
    x = patch->node[ijk]->x[0]-BH_center_x;
    y = patch->node[ijk]->x[1]-BH_center_y;
    z = patch->node[ijk]->x[2]-BH_center_z;
    KS_set_args

    _dddgamma_D0D1D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D0D1D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D1D2D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D1D1D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk1_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D2D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro );

    _dddgamma_D0D1D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D1D1D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro );

    _dddgamma_D0D2D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D1D1D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D1D1D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D0D2D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D1D1D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 4*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D1 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk0_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D1D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 4*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D1 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _dddgamma_D2D2D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk2_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D1D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro );

    _dddgamma_D2D2D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro );

    _dddgamma_D2D2D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
pow((bbn_ks_dk2_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D1D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D1D1D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D2D2D2 KS_func_pass_args_macro ) + 6*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 6*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 6*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D2 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D2 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D2D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D1D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D1D2D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D2D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D2D2D2 KS_func_pass_args_macro ) + 3*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 3*
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 3*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 3*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro );

    _dddgamma_D1D2D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D2D2D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D1D1D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro );

    _dddgamma_D1D2D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D1D2D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D0D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro );

    _dddgamma_D0D1D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D1D2D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D1D2D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D1D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D1D1D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk1_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D1D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D2D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D2D2D2 KS_func_pass_args_macro ) + 3*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 3*
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 3*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 3*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
3*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
3*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro );

    _dddgamma_D0D0D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro );

    _dddgamma_D0D0D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk0_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D1D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 4*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D1 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D2D2D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
pow((bbn_ks_dk2_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D0D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro );

    _dddgamma_D1D1D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk1_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D1D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D1D1D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
pow((bbn_ks_dk1_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D1D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk1_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D2D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D0D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk0_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D2D2D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk2_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D2D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro );

    _dddgamma_D0D2D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D0D0 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D0D0 KS_func_pass_args_macro ) + 3*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 3*
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 3*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 3*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D1D1D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
pow((bbn_ks_dk1_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D2D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D1D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro );

    _dddgamma_D1D2D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D2D2D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D1D1D1 KS_func_pass_args_macro ) + 6*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 6*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 6*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D0D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D0D0 KS_func_pass_args_macro ) + 6*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 6*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 6*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D0D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk0_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D1D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 4*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D1 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D1D1D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D1D1D1 KS_func_pass_args_macro ) + 6*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 6*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 6*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D1D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D0D1D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 4*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D1 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D1D2D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D0D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk0_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D2D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D0D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
pow((bbn_ks_dk0_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D2D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D1D1D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk1_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D0D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro );

    _dddgamma_D2D2D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
pow((bbn_ks_dk2_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D1D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro );

    _dddgamma_D0D1D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D1D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D1D1D1 KS_func_pass_args_macro ) + 3*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 3*
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 3*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 3*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
3*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
3*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D1D2D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D0D0 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D0D0 KS_func_pass_args_macro ) + 3*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 3*
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 3*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 3*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D1D2D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D1D1D1 KS_func_pass_args_macro ) + 6*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 6*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 6*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D2D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D1D2D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D0D1D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D0D0D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _dddgamma_D1D2D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D1D2D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D0D1D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D0D0 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D0D0 KS_func_pass_args_macro ) + 3*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 3*
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 3*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 3*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
3*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
3*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D1D2D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
pow((bbn_ks_dk0_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D0D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
pow((bbn_ks_dk0_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D0D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro );

    _dddgamma_D1D2D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk2_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D2D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D1D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D1D1D1 KS_func_pass_args_macro ) + 3*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 3*
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 3*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 3*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 4*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D1 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _dddgamma_D0D1D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D1D1D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro );

    _dddgamma_D0D1D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D0D1D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D1D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D1D2D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D2D2D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 4*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D1 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D1D1D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro );

    _dddgamma_D0D1D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D0D1D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D1D2D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D0D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro );

    _dddgamma_D0D1D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D0D1D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D1D1D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D2D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _dddgamma_D1D2D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D2D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D0D1D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D1D2D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D1D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D1D2D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D1D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 4*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D1 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D2D2D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro );

    _dddgamma_D1D1D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
pow((bbn_ks_dk1_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D1D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D2D2D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D2D2D2 KS_func_pass_args_macro ) + 6*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 6*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 6*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D2 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D2 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D2D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k0(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D2D2D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk2_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D2D2D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk2_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D2D2D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk2_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D1D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D2D2D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D0D0 KS_func_pass_args_macro ) + 6*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 6*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 6*
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D2D2D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro );

    _dddgamma_D0D2D2D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D2D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D2D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D2D2D2 KS_func_pass_args_macro ) + 3*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 3*
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 3*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 3*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro );

    _dddgamma_D1D1D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D0D0 KS_func_pass_args_macro ) + 6*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 6*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 6*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk1_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D2D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k1(x, y, z)*
(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + bbn_ks_k2(x, y, z)*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + (bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + (bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro );

    _dddgamma_D0D0D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 4*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D1 KS_func_pass_args_macro ), 2) + 4*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _dddgamma_D1D1D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D2D2D1D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*bbn_ks_c(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k2(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k2(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D0D2D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D1D1D1 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D1D1D1 KS_func_pass_args_macro ) + 3*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 3*
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 3*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 3*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
3*bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D0D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk0_D0 KS_func_pass_args_macro ), 2);

    _dddgamma_D0D1D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k1(x, y, z)*(bbn_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D2D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D0 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D0 KS_func_pass_args_macro )*
(bbn_ks_dk0_D1 KS_func_pass_args_macro );

    _dddgamma_D0D2D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*bbn_ks_c(x, y, z)*
(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 
bbn_ks_k0(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k0(x, y, z)*(bbn_ks_dk2_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D1 KS_func_pass_args_macro )*
(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D1 KS_func_pass_args_macro )*
(bbn_ks_dk2_D1 KS_func_pass_args_macro );

    _dddgamma_D0D0D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*bbn_ks_c(x, y, z)*bbn_ks_k0(x, y, z)*(bbn_ks_dddk0_D2D2D2 KS_func_pass_args_macro ) + 6*
bbn_ks_c(x, y, z)*(bbn_ks_dk0_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k0(x, y, z), 2)*(bbn_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 6*bbn_ks_k0(x, y, z)*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 6*
bbn_ks_k0(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk0_D2 KS_func_pass_args_macro ) + 
6*(bbn_ks_dc_D2 KS_func_pass_args_macro )*pow((bbn_ks_dk0_D2 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D1D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 4*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(bbn_ks_k1(x, y, z), 2)*(bbn_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*bbn_ks_k1(x, y, z)*
(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D1D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*bbn_ks_k1(x, y, z)*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*(bbn_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(bbn_ks_dc_D1 KS_func_pass_args_macro )*(bbn_ks_dk1_D1 KS_func_pass_args_macro )*
(bbn_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
pow((bbn_ks_dk1_D1 KS_func_pass_args_macro ), 2);

    _dddgamma_D1D2D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
bbn_ks_c(x, y, z)*bbn_ks_k1(x, y, z)*(bbn_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
bbn_ks_k2(x, y, z)*(bbn_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + bbn_ks_c(x, y, z)*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
bbn_ks_c(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_c(x, y, z)*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_c(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*bbn_ks_k2(x, y, z)*(bbn_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k1(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k1(x, y, z)*(bbn_ks_dk2_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
bbn_ks_k2(x, y, z)*(bbn_ks_ddc_D2D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*bbn_ks_k2(x, y, z)*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*(bbn_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(bbn_ks_dc_D0 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(bbn_ks_dc_D2 KS_func_pass_args_macro )*
(bbn_ks_dk1_D0 KS_func_pass_args_macro )*(bbn_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(bbn_ks_dc_D2 KS_func_pass_args_macro )*(bbn_ks_dk1_D2 KS_func_pass_args_macro )*
(bbn_ks_dk2_D0 KS_func_pass_args_macro );
  }
}


