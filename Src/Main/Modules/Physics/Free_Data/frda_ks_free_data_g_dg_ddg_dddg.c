#include "frda_header.h"
#include "frda_KerrSchild_header.h"
#undef x
#undef y
#undef z

#define STR_LEN (99)

#define get_v(stem,indices) \
 (sprintf(f_name,"%s_%s",stem,indices) && !free_coeffs(patch->fields[Ind(f_name)])?patch->fields[Ind(f_name)]->v:0)
 
#define KS_func_pass_args_sed KS_func_pass_args_macro

#define KS_set_args \
  struct Analytic_Func_Arg_S farg[1];\
  farg->x = x;\
  farg->y = y;\
  farg->z = z;\
  farg->X = frda_ks_X KS_func_pass_args_macro;\
  farg->Y = frda_ks_Y KS_func_pass_args_macro;\
  farg->Z = frda_ks_Z KS_func_pass_args_macro;\
  farg->R = frda_ks_R KS_func_pass_args_macro;\
  farg->dX_D0 = frda_ks_dX_D0 KS_func_pass_args_sed;\
  farg->dX_D1 = frda_ks_dX_D1 KS_func_pass_args_sed;\
  farg->dX_D2 = frda_ks_dX_D2 KS_func_pass_args_sed;\
  farg->dY_D0 = frda_ks_dY_D0 KS_func_pass_args_sed;\
  farg->dY_D1 = frda_ks_dY_D1 KS_func_pass_args_sed;\
  farg->dY_D2 = frda_ks_dY_D2 KS_func_pass_args_sed;\
  farg->dZ_D0 = frda_ks_dZ_D0 KS_func_pass_args_sed;\
  farg->dZ_D1 = frda_ks_dZ_D1 KS_func_pass_args_sed;\
  farg->dZ_D2 = frda_ks_dZ_D2 KS_func_pass_args_sed;

void frda_kerr_schild_g_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);
        
void frda_kerr_schild_dg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);
	
void frda_kerr_schild_ddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);

void frda_kerr_schild_dddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);

void frda_kerr_schild_g_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem)

{
  char f_name[STR_LEN];
  const unsigned nn = patch->nn;
  unsigned ijk;
    double *const _gamma_D2D2 = get_v(stem,"D2D2");
    double *const _gamma_D0D2 = get_v(stem,"D0D2");
    double *const _gamma_D0D0 = get_v(stem,"D0D0");
    double *const _gamma_D0D1 = get_v(stem,"D0D1");
    double *const _gamma_D1D2 = get_v(stem,"D1D2");
    double *const _gamma_D1D1 = get_v(stem,"D1D1");
    
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
frda_ks_c(x, y, z)*pow(frda_ks_k0(x, y, z), 2) + 1.0;
      _gamma_D0D1[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k1
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z);
      _gamma_D0D2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k0
// k2
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z);
      _gamma_D1D1[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k1
frda_ks_c(x, y, z)*pow(frda_ks_k1(x, y, z), 2) + 1.0;
      _gamma_D1D2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k1
// k2
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z);
      _gamma_D2D2[ijk] =
/* mcode in progress ... */
// Not supported in C:
// c
// k2
frda_ks_c(x, y, z)*pow(frda_ks_k2(x, y, z), 2) + 1.0;
    }
}
void frda_kerr_schild_dg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem)

{

  char f_name[STR_LEN];
  const unsigned nn = patch->nn;
  unsigned ijk;
  
double *const _dgamma_D1D2D2 = get_v(stem,"D1D2D2");
double *const _dgamma_D0D0D1 = get_v(stem,"D0D0D1");
double *const _dgamma_D0D2D1 = get_v(stem,"D0D2D1");
double *const _dgamma_D0D1D0 = get_v(stem,"D0D1D0");
double *const _dgamma_D1D2D1 = get_v(stem,"D1D2D1");
double *const _dgamma_D2D2D0 = get_v(stem,"D2D2D0");
double *const _dgamma_D0D0D0 = get_v(stem,"D0D0D0");
double *const _dgamma_D0D0D2 = get_v(stem,"D0D0D2");
double *const _dgamma_D0D2D2 = get_v(stem,"D0D2D2");
double *const _dgamma_D2D2D1 = get_v(stem,"D2D2D1");
double *const _dgamma_D0D1D1 = get_v(stem,"D0D1D1");
double *const _dgamma_D0D2D0 = get_v(stem,"D0D2D0");
double *const _dgamma_D1D2D0 = get_v(stem,"D1D2D0");
double *const _dgamma_D1D1D1 = get_v(stem,"D1D1D1");
double *const _dgamma_D0D1D2 = get_v(stem,"D0D1D2");
double *const _dgamma_D1D1D0 = get_v(stem,"D1D1D0");
double *const _dgamma_D1D1D2 = get_v(stem,"D1D1D2");
double *const _dgamma_D2D2D2 = get_v(stem,"D2D2D2");

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro );

    _dgamma_D0D0D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k0
(2*frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro ))*frda_ks_k0(x, y, z);

    _dgamma_D0D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro );

    _dgamma_D0D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro );

    _dgamma_D1D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro );

    _dgamma_D2D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k2
(2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro ))*frda_ks_k2(x, y, z);

    _dgamma_D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k0
(2*frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro ))*frda_ks_k0(x, y, z);

    _dgamma_D0D0D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k0
(2*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro ))*frda_ks_k0(x, y, z);

    _dgamma_D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro );

    _dgamma_D2D2D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k2
(2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro ))*frda_ks_k2(x, y, z);

    _dgamma_D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro );

    _dgamma_D0D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k2
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro );

    _dgamma_D1D2D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k1
// k2
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro );

    _dgamma_D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k1
(2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro ))*frda_ks_k1(x, y, z);

    _dgamma_D0D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// c
// k0
// k1
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro );

    _dgamma_D1D1D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k1
(2*frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro ))*frda_ks_k1(x, y, z);

    _dgamma_D1D1D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k1
(2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro ))*frda_ks_k1(x, y, z);

    _dgamma_D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// c
// k2
(2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro ))*frda_ks_k2(x, y, z);

  }
}

void frda_kerr_schild_ddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem)

{

  char f_name[STR_LEN];
  const unsigned nn = patch->nn;
  unsigned ijk;

double *const _ddgamma_D2D2D0D2 = get_v(stem,"D2D2D0D2");
double *const _ddgamma_D0D1D0D1 = get_v(stem,"D0D1D0D1");
double *const _ddgamma_D1D1D0D1 = get_v(stem,"D1D1D0D1");
double *const _ddgamma_D1D2D0D1 = get_v(stem,"D1D2D0D1");
double *const _ddgamma_D0D0D1D0 = get_v(stem,"D0D0D1D0");
double *const _ddgamma_D1D1D2D0 = get_v(stem,"D1D1D2D0");
double *const _ddgamma_D0D0D0D1 = get_v(stem,"D0D0D0D1");
double *const _ddgamma_D1D1D0D0 = get_v(stem,"D1D1D0D0");
double *const _ddgamma_D0D0D2D2 = get_v(stem,"D0D0D2D2");
double *const _ddgamma_D2D2D1D1 = get_v(stem,"D2D2D1D1");
double *const _ddgamma_D2D2D2D2 = get_v(stem,"D2D2D2D2");
double *const _ddgamma_D0D0D2D0 = get_v(stem,"D0D0D2D0");
double *const _ddgamma_D0D2D2D1 = get_v(stem,"D0D2D2D1");
double *const _ddgamma_D1D2D0D2 = get_v(stem,"D1D2D0D2");
double *const _ddgamma_D1D2D1D0 = get_v(stem,"D1D2D1D0");
double *const _ddgamma_D2D2D2D0 = get_v(stem,"D2D2D2D0");
double *const _ddgamma_D0D0D0D0 = get_v(stem,"D0D0D0D0");
double *const _ddgamma_D0D0D1D1 = get_v(stem,"D0D0D1D1");
double *const _ddgamma_D0D1D0D2 = get_v(stem,"D0D1D0D2");
double *const _ddgamma_D1D2D1D1 = get_v(stem,"D1D2D1D1");
double *const _ddgamma_D0D2D1D1 = get_v(stem,"D0D2D1D1");
double *const _ddgamma_D1D2D2D2 = get_v(stem,"D1D2D2D2");
double *const _ddgamma_D2D2D1D2 = get_v(stem,"D2D2D1D2");
double *const _ddgamma_D0D1D2D0 = get_v(stem,"D0D1D2D0");
double *const _ddgamma_D0D2D2D0 = get_v(stem,"D0D2D2D0");
double *const _ddgamma_D0D1D1D0 = get_v(stem,"D0D1D1D0");
double *const _ddgamma_D1D2D2D0 = get_v(stem,"D1D2D2D0");
double *const _ddgamma_D1D1D1D0 = get_v(stem,"D1D1D1D0");
double *const _ddgamma_D1D2D1D2 = get_v(stem,"D1D2D1D2");
double *const _ddgamma_D1D2D0D0 = get_v(stem,"D1D2D0D0");
double *const _ddgamma_D0D0D0D2 = get_v(stem,"D0D0D0D2");
double *const _ddgamma_D2D2D2D1 = get_v(stem,"D2D2D2D1");
double *const _ddgamma_D0D2D0D0 = get_v(stem,"D0D2D0D0");
double *const _ddgamma_D0D2D0D1 = get_v(stem,"D0D2D0D1");
double *const _ddgamma_D1D1D1D1 = get_v(stem,"D1D1D1D1");
double *const _ddgamma_D0D1D0D0 = get_v(stem,"D0D1D0D0");
double *const _ddgamma_D1D1D2D2 = get_v(stem,"D1D1D2D2");
double *const _ddgamma_D0D2D1D2 = get_v(stem,"D0D2D1D2");
double *const _ddgamma_D0D1D2D1 = get_v(stem,"D0D1D2D1");
double *const _ddgamma_D0D1D1D2 = get_v(stem,"D0D1D1D2");
double *const _ddgamma_D1D1D2D1 = get_v(stem,"D1D1D2D1");
double *const _ddgamma_D0D1D1D1 = get_v(stem,"D0D1D1D1");
double *const _ddgamma_D0D2D1D0 = get_v(stem,"D0D2D1D0");
double *const _ddgamma_D2D2D1D0 = get_v(stem,"D2D2D1D0");
double *const _ddgamma_D0D1D2D2 = get_v(stem,"D0D1D2D2");
double *const _ddgamma_D2D2D0D0 = get_v(stem,"D2D2D0D0");
double *const _ddgamma_D0D0D2D1 = get_v(stem,"D0D0D2D1");
double *const _ddgamma_D0D2D2D2 = get_v(stem,"D0D2D2D2");
double *const _ddgamma_D1D1D1D2 = get_v(stem,"D1D1D1D2");
double *const _ddgamma_D1D2D2D1 = get_v(stem,"D1D2D2D1");
double *const _ddgamma_D0D0D1D2 = get_v(stem,"D0D0D1D2");
double *const _ddgamma_D0D2D0D2 = get_v(stem,"D0D2D0D2");
double *const _ddgamma_D1D1D0D2 = get_v(stem,"D1D1D0D2");
double *const _ddgamma_D2D2D0D1 = get_v(stem,"D2D2D0D1");

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D1D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*pow((frda_ks_dk1_D0 KS_func_pass_args_macro ), 2) + pow(frda_ks_k1(x, y, z), 2)*
(frda_ks_ddc_D0D0 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

    _ddgamma_D0D0D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*pow((frda_ks_dk0_D2 KS_func_pass_args_macro ), 2) + pow(frda_ks_k0(x, y, z), 2)*
(frda_ks_ddc_D2D2 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro );

    _ddgamma_D2D2D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*pow((frda_ks_dk2_D1 KS_func_pass_args_macro ), 2) + pow(frda_ks_k2(x, y, z), 2)*
(frda_ks_ddc_D1D1 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

    _ddgamma_D2D2D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*pow((frda_ks_dk2_D2 KS_func_pass_args_macro ), 2) + pow(frda_ks_k2(x, y, z), 2)*
(frda_ks_ddc_D2D2 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

    _ddgamma_D0D0D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*pow((frda_ks_dk0_D0 KS_func_pass_args_macro ), 2) + pow(frda_ks_k0(x, y, z), 2)*
(frda_ks_ddc_D0D0 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D0D0D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k0
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*pow((frda_ks_dk0_D1 KS_func_pass_args_macro ), 2) + pow(frda_ks_k0(x, y, z), 2)*
(frda_ks_ddc_D1D1 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D1D1D1[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*pow((frda_ks_dk1_D1 KS_func_pass_args_macro ), 2) + pow(frda_ks_k1(x, y, z), 2)*
(frda_ks_ddc_D1D1 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

    _ddgamma_D1D1D2D2[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k1
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*pow((frda_ks_dk1_D2 KS_func_pass_args_macro ), 2) + pow(frda_ks_k1(x, y, z), 2)*
(frda_ks_ddc_D2D2 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro );

    _ddgamma_D2D2D0D0[ijk]=
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// c
// k2
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*pow((frda_ks_dk2_D0 KS_func_pass_args_macro ), 2) + pow(frda_ks_k2(x, y, z), 2)*
(frda_ks_ddc_D0D0 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );
  }	
}

void frda_kerr_schild_dddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem)

{

  char f_name[STR_LEN];
  const unsigned nn = patch->nn;
  unsigned ijk;

double *const _dddgamma_D0D1D2D2D1 = get_v(stem,"D0D1D2D2D1");
double *const _dddgamma_D0D1D2D0D1 = get_v(stem,"D0D1D2D0D1");
double *const _dddgamma_D1D2D0D2D0 = get_v(stem,"D1D2D0D2D0");
double *const _dddgamma_D1D1D1D2D1 = get_v(stem,"D1D1D1D2D1");
double *const _dddgamma_D0D2D0D0D1 = get_v(stem,"D0D2D0D0D1");
double *const _dddgamma_D0D2D0D2D2 = get_v(stem,"D0D2D0D2D2");
double *const _dddgamma_D2D2D2D1D2 = get_v(stem,"D2D2D2D1D2");
double *const _dddgamma_D0D1D0D2D0 = get_v(stem,"D0D1D0D2D0");
double *const _dddgamma_D1D1D2D2D0 = get_v(stem,"D1D1D2D2D0");
double *const _dddgamma_D0D2D0D2D1 = get_v(stem,"D0D2D0D2D1");
double *const _dddgamma_D1D1D2D1D0 = get_v(stem,"D1D1D2D1D0");
double *const _dddgamma_D1D1D0D1D2 = get_v(stem,"D1D1D0D1D2");
double *const _dddgamma_D0D2D1D1D0 = get_v(stem,"D0D2D1D1D0");
double *const _dddgamma_D1D1D1D0D1 = get_v(stem,"D1D1D1D0D1");
double *const _dddgamma_D0D0D2D1D1 = get_v(stem,"D0D0D2D1D1");
double *const _dddgamma_D0D1D2D1D2 = get_v(stem,"D0D1D2D1D2");
double *const _dddgamma_D0D0D0D1D1 = get_v(stem,"D0D0D0D1D1");
double *const _dddgamma_D2D2D1D2D1 = get_v(stem,"D2D2D1D2D1");
double *const _dddgamma_D1D1D1D2D2 = get_v(stem,"D1D1D1D2D2");
double *const _dddgamma_D2D2D2D2D1 = get_v(stem,"D2D2D2D2D1");
double *const _dddgamma_D2D2D1D0D0 = get_v(stem,"D2D2D1D0D0");
double *const _dddgamma_D0D1D2D1D1 = get_v(stem,"D0D1D2D1D1");
double *const _dddgamma_D1D1D2D2D2 = get_v(stem,"D1D1D2D2D2");
double *const _dddgamma_D1D2D1D0D1 = get_v(stem,"D1D2D1D0D1");
double *const _dddgamma_D0D1D1D1D0 = get_v(stem,"D0D1D1D1D0");
double *const _dddgamma_D2D2D2D0D1 = get_v(stem,"D2D2D2D0D1");
double *const _dddgamma_D1D2D2D2D2 = get_v(stem,"D1D2D2D2D2");
double *const _dddgamma_D1D2D1D2D1 = get_v(stem,"D1D2D1D2D1");
double *const _dddgamma_D2D2D0D2D1 = get_v(stem,"D2D2D0D2D1");
double *const _dddgamma_D1D1D2D1D2 = get_v(stem,"D1D1D2D1D2");
double *const _dddgamma_D1D2D2D1D0 = get_v(stem,"D1D2D2D1D0");
double *const _dddgamma_D1D2D0D0D1 = get_v(stem,"D1D2D0D0D1");
double *const _dddgamma_D0D2D0D1D1 = get_v(stem,"D0D2D0D1D1");
double *const _dddgamma_D0D2D1D0D1 = get_v(stem,"D0D2D1D0D1");
double *const _dddgamma_D0D0D0D1D2 = get_v(stem,"D0D0D0D1D2");
double *const _dddgamma_D0D0D2D1D2 = get_v(stem,"D0D0D2D1D2");
double *const _dddgamma_D0D1D0D1D0 = get_v(stem,"D0D1D0D1D0");
double *const _dddgamma_D1D2D1D0D0 = get_v(stem,"D1D2D1D0D0");
double *const _dddgamma_D0D2D1D0D0 = get_v(stem,"D0D2D1D0D0");
double *const _dddgamma_D1D2D0D1D1 = get_v(stem,"D1D2D0D1D1");
double *const _dddgamma_D0D1D2D0D2 = get_v(stem,"D0D1D2D0D2");
double *const _dddgamma_D1D1D2D0D0 = get_v(stem,"D1D1D2D0D0");
double *const _dddgamma_D0D1D2D2D2 = get_v(stem,"D0D1D2D2D2");
double *const _dddgamma_D0D0D2D0D2 = get_v(stem,"D0D0D2D0D2");
double *const _dddgamma_D0D0D0D0D2 = get_v(stem,"D0D0D0D0D2");
double *const _dddgamma_D1D1D0D1D1 = get_v(stem,"D1D1D0D1D1");
double *const _dddgamma_D2D2D0D0D1 = get_v(stem,"D2D2D0D0D1");
double *const _dddgamma_D0D0D1D2D2 = get_v(stem,"D0D0D1D2D2");
double *const _dddgamma_D1D1D0D2D0 = get_v(stem,"D1D1D0D2D0");
double *const _dddgamma_D0D1D1D2D0 = get_v(stem,"D0D1D1D2D0");
double *const _dddgamma_D1D1D0D0D1 = get_v(stem,"D1D1D0D0D1");
double *const _dddgamma_D1D1D0D0D2 = get_v(stem,"D1D1D0D0D2");
double *const _dddgamma_D1D2D0D1D2 = get_v(stem,"D1D2D0D1D2");
double *const _dddgamma_D0D0D1D0D2 = get_v(stem,"D0D0D1D0D2");
double *const _dddgamma_D0D0D0D2D0 = get_v(stem,"D0D0D0D2D0");
double *const _dddgamma_D2D2D2D1D1 = get_v(stem,"D2D2D2D1D1");
double *const _dddgamma_D0D2D0D0D2 = get_v(stem,"D0D2D0D0D2");
double *const _dddgamma_D2D2D2D2D0 = get_v(stem,"D2D2D2D2D0");
double *const _dddgamma_D0D2D0D0D0 = get_v(stem,"D0D2D0D0D0");
double *const _dddgamma_D2D2D1D2D0 = get_v(stem,"D2D2D1D2D0");
double *const _dddgamma_D1D1D0D1D0 = get_v(stem,"D1D1D0D1D0");
double *const _dddgamma_D1D2D2D0D2 = get_v(stem,"D1D2D2D0D2");
double *const _dddgamma_D0D1D0D1D1 = get_v(stem,"D0D1D0D1D1");
double *const _dddgamma_D2D2D0D2D2 = get_v(stem,"D2D2D0D2D2");
double *const _dddgamma_D1D2D2D1D2 = get_v(stem,"D1D2D2D1D2");
double *const _dddgamma_D2D2D1D1D1 = get_v(stem,"D2D2D1D1D1");
double *const _dddgamma_D0D0D0D0D0 = get_v(stem,"D0D0D0D0D0");
double *const _dddgamma_D0D0D1D2D1 = get_v(stem,"D0D0D1D2D1");
double *const _dddgamma_D1D1D1D1D0 = get_v(stem,"D1D1D1D1D0");
double *const _dddgamma_D1D1D1D1D1 = get_v(stem,"D1D1D1D1D1");
double *const _dddgamma_D1D1D1D2D0 = get_v(stem,"D1D1D1D2D0");
double *const _dddgamma_D0D1D0D2D1 = get_v(stem,"D0D1D0D2D1");
double *const _dddgamma_D2D2D1D1D0 = get_v(stem,"D2D2D1D1D0");
double *const _dddgamma_D1D2D2D2D0 = get_v(stem,"D1D2D2D2D0");
double *const _dddgamma_D0D0D1D1D2 = get_v(stem,"D0D0D1D1D2");
double *const _dddgamma_D0D2D0D2D0 = get_v(stem,"D0D2D0D2D0");
double *const _dddgamma_D0D0D1D0D0 = get_v(stem,"D0D0D1D0D0");
double *const _dddgamma_D0D2D2D1D0 = get_v(stem,"D0D2D2D1D0");
double *const _dddgamma_D0D2D2D0D0 = get_v(stem,"D0D2D2D0D0");
double *const _dddgamma_D1D1D1D1D2 = get_v(stem,"D1D1D1D1D2");
double *const _dddgamma_D0D0D2D2D1 = get_v(stem,"D0D0D2D2D1");
double *const _dddgamma_D2D2D0D1D0 = get_v(stem,"D2D2D0D1D0");
double *const _dddgamma_D1D1D2D0D2 = get_v(stem,"D1D1D2D0D2");
double *const _dddgamma_D0D1D1D1D1 = get_v(stem,"D0D1D1D1D1");
double *const _dddgamma_D1D2D0D0D0 = get_v(stem,"D1D2D0D0D0");
double *const _dddgamma_D1D2D0D0D2 = get_v(stem,"D1D2D0D0D2");
double *const _dddgamma_D2D2D0D1D2 = get_v(stem,"D2D2D0D1D2");
double *const _dddgamma_D0D0D1D1D1 = get_v(stem,"D0D0D1D1D1");
double *const _dddgamma_D0D2D0D1D0 = get_v(stem,"D0D2D0D1D0");
double *const _dddgamma_D0D2D1D2D1 = get_v(stem,"D0D2D1D2D1");
double *const _dddgamma_D1D2D1D2D2 = get_v(stem,"D1D2D1D2D2");
double *const _dddgamma_D0D1D1D0D0 = get_v(stem,"D0D1D1D0D0");
double *const _dddgamma_D0D0D1D2D0 = get_v(stem,"D0D0D1D2D0");
double *const _dddgamma_D1D2D1D1D0 = get_v(stem,"D1D2D1D1D0");
double *const _dddgamma_D1D2D2D2D1 = get_v(stem,"D1D2D2D2D1");
double *const _dddgamma_D0D1D0D0D0 = get_v(stem,"D0D1D0D0D0");
double *const _dddgamma_D1D2D2D1D1 = get_v(stem,"D1D2D2D1D1");
double *const _dddgamma_D0D0D0D1D0 = get_v(stem,"D0D0D0D1D0");
double *const _dddgamma_D0D0D0D0D1 = get_v(stem,"D0D0D0D0D1");
double *const _dddgamma_D0D0D0D2D2 = get_v(stem,"D0D0D0D2D2");
double *const _dddgamma_D1D2D0D2D1 = get_v(stem,"D1D2D0D2D1");
double *const _dddgamma_D2D2D1D1D2 = get_v(stem,"D2D2D1D1D2");
double *const _dddgamma_D1D2D1D1D1 = get_v(stem,"D1D2D1D1D1");
double *const _dddgamma_D0D0D1D0D1 = get_v(stem,"D0D0D1D0D1");
double *const _dddgamma_D0D1D0D1D2 = get_v(stem,"D0D1D0D1D2");
double *const _dddgamma_D0D2D1D2D2 = get_v(stem,"D0D2D1D2D2");
double *const _dddgamma_D1D1D0D2D2 = get_v(stem,"D1D1D0D2D2");
double *const _dddgamma_D0D1D0D0D2 = get_v(stem,"D0D1D0D0D2");
double *const _dddgamma_D0D1D0D2D2 = get_v(stem,"D0D1D0D2D2");
double *const _dddgamma_D0D2D1D2D0 = get_v(stem,"D0D2D1D2D0");
double *const _dddgamma_D0D1D0D0D1 = get_v(stem,"D0D1D0D0D1");
double *const _dddgamma_D1D2D1D1D2 = get_v(stem,"D1D2D1D1D2");
double *const _dddgamma_D2D2D0D1D1 = get_v(stem,"D2D2D0D1D1");
double *const _dddgamma_D1D1D2D2D1 = get_v(stem,"D1D1D2D2D1");
double *const _dddgamma_D0D1D2D1D0 = get_v(stem,"D0D1D2D1D0");
double *const _dddgamma_D0D2D1D0D2 = get_v(stem,"D0D2D1D0D2");
double *const _dddgamma_D2D2D2D1D0 = get_v(stem,"D2D2D2D1D0");
double *const _dddgamma_D0D1D1D1D2 = get_v(stem,"D0D1D1D1D2");
double *const _dddgamma_D1D2D1D0D2 = get_v(stem,"D1D2D1D0D2");
double *const _dddgamma_D0D0D2D2D0 = get_v(stem,"D0D0D2D2D0");
double *const _dddgamma_D0D1D1D0D1 = get_v(stem,"D0D1D1D0D1");
double *const _dddgamma_D0D1D2D2D0 = get_v(stem,"D0D1D2D2D0");
double *const _dddgamma_D0D2D0D1D2 = get_v(stem,"D0D2D0D1D2");
double *const _dddgamma_D0D2D2D2D1 = get_v(stem,"D0D2D2D2D1");
double *const _dddgamma_D1D1D0D2D1 = get_v(stem,"D1D1D0D2D1");
double *const _dddgamma_D0D0D2D1D0 = get_v(stem,"D0D0D2D1D0");
double *const _dddgamma_D1D2D0D1D0 = get_v(stem,"D1D2D0D1D0");
double *const _dddgamma_D0D2D2D1D2 = get_v(stem,"D0D2D2D1D2");
double *const _dddgamma_D0D1D1D2D1 = get_v(stem,"D0D1D1D2D1");
double *const _dddgamma_D1D2D2D0D0 = get_v(stem,"D1D2D2D0D0");
double *const _dddgamma_D0D1D2D0D0 = get_v(stem,"D0D1D2D0D0");
double *const _dddgamma_D1D2D1D2D0 = get_v(stem,"D1D2D1D2D0");
double *const _dddgamma_D2D2D1D0D1 = get_v(stem,"D2D2D1D0D1");
double *const _dddgamma_D2D2D1D2D2 = get_v(stem,"D2D2D1D2D2");
double *const _dddgamma_D1D1D1D0D0 = get_v(stem,"D1D1D1D0D0");
double *const _dddgamma_D0D1D1D0D2 = get_v(stem,"D0D1D1D0D2");
double *const _dddgamma_D0D2D2D2D0 = get_v(stem,"D0D2D2D2D0");
double *const _dddgamma_D2D2D2D2D2 = get_v(stem,"D2D2D2D2D2");
double *const _dddgamma_D0D2D2D0D1 = get_v(stem,"D0D2D2D0D1");
double *const _dddgamma_D0D2D1D1D2 = get_v(stem,"D0D2D1D1D2");
double *const _dddgamma_D2D2D0D0D2 = get_v(stem,"D2D2D0D0D2");
double *const _dddgamma_D2D2D2D0D0 = get_v(stem,"D2D2D2D0D0");
double *const _dddgamma_D2D2D0D2D0 = get_v(stem,"D2D2D0D2D0");
double *const _dddgamma_D1D1D2D0D1 = get_v(stem,"D1D1D2D0D1");
double *const _dddgamma_D2D2D0D0D0 = get_v(stem,"D2D2D0D0D0");
double *const _dddgamma_D2D2D2D0D2 = get_v(stem,"D2D2D2D0D2");
double *const _dddgamma_D0D2D2D0D2 = get_v(stem,"D0D2D2D0D2");
double *const _dddgamma_D0D2D2D2D2 = get_v(stem,"D0D2D2D2D2");
double *const _dddgamma_D1D1D0D0D0 = get_v(stem,"D1D1D0D0D0");
double *const _dddgamma_D1D2D2D0D1 = get_v(stem,"D1D2D2D0D1");
double *const _dddgamma_D0D0D1D1D0 = get_v(stem,"D0D0D1D1D0");
double *const _dddgamma_D1D1D1D0D2 = get_v(stem,"D1D1D1D0D2");
double *const _dddgamma_D2D2D1D0D2 = get_v(stem,"D2D2D1D0D2");
double *const _dddgamma_D0D2D1D1D1 = get_v(stem,"D0D2D1D1D1");
double *const _dddgamma_D0D0D2D0D0 = get_v(stem,"D0D0D2D0D0");
double *const _dddgamma_D0D1D1D2D2 = get_v(stem,"D0D1D1D2D2");
double *const _dddgamma_D0D0D2D0D1 = get_v(stem,"D0D0D2D0D1");
double *const _dddgamma_D0D0D0D2D1 = get_v(stem,"D0D0D0D2D1");
double *const _dddgamma_D0D2D2D1D1 = get_v(stem,"D0D2D2D1D1");
double *const _dddgamma_D0D0D2D2D2 = get_v(stem,"D0D0D2D2D2");
double *const _dddgamma_D1D1D2D1D1 = get_v(stem,"D1D1D2D1D1");
double *const _dddgamma_D1D2D0D2D2 = get_v(stem,"D1D2D0D2D2");


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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk1_D1 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 4*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk1_D1 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk0_D1 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 4*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk0_D1 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk2_D1 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
pow((frda_ks_dk2_D0 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D2D2D2 KS_func_pass_args_macro ) + 6*
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 6*frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 6*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D2 KS_func_pass_args_macro )*pow((frda_ks_dk1_D2 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D2D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D2D2D2 KS_func_pass_args_macro ) + 3*frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 3*
frda_ks_c(x, y, z)*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 3*frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 3*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk1_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk1_D0 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D2D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D2D2D2 KS_func_pass_args_macro ) + 3*frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 3*
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 3*frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 3*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
3*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
3*frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk0_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 4*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk1_D1 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
pow((frda_ks_dk2_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk1_D0 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
pow((frda_ks_dk1_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk1_D0 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk0_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk2_D1 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D0D0 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D0D0 KS_func_pass_args_macro ) + 3*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 3*
frda_ks_c(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 3*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 3*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
pow((frda_ks_dk1_D0 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D1D1D1 KS_func_pass_args_macro ) + 6*
frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 6*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 6*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk2_D1 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D0D0 KS_func_pass_args_macro ) + 6*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 6*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 6*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk0_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk0_D1 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 4*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk1_D1 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D1D1D1 KS_func_pass_args_macro ) + 6*
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 6*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 6*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk1_D1 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 4*
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk2_D1 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk0_D1 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
pow((frda_ks_dk0_D0 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk1_D1 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
pow((frda_ks_dk2_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D1D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D1D1D1 KS_func_pass_args_macro ) + 3*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 3*
frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 3*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 3*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
3*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
3*frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D0D0 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D0D0 KS_func_pass_args_macro ) + 3*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 3*
frda_ks_c(x, y, z)*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 3*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 3*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D1D1D1 KS_func_pass_args_macro ) + 6*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 6*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 6*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk0_D1 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D0D0 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D0D0 KS_func_pass_args_macro ) + 3*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 3*
frda_ks_c(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 3*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 3*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
3*frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
3*frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk1_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
pow((frda_ks_dk0_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
pow((frda_ks_dk0_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk2_D1 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D1D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D1D1D1 KS_func_pass_args_macro ) + 3*frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 3*
frda_ks_c(x, y, z)*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 3*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 3*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 4*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk0_D1 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk1_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 4*
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk2_D1 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk1_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk0_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 4*
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk2_D1 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D1D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*pow((frda_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D0D1 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D0D1 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
pow((frda_ks_dk1_D0 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D2D2D2 KS_func_pass_args_macro ) + 6*
frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 6*frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 6*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D2 KS_func_pass_args_macro )*pow((frda_ks_dk2_D2 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k0(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk2_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk2_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D0D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk2_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D0D0 KS_func_pass_args_macro ) + 6*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 6*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D0 KS_func_pass_args_macro ) + 6*
frda_ks_k2(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk2_D0 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 4*
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk2_D2 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D2D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D2D2D2 KS_func_pass_args_macro ) + 3*frda_ks_c(x, y, z)*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 3*
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 3*frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 3*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D0D0 KS_func_pass_args_macro ) + 6*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D0D0 KS_func_pass_args_macro ) + 6*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D0 KS_func_pass_args_macro ) + 6*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk1_D0 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k1(x, y, z)*
(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + frda_ks_k2(x, y, z)*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + (frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro ) + (frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D1D1 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D1D1 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 4*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*pow((frda_ks_dk0_D1 KS_func_pass_args_macro ), 2) + 4*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k1(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddk2_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k2(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k2(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D1D1D1 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D1D1D1 KS_func_pass_args_macro ) + 3*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 3*
frda_ks_c(x, y, z)*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D1D1 KS_func_pass_args_macro ) + 3*frda_ks_k0(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 3*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
3*frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D0D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D0D2 KS_func_pass_args_macro ) + 4*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D0D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D0 KS_func_pass_args_macro ) + 
4*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk0_D0 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk1_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k1(x, y, z)*(frda_ks_dddk0_D1D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddc_D1D2D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D0D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D0D1D2 KS_func_pass_args_macro ) + 2*frda_ks_k0(x, y, z)*
(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D0D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D1 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D0 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D0 KS_func_pass_args_macro )*
(frda_ks_dk0_D1 KS_func_pass_args_macro );

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
frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk2_D1D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk0_D1D1D2 KS_func_pass_args_macro ) + 2*frda_ks_c(x, y, z)*
(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_ddk0_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k0(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 
frda_ks_k0(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k0(x, y, z)*(frda_ks_dk2_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D1D1 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk0_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D1 KS_func_pass_args_macro )*
(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_dk2_D1 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk0_D1 KS_func_pass_args_macro )*
(frda_ks_dk2_D1 KS_func_pass_args_macro );

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
2*frda_ks_c(x, y, z)*frda_ks_k0(x, y, z)*(frda_ks_dddk0_D2D2D2 KS_func_pass_args_macro ) + 6*
frda_ks_c(x, y, z)*(frda_ks_dk0_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k0(x, y, z), 2)*(frda_ks_dddc_D2D2D2 KS_func_pass_args_macro ) + 6*frda_ks_k0(x, y, z)*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk0_D2D2 KS_func_pass_args_macro ) + 6*
frda_ks_k0(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk0_D2 KS_func_pass_args_macro ) + 
6*(frda_ks_dc_D2 KS_func_pass_args_macro )*pow((frda_ks_dk0_D2 KS_func_pass_args_macro ), 2);

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
2*frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk1_D1D1D2 KS_func_pass_args_macro ) + 4*
frda_ks_c(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
pow(frda_ks_k1(x, y, z), 2)*(frda_ks_dddc_D1D1D2 KS_func_pass_args_macro ) + 4*frda_ks_k1(x, y, z)*
(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D2 KS_func_pass_args_macro ) + 2*
frda_ks_k1(x, y, z)*(frda_ks_ddc_D1D1 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D1D1 KS_func_pass_args_macro ) + 
4*frda_ks_k1(x, y, z)*(frda_ks_dk1_D1 KS_func_pass_args_macro )*(frda_ks_ddc_D1D2 KS_func_pass_args_macro ) + 
4*(frda_ks_dc_D1 KS_func_pass_args_macro )*(frda_ks_dk1_D1 KS_func_pass_args_macro )*
(frda_ks_dk1_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
pow((frda_ks_dk1_D1 KS_func_pass_args_macro ), 2);

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
frda_ks_c(x, y, z)*frda_ks_k1(x, y, z)*(frda_ks_dddk2_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
frda_ks_k2(x, y, z)*(frda_ks_dddk1_D0D2D2 KS_func_pass_args_macro ) + frda_ks_c(x, y, z)*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 2*
frda_ks_c(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_c(x, y, z)*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_c(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*frda_ks_k2(x, y, z)*(frda_ks_dddc_D0D2D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk2_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk2_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k1(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk2_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k1(x, y, z)*(frda_ks_dk2_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_ddk1_D2D2 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_ddk1_D0D2 KS_func_pass_args_macro ) + 
frda_ks_k2(x, y, z)*(frda_ks_ddc_D2D2 KS_func_pass_args_macro )*(frda_ks_dk1_D0 KS_func_pass_args_macro ) + 
2*frda_ks_k2(x, y, z)*(frda_ks_dk1_D2 KS_func_pass_args_macro )*(frda_ks_ddc_D0D2 KS_func_pass_args_macro ) + 
2*(frda_ks_dc_D0 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*(frda_ks_dc_D2 KS_func_pass_args_macro )*
(frda_ks_dk1_D0 KS_func_pass_args_macro )*(frda_ks_dk2_D2 KS_func_pass_args_macro ) + 2*
(frda_ks_dc_D2 KS_func_pass_args_macro )*(frda_ks_dk1_D2 KS_func_pass_args_macro )*
(frda_ks_dk2_D0 KS_func_pass_args_macro );
  }
}


