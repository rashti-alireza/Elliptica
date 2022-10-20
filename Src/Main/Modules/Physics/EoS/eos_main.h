#include "eos_header.h"
#include "maths_spectral_methods_lib.h"

#define MAX_STR (400)

double EoS_rho0_h_pwp(EoS_T *const eos);
double EoS_p_h_pwp(EoS_T *const eos);
double EoS_e_h_pwp(EoS_T *const eos);
double EoS_e0_h_pwp(EoS_T *const eos);
double EoS_rho0_h_p(EoS_T *const eos);
double EoS_p_h_p(EoS_T *const eos);
double EoS_e_h_p(EoS_T *const eos);
double EoS_e0_h_p(EoS_T *const eos);
double EoS_drho0_dh_h_pwp(EoS_T *const eos);
double EoS_drho0_dh_h_p(EoS_T *const eos);
double EoS_de_dh_h_pwp(EoS_T *const eos);
double EoS_de_dh_h_p(EoS_T *const eos);
double EoS_p_h_pwp_ncs(EoS_T *const eos);
double EoS_e_h_pwp_ncs(EoS_T *const eos);
double EoS_rho0_h_pwp_ncs(EoS_T *const eos);
EoS_T *init_EoS(Physics_T *const phys);
void free_EoS(EoS_T *eos);
static void populate_EoS(EoS_T *const eos);
static void fill_h_th(EoS_T *const eos);
static void fill_a(EoS_T *const eos);
static void fill_n(EoS_T *const eos);
static void fill_K(EoS_T *const eos);
static double *read_EoS_in_parameter_file(const char *const par,Uint *const N);





