#include "eos_header.h"
#include "maths_spectral_methods_lib.h"

static Uint find_threshold_number_h(const EoS_T *const eos);
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

