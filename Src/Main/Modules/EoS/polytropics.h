#include "core_lib.h"
#include "EoS_lib.h"

static unsigned find_threshold_number_h(const EoS_T *const eos);
double EoS_rho_h_pwp(EoS_T *const eos);
double EoS_p_h_pwp(EoS_T *const eos);
double EoS_e_h_pwp(EoS_T *const eos);
double EoS_rho_h_p(EoS_T *const eos);
double EoS_p_h_p(EoS_T *const eos);
double EoS_e_h_p(EoS_T *const eos);
double EoS_drho_dh_h_pwp(EoS_T *const eos);
double EoS_drho_dh_h_p(EoS_T *const eos);
double EoS_de_dh_h_pwp(EoS_T *const eos);
double EoS_de_dh_h_p(EoS_T *const eos);
