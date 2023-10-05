#include "eos_header.h"
#include "maths_spectral_methods_lib.h"
#include "maths_equation_solvings_lib.h"
//#include "maths_calculus_lib.h"

// max number of eos table's row
#define EOS_MAX_NUM_ROW_TABLE (4000)

// max number of eos table's column char
#define EOS_MAX_NUM_COL_TABLE (1000)

void eos_tab_read_table(EoS_T* const eos);
void eos_tab_set_hermite(EoS_T* const eos);
void eos_tab_set_hermite_log(EoS_T* const eos);
static double logy_of_logh_hermite(EoS_T* const eos, 
                                 Interpolation_T *const interp_s,
                                 const double y_shift/* shifting constant */,
                                 const double y_floor);
static double p_of_h_hermite_log(EoS_T* const eos);
static double e_of_h_hermite_log(EoS_T* const eos);
static double rho0_of_h_hermite_log(EoS_T* const eos) __attribute__((unused));
static double rho0_e_p_h(EoS_T* const eos) __attribute__((unused));
static double e0_of_e_and_rho0(EoS_T* const eos);
static double dp_dh_hermite_log(EoS_T* const eos);
static double de_dh_hermite_log(EoS_T* const eos);
static double drho0_dh_e_p_h(EoS_T* const eos);
