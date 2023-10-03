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

//////////////////////////////////Interpolation method
double EoS_rho0_h_tab(EoS_T *const eos);
double EoS_p_h_tab(EoS_T *const eos);
double EoS_e_h_tab(EoS_T *const eos);
double EoS_e0_h_tab(EoS_T *const eos);
double EoS_rho0_h_tab(EoS_T *const eos);
double EoS_p_h_tab(EoS_T *const eos);
double EoS_e_h_tab(EoS_T *const eos);
double EoS_e0_h_tab(EoS_T *const eos);
double EoS_drho0_dh_h_tab(EoS_T *const eos);
double EoS_drho0_dh_h_tab(EoS_T *const eos);
double EoS_de_dh_h_tab(EoS_T *const eos);
double EoS_de_dh_h_tab(EoS_T *const eos);

/////////////////////////////////Root finder method
double EoS_enthalpy_def(void* const eos, const double* const params);
double EoS_rho0_RF(EoS_T *const eos);
double EoS_p_rho0_tab(EoS_T *const eos);
double EoS_e_rho0_tab(EoS_T *const eos);
double EoS_e0_rho0_tab(EoS_T *const eos);
double EoS_de_dh_RF(EoS_T *const eos);
double EoS_drho0_dh_RF(EoS_T *const eos);
/////////////////////////////////////////////////////
Uint get_sample_size(const char* const eos_file_name);
