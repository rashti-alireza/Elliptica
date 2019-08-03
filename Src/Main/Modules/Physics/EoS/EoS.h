#include "core_lib.h"
#include "memory_managing_lib.h"
#include "text_tools_lib.h"
#include "physics_EoS_lib.h"

#define MAX_STR 400

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
EoS_T *initialize_EoS(void);
void free_EoS(EoS_T *eos);
static void populate_EoS(EoS_T *const eos);
static void fill_h_th(EoS_T *const eos);
static void fill_a(EoS_T *const eos);
static void fill_n(EoS_T *const eos);
static double *read_EoS_in_parameter_file(const char *const par,unsigned *const N);





