#include "core_lib.h"
#include "memory_managing_lib.h"

double EoS_rho_h_pwp(EoS_T *const eos);
double EoS_p_h_pwp(EoS_T *const eos);
double EoS_e_h_pwp(EoS_T *const eos);

EoS_T *initialize_EoS(void);
void clean_EoS(EoS_T **eos);
static void populate_EoS(EoS_T *const eos);
static void plan_EoS(EoS_T *const eos);


