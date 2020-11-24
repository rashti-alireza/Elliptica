#ifndef Tij_LIB_H
#define Tij_LIB_H

#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_EoS_lib.h"
#include "physics_lib.h"
#include "fields_lib.h"

int Tij_main(Physics_T *const phys);
int Tij_mount(Grid_T *const grid);
void Tij_NS_idealfluid_CTS_gConf_add_fields(Grid_T *const grid);
void Tij_NS_idealfluid_CTS_gConf_update(Physics_T *const phys);
void Tij_NS_IF_CTS_gConf_enthalpy(Patch_T *const patch,const double Euler_C);
void Tij_NS_IF_CTS_gConf_u0(Patch_T *const patch);
void Tij_NS_IF_CTS_gConf_psi6J_Ui(Patch_T *const patch);
void Tij_NS_IF_CTS_gConf_psi6E(Patch_T *const patch);
void Tij_NS_IF_CTS_gConf_psi6S(Patch_T *const patch);
void Tij_NS_IF_CTS_gConf_derives(Patch_T *const patch);
void Tij_NS_eos_update_rho0(Patch_T *const patch);
void Tij_NS_neat_enthalpy(Patch_T *const patch);


#endif

