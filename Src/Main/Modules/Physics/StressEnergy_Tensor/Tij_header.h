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

int Tij_tune(Obj_Man_T *const obj);
int Tij_mount(Grid_T *const grid);
void Tij_idealfluid_CTS_nonflat_add_fields(Grid_T *const grid);
void Tij_idealfluid_CTS_nonflat_update(Obj_Man_T *const obj);
void Tij_IF_CTS_nonflat_enthalpy(Patch_T *const patch,const double Euler_C);
void Tij_IF_CTS_nonflat_u0(Patch_T *const patch);
void Tij_IF_CTS_nonflat_psi6J_Ui(Patch_T *const patch);
void Tij_IF_CTS_nonflat_psi6E(Patch_T *const patch);
void Tij_IF_CTS_nonflat_psi6S(Patch_T *const patch);
void Tij_IF_CTS_nonflat_derives(Patch_T *const patch);
void Tij_eos_update_rho0(Patch_T *const patch);
void Tij_neat_enthalpy(Patch_T *const patch);


#endif

