#ifndef physics_Stress_Energy_LIB_H
#define physics_Stress_Energy_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PHYSICS_T;


int Tij_main(struct PHYSICS_T *const phys);
void Tij_NS_IF_XCTS_gConf_u0(Patch_T *const patch);

#endif


