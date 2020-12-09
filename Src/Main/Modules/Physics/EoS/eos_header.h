#ifndef EOS_LIB_H
#define EOS_LIB_H

#include "core_lib.h"
#include "physics_lib.h"
#include "physics_EoS_lib.h"

/* prefix for this module */
#define P_ "EOS_"


EoS_T *init_EoS(Physics_T *const phys);
void free_EoS(EoS_T *eos);

#endif

