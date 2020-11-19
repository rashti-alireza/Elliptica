#ifndef Tij_LIB_H
#define Tij_LIB_H

#include "core_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_EoS_lib.h"
#include "managers_lib.h"
#include "fields_lib.h"

int Tij_update(Obj_Man_T *const obj);
int Tij_mount(Obj_Man_T *const obj);
void Tij_ideal_fluid_CTS_add_fields(Obj_Man_T *const obj);

#endif

