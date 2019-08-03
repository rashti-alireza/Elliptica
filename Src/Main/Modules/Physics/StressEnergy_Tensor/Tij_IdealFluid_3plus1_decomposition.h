#include "core_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"
#include "memory_managing_lib.h"
#include "utilities_lib.h"
#include "EoS_lib.h"

void Tij_IF_u0(Patch_T *const patch);
void Tij_IF_build_psi6J_Ui(Patch_T *const patch);
void Tij_IF_build_psi6E(Patch_T *const patch);
void Tij_IF_build_psi6S(Patch_T *const patch);
void Tij_IF_build_psi6Sources(Grid_T *const grid);

