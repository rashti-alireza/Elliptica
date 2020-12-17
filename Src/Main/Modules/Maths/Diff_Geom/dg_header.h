#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "maths_special_functions_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"

void Christoffel_symbol_3d(Patch_T *const patch,const char *const ig,const char *const dg,const char *const Chris);
void dChristoffel_symbol_3d(Patch_T *const patch,const char *const dChris);

void Ricci_3d(Patch_T *const patch,const char *const ig,
              const char *const Chris,const char *const dChris,
              const char *const Ricci,const char *const trRicci);

void scale_to_BSSN_metric_3d(Patch_T *const patch,
                              const char *const given_g,
                              const char *const BSSN_g,
                              const char *const iBSSN_g,
                              const char *const Psi);

