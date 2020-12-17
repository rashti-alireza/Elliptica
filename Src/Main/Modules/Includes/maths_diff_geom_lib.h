#ifndef maths_diff_geom_LIB_H
#define maths_diff_geom_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PATCH_T;

void Christoffel_symbol_3d(struct PATCH_T *const patch,const char *const ig,const char *const dg,const char *const Chris);
void dChristoffel_symbol_3d(Patch_T *const patch,const char *const dChris);
void Ricci_3d(Patch_T *const patch,const char *const ig,
              const char *const Chris,const char *const dChris,
              const char *const Ricci,const char *const trRicci);


void scale_to_BSSN_metric_3d(struct PATCH_T *const patch,
                              const char *const given_g,
                              const char *const BSSN_g,
                              const char *const iBSSN_g,
                              const char *const Psi);


#endif

