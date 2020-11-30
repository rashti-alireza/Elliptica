#ifndef maths_diff_geom_LIB_H
#define maths_diff_geom_LIB_H

/* forward declaration */

struct PATCH_T;

void Christoffel_symbol_3d(struct PATCH_T *const patch,const char *const ig,const char *const dg,const char *const Chris);
void dChristoffel_symbol_3d(Patch_T *const patch,const char *const dChris);
void Ricci_3d(Patch_T *const patch,const char *const ig,
              const char *const Chris,const char *const dChris,
              const char *const Ricci,const char *const trRicci);



#endif

