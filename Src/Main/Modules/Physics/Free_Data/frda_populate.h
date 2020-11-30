#include "frda_header.h"
#include "maths_linear_algebra_lib.h"

void frda_populate_gConf_dgConf_igConf_KerrSchild(Physics_T *const phys);
void frda_compatible_Christoffel_symbol(Physics_T *const phys,const char *const ig,const char *const dg, const char *const Chris);
void frda_1st_derivative_Christoffel_symbol(Physics_T *const phys,const char *const dChris);

