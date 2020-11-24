#include "star_header.h"


int star_NS_idealfluid_gConf_find_Euler_const(Physics_T *const phys);
int star_NS_idealfluid_extrapolate_matter_fields(Physics_T *const phys);
static double Euler_eq_const_gConf_rootfinder_eq(void *params,const double *const x);
static void W_spin_vector_idealfluid(Patch_T *const patch,const double Omega_NS[3],const double C_NS[3]);
void star_W_spin_vector_idealfluid_update(Physics_T *const phys);



