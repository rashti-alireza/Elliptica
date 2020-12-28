#ifndef physics_star_LIB_H
#define physics_star_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PHYSICS_T;

int star_main(struct PHYSICS_T *const phys);
double star_NS_baryonic_gConf_mass(Physics_T *const phys,const double Euler_C);

#endif


