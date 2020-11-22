#ifndef star_LIB_H
#define star_LIB_H

#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_EoS_lib.h"
#include "physics_observe_lib.h"
#include "physics_lib.h"
#include "fields_lib.h"
#include "maths_equation_solvings_lib.h"

/* root finder struc for Euler eq const */  
struct NS_Euler_eq_const_RootFinder_S
{
  Physics_T *phys;
  double NS_baryonic_mass;
};

void star_idealfluid_NS_gConf_find_Euler_const(Physics_T *const phys);
double star_NS_baryonic_gConf_mass(Physics_T *const phys,const double Euler_C);
void star_NS_idealfluid_gConf_add_fields(Grid_T *const grid);

#endif

