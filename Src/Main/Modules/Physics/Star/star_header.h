#ifndef star_LIB_H
#define star_LIB_H

#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_EoS_lib.h"
#include "physics_observables_lib.h"
#include "managers_lib.h"
#include "fields_lib.h"
#include "maths_equation_solvings_lib.h"

/* root finder struc for Euler eq const */  
struct NS_Euler_eq_const_RootFinder_S
{
  Obj_Man_T *obj;
  double NS_baryonic_mass;
};

void star_idealfluid_NS_nonflat_find_Euler_const(Obj_Man_T *const obj);
double star_NS_baryonic_nonflat_mass(Obj_Man_T *const obj,const double Euler_C);

#endif

