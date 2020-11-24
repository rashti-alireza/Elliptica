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

int star_NS_idealfluid_gConf_find_Euler_const(Physics_T *const phys);
double star_NS_baryonic_gConf_mass(Physics_T *const phys,const double Euler_C);
void star_NS_idealfluid_gConf_add_fields(Grid_T *const grid);
int star_NS_idealfluid_extrapolate_matter_fields(Physics_T *const phys);

int 
star_extrapolate
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  );
  

double star_NS_mass_shedding_indicator(Physics_T *const phys);
int star_NS_find_star_surface(Physics_T *const phys);

#endif

