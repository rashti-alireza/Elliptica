#ifndef physics_star_LIB_H
#define physics_star_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PHYSICS_T;

int star_main(struct PHYSICS_T *const phys);
double star_NS_baryonic_gConf_mass(Physics_T *const phys,const double Euler_C);

void star_populate_psi_alphaPsi_matter_fields_TOV
      (Physics_T *const phys,const char *const region,
      const char *const Psi,const char *const AlphaPsi,
      const char *const Enthalpy,const char *const Rho0,
      const char *const Phi,const char *const W);
double star_NS_current_Euler_eq_const(Physics_T *const phys);

void 
star_print_properties
  (Physics_T *const phys,
  const char *const params,
  FILE *const file,
  const int pr_screen);

void star_W_spin_vector_idealfluid_update(Physics_T *const phys,const char *const region);
double star_NS_mass_shedding_indicator(Physics_T *const phys);


#endif


