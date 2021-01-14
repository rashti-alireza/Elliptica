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
#include "maths_spectral_methods_lib.h"
#include "TOV_lib.h"

/* parameter prefix 
// PLEASE keep it capitalized. */
#define P_ "STAR_"


/* root finder struc for Euler eq const */  
struct NS_Euler_eq_const_RootFinder_S
{
  Physics_T *phys;
  double NS_baryonic_mass;
};

/* root finder structure for NS center */
struct NC_Center_RootFinder_S
{
  Patch_T *patch;
  Patch_T **patches;
  Uint Np;/* number of patches */
  Root_Finder_T *root_finder;
};

/* root finder struc for force balance equation */  
struct Force_Balance_RootFinder_S
{
  Patch_T *patch;
  //Root_Finder_T *root_finder;
  double dLnGamma;
  double y_CM;
  double x_CM;
  double Vr;
  double D;
  double Omega;
  const double *X;
  const double *V2CM;
  Uint find_y_CM: 1;
  Uint find_x_CM: 1;
  Uint find_Omega: 1;
  int dir;/* direction of derivative */
};



int star_NS_ifluid_gConf_find_EulerC_fix_baryon_mass(Physics_T *const phys);
double star_NS_baryonic_gConf_mass(Physics_T *const phys,const double Euler_C);
void star_NS_idealfluid_gConf_add_fields(Grid_T *const grid);
int star_NS_idealfluid_extrapolate_matter_fields(Physics_T *const phys);

int 
star_NS_extrapolate
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  );
  

double star_NS_mass_shedding_indicator(Physics_T *const phys);
int star_NS_find_star_surface(Physics_T *const phys);
int star_NS_idealfluid_gConf_force_balance(Physics_T *const phys);
double star_NS_idealfluid_gConf_dLnGamma_force_bal(Patch_T *const patch,const double *const NS_centerX,const int dir);
double star_NS_idealfluid_gConf_root_force_bal(void *params,const double *const x);
void star_NS_find_where_denthalpy_is_0(Physics_T *const phys,double xdh0[3]);
int star_NS_keep_center_fixed(Physics_T *const phys);
void star_start_off_TOV (Physics_T *const phys);
double star_NS_current_Euler_eq_const(Physics_T *const phys);
void star_W_spin_vector_idealfluid_update(Physics_T *const phys,const char *const region);

#endif

