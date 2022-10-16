#ifndef physics_EoS_LIB_H
#define physics_EoS_LIB_H
#include "elliptica_system_lib.h"

#define EOS_MAX_STR (400)

/* forward declaration */
struct PHYSICS_T;

/* struct for equation of states */
typedef struct EquationOfState_T
{
 struct PHYSICS_T *phys;
 char description[EOS_MAX_STR];
 char type[EOS_MAX_STR];
 char unit[EOS_MAX_STR];
 double *K;/* polytropic constant */
 double *rho0_th;/* thresholds of rest mass density */
 double *h_th;/* enthalpy thresholds */
 double *n;/* polytropic index n = 1/(gamma-1) */
 double *gamma;/* polytropic index */
 double *a;/* constant coefficient to ensure continuity */
 double h;/* enthalpy */
 Uint N;/* number of intervals i.e number of pieces */
 /* some flags for extra checks:
 // the default value is 0, so initially they are inactive, 
 // but when someone puts them 1, they become active. */
 Uint enthalpy_fatal: 1;/* checks for enthalpy values */
 double (*pressure)(struct EquationOfState_T *const eos);
 double (*energy_density)(struct EquationOfState_T *const eos);/* total energy density */
 double (*rest_mass_density)(struct EquationOfState_T *const eos);
 double (*specific_internal_energy)(struct EquationOfState_T *const eos);
 double (*de_dh)(struct EquationOfState_T *const eos);/* d(energy_density)/dh */
 double (*drho0_dh)(struct EquationOfState_T *const eos);/* d(rest_mass_density)/dh */
 /* save coeffs for spline interpolation */
 struct
 {
   double h_floor;/* set some themodynamics var to 0 if h < h_floor. */
   Uint sample_size;/* the number of points for a spline fit */
   double *h_sample;/* enthalpy sample values */
   double *p_sample;/* pressure sample values */
   double *e_sample;/* energy_density sample values */
   double *rho0_sample;/* rest_mass_density sample values */
   void *interp_p;/* spline interpolation struct for p */
   void *interp_e;/* spline interpolation struct for e */
   void *interp_rho0;/* spline interpolation struct for rho0 */
 }cubic_spline[1];
 
}EoS_T;
#undef EOS_MAX_STR

EoS_T *init_EoS(struct PHYSICS_T *const phys);
void free_EoS(EoS_T *eos);
void test_EoS(struct PHYSICS_T *const phys);


#endif


