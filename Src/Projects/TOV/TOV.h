#include "core_lib.h"
#include "memory_managing_lib.h"
#include "EoS_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"

TOV_T *TOV_solution(TOV_T *const TOV);
static void solve_ODE_enthalpy_approach(TOV_T *const TOV);
static double dr_dh(const double h,const double r, const double m);
static double dm_dh(const double h,const double r, const double m);
static double calculate_baryonic_mass(const TOV_T *const TOV);
static double *baryonic_mass_integrand(const TOV_T *const TOV);
static void calculate_phi(TOV_T *const TOV);

