#include "core_lib.h"
#include "physics_EoS_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "maths_spectral_methods_lib.h"
#include "physics_lib.h"
#include "TOV_lib.h"

/* global variable(s) */
static EoS_T *tov_eos = 0;

TOV_T *TOV_solution(TOV_T *const TOV);
TOV_T *TOV_init(void);
void TOV_free(TOV_T *TOV);
static void solve_ODE_enthalpy_approach(TOV_T *const TOV);
static double dr_dh(const double h,const double r, const double m);
static double dm_dh(const double h,const double r, const double m);
static double r_approx(const double h,const double h_c);
static double m_approx(const double h,const double h_c);
static void isotropic_coords_transformation(TOV_T *const TOV);
static double drbar_dh(const double h,const double rbar,const double r, const double m);
static double calculate_baryonic_mass(const TOV_T *const TOV);
static double *baryonic_mass_integrand(const TOV_T *const TOV);
static void calculate_phi(TOV_T *const TOV);
static void calculate_ADM_and_Komar_mass(TOV_T *const TOV);
static double *Komar_mass_integrand(const TOV_T *const TOV);
static double *ADM_mass_integrand(const TOV_T *const TOV);
static double c_rbar(TOV_T *const TOV);

