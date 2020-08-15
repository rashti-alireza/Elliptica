#ifndef bbn_headers_LIB_H
#define bbn_headers_LIB_H



#include "core_lib.h"
#include "TOV_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "maths_approximation_lib.h"
#include "maths_calculus_lib.h"
#include "physics_EoS_lib.h"
#include "physics_StressEnergyTensor_lib.h"
#include "physics_transformation_lib.h"
#include "fields_lib.h"
#include "physics_observables_lib.h"

/*----------------------------------------------------------------------*/
/* some sturcture for different purposes                                */
/*----------------------------------------------------------------------*/

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
  double Omega_BHNS;
  const double *X;
  unsigned find_y_CM: 1;
  unsigned find_x_CM: 1;
  unsigned find_Omega: 1;
  int dir;/* direction of derivative */
};

/* items needed to calculate Ps and Js ADMs */
struct items_S
{
  Patch_T *patch;/* the patch in which the following variables are defined */
  /* physical metric components */
  double *g00;
  double *g01;
  double *g02;
  double *g11;
  double *g12;
  double *g22;
  /* normal vector at the surface S, outward */
  double *n_U0;
  double *n_U1;
  double *n_U2;
  /* integration flags */
  unsigned surface_integration_flg: 1;/* if 1 means it measn 
                                      // we need surface integration 
                                      // on this patch as well, 
                                      // 0 means, no need */
  /* which hypersurface the surface integral is carried out */
  unsigned X_surface: 1;
  unsigned Y_surface: 1;
  unsigned Z_surface: 1;
  /* index of hypersurface for each X,Y and Z respectively */
  unsigned I;
  unsigned J;
  unsigned K;
};

void bbn_free_data_Gamma_patch(Patch_T *const patch);
void bbn_study_initial_data(Grid_T *const grid);
void bbn_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder);
void bbn_update_psi10A_UiUj(Patch_T *const patch);
void bbn_populate_free_data(Grid_T *const grid);
void bbn_free_data_gammas(Grid_T *const grid);
void bbn_free_conformal_metric_derivatives(Patch_T *const patch);
void bbn_preparing_conformal_metric_derivatives(Patch_T *const patch);
void bbn_free_data_Gamma(Grid_T *const grid);
void bbn_free_data_Ricci(Grid_T *const grid);
void bbn_free_data_tr_KSKij(Grid_T *const grid);
void bbn_free_data_KS_trKij(Patch_T *const patch);
void bbn_free_data_dGamma(Grid_T *const grid);
double bbn_NS_baryonic_mass(Grid_T *const grid,const double Euler_C);
double force_balance_root_finder_eq(void *params,const double *const x);
double dLnGamma_in_force_balance_eq(Patch_T *const patch,const double *const NS_centerX,const int dir);
void bbn_calculate_constraints_1st(Grid_T *const grid);
void bbn_calculate_constraints_2nd(Grid_T *const grid);
void bbn_update_enthalpy_and_denthalpy(Grid_T *const grid);
void bbn_update_stress_energy_tensor(Grid_T *const grid,const int flag);
void bbn_extrapolate_metric_fields_insideBH(Grid_T *const grid);
double bbn_BH_irreducible_mass(Grid_T *const grid);
void bbn_free_metric_and_Gamma_and_derivatives(Grid_T *const grid);
void bbn_make_metric_and_Gamma_and_derivatives(Grid_T *const grid);
void bbn_make_K_UiUj_and_dK_UiUj(Grid_T *const grid);
void bbn_free_K_UiUj_and_dK_UiUj(Grid_T *const grid);
void bbn_make_normal_vector_on_BH_horizon(Grid_T *const grid);
void bbn_add_fields(Grid_T *const grid);
void bbn_partial_derivatives_fields(Grid_T *const grid);
void bbn_update_derivative_phi(Patch_T *const patch);
void bbn_update_derivative_enthalpy(Patch_T *const patch);
void bbn_update_derivative_phi(Patch_T *const patch);
void bbn_update_rho0(Patch_T *const patch);
void bbn_update_derivative_rho0(Patch_T *const patch);
void bbn_update_derivative_u0(Patch_T *const patch);
void bbn_update_derivative_HS(Patch_T *const patch);
void bbn_update_derivative_K(Patch_T *const patch);
void bbn_update_derivative_psi(Patch_T *const patch);
void bbn_update_derivative_eta(Patch_T *const patch);
void bbn_update_derivative_B0_U0(Patch_T *const patch);
void bbn_update_derivative_B0_U1(Patch_T *const patch);
void bbn_update_derivative_B0_U2(Patch_T *const patch);
void bbn_update_derivative_Beta_U0(Patch_T *const patch);
void bbn_update_derivative_Beta_U1(Patch_T *const patch);
void bbn_update_derivative_Beta_U2(Patch_T *const patch);
void bbn_update_Beta_U0(Patch_T *const patch);
void bbn_update_Beta_U1(Patch_T *const patch);
void bbn_update_Beta_U2(Patch_T *const patch);
void bbn_update_enthalpy_and_denthalpy(Grid_T *const grid);
void bbn_update_stress_energy_tensor(Grid_T *const grid,const int flag);
void bbn_update_derivative_B1_U0(Patch_T *const patch);
void bbn_update_derivative_B1_U1(Patch_T *const patch);
void bbn_update_derivative_B1_U2(Patch_T *const patch);
void bbn_update_B1_U012(Patch_T *const patch);
void bbn_update_Aij(Grid_T *const grid);
Grid_T *bbn_init_from_checkpoint(FILE *const file);
void bbn_write_checkpoint(Grid_T *const grid);
Parameter_T *bbn_parameter_query_from_checkpoint_file(const char *const par_name,FILE *const file);
int bbn_IsCheckpointFileCompleted(const char *const file_path);
void bbn_add_fields_in_patch(Patch_T *const patch);
void bbn_bam_set_bam_fields(Grid_T *const grid);
void bbn_print_properties(Grid_T *const grid,const unsigned iteration, const char *const folder,const char *const open_file_mode,const int pr_flg);
void bbn_plan_obs_CS(Observable_T *obs);
void bbn_free_obs_CS(Observable_T *obs);
void bbn_measures(Grid_T *const grid);
void bbn_Rc_NS(double Rc[3],Grid_T *const grid);
void bbn_Rc_BH(double Rc[3],Grid_T *const grid);
void bbn_define_spin_JRP(double S[3],Grid_T *const grid,const char *const kind);
void bbn_define_spin_integral(double S[3],Grid_T *const grid,const char *const kind);
double bbn_BH_ADM_mass(Observable_T *const obs);
void bbn_transform_populate_boost_rotation(Transformation_T *const tB,Transformation_T *const tR);
void bbn_transform_get_k_and_H_KerrSchild(const double x,const double y,const double z,const double a,const double m,Transformation_T *const tB,Transformation_T *const tR,double *const kt,double *const k0,double *const k1,double *const k2,double *const H);
void bbn_populate_spin_integrands_Campanelli(Patch_T *const patch,const double xc[3],const double *const normal[3]);
double bbn_mass_shedding_indicator(Grid_T *const grid);
void bbn_test_induced_metric_algorithm(Grid_T *const grid);

void
bbn_compute_AKV_from_z
  (
  Grid_T *const grid/* grid */,
  const double *const akv/* akv scalar values */,
  const char *const dakv_D0/* d/dx akv name */,
  const char *const dakv_D1/* d/dy akv name */,
  const char *const dakv_D2/* d/dz akv name */,
  const char *const type/* NS or BH */,
  const unsigned Ntheta/* number of points in theta direction */,
  const unsigned Nphi/* number of points in theta direction */,
  const unsigned lmax/* l max in Ylm, if asked for spherical harmonic */,
  const int interpolation_type/* 1 double fourier, 0: spherical harmonic */
  );

void
bbn_compute_induced_metric_on_S2_CS_CTS
  (
  Grid_T *const grid/* grid */,
  const char *const type/* NS or BH */,
  const unsigned Ntheta/* number of points in theta direction */,
  const unsigned Nphi/* number of points in phi direction */,
  const unsigned lmax/* l max in Ylm */,
  double **const ph_D0D0/* induced h00  pointer */,
  double **const ph_D0D1/* induced h01  pointer */,
  double **const ph_D1D1/* induced h11  pointer */,
  const int expansion_type/* 1 double fourier, 0: spherical harmonic */
  );
  
void
bbn_inclusion_map_S2_to_M_CS
  (
  Grid_T *const grid/* grid */,
  const char *const type/* NS or BH */,
  const unsigned Ntheta/* number of points in theta direction */,
  const unsigned Nphi/* number of points in phi direction */,
  const unsigned lmax/* l max in Ylm */,
  const int expansion_type/* 1 double fourier, 0: spherical harmonic */,
  const double *const S2akv_U0/* akv on S2 */,
  const double *const S2akv_U1/* akv on S2 */,
  const char *const name_akv_U0/* inclusion akv vector v^0 name */,
  const char *const name_akv_U1/* inclusion akv vector v^1 name */,
  const char *const name_akv_U2/* inclusion akv vector v^2 name */
  );

void bbn_populate_spin_integrands_akv(Patch_T *const patch,const double *const normal[3]);

void 
bbn_define_spin_akv
  (
  double S[3]/* spin Sx,Sy,Sz */,
  Grid_T *const grid/* grid */,
  const char *const kind/* "NS" or "BH" */
  );



#endif






