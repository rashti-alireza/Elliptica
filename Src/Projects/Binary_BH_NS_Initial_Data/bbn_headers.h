#ifndef bbn_headers_LIB_H
#define bbn_headers_LIB_H



#include "core_lib.h"
#include "TOV_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_calculus_lib.h"
#include "physics_EoS_lib.h"
#include "physics_stress_energy_lib.h"
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
  const double *V2CM;
  Uint find_y_CM: 1;
  Uint find_x_CM: 1;
  Uint find_Omega: 1;
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
  Uint surface_integration_flg: 1;/* if 1 means it measn 
                                      // we need surface integration 
                                      // on this patch as well, 
                                      // 0 means, no need */
  /* which hypersurface the surface integral is carried out */
  Uint X_surface: 1;
  Uint Y_surface: 1;
  Uint Z_surface: 1;
  /* index of hypersurface for each X,Y and Z respectively */
  Uint I;
  Uint J;
  Uint K;
};

/* alpha and beta (gauges) for exporting of ID */
struct IDGauge_S
{
  Grid_T *grid;
  const char *lapse_type;/* o. XCTS: whatever made by XCTS eqs.
                         // o. one: set all to 1.
                         // o. puncture1: alpha = 1/(bssn_psi)^2. 
                         // o. puncture2: alpha = 1/(1+bssn_psi^4). */
  const char *shift_type;/* o. XCTS: whatever made by XCTS eqs.
                         // o. zero: set all to 0. */
  double Mb;/* bare mass */
  double r_CutOff;/* smalles R for puncture */
  double rfill;/* where the BH filling started, generally BH radius */
  double rmin;/* how much inside the BH, rmin <= rfill */
  /* puncture like for psi ~ 1+Mb/(r+r_CutOff) */
  double (*psi_punc0)(const double r, const double Mb,const double r_CutOff);
};

void bbn_free_data_Gamma_patch(Patch_T *const patch);
void bbn_study_initial_data(Grid_T *const grid);
void bbn_print_fields(Grid_T *const grid,const Uint iteration, const char *const folder);
void bbn_update_psi10A_UiUj(Patch_T *const patch);
void bbn_populate_free_data(Grid_T *const grid);
void bbn_free_data_gammas(Grid_T *const grid);
void bbn_rm_1st_derivatives_conformal_metric(Patch_T *const patch);
void bbn_1st_derivatives_conformal_metric(Patch_T *const patch);
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
void bbn_print_properties(Grid_T *const grid,const Uint iteration, const char *const folder,const char *const open_file_mode,const int pr_flg);
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
void bbn_create_grid_prototype_BC(Grid_T *const grid);
void bbn_1st_2nd_derivatives_conformal_metric(Patch_T *const patch);
void bbn_rm_1st_2nd_derivatives_conformal_metric(Patch_T *const patch);
void bbn_add_and_take_2nd_derivatives_K(Patch_T *const patch);

int 
bbn_bhfiller
  (
  Grid_T *const grid/* the whole grid */,
  const char *const method/* the method to be used for extrapolating */
  );

double bbn_bhf_smoother(const double r, const double rmax,const double rmin);
void bbn_bhf_ChebTn_extrapolate (double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill,const Uint N);
void bbn_rm_1st_2nd_derivatives_Kij(Patch_T *const patch);
void bbn_1st_2nd_derivatives_adm_Kij(Patch_T *const patch);
void bbn_adm_Kij(Patch_T *const patch);
void bbn_populate_spin_integrands_akv(Patch_T *const patch,const double *const normal[3]);

void 
bbn_define_spin_akv
  (
  double S[3]/* spin Sx,Sy,Sz */,
  Grid_T *const grid/* grid */,
  const char *const kind/* "NS" or "BH" */
  );

void bbn_free_data_g_gI_analytic(
        Patch_T *const patch,
        double *(*get_v)(const char *const fname,void *params),
        void *params);
void bbn_free_data_dg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params);
void bbn_free_data_ddg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params);
void bbn_free_data_dddg_analytic(
	Patch_T *const patch, 
	double *(*get_v)(const char *const fname,void *params),
	void *params);
	
void bbn_ks_free_data_set_params(Grid_T *const grid);
double *bbn_ks_read_analytic(const char *const name, void *params);
void bbn_bam_set_gauges(struct IDGauge_S *const gauge);
void bbn_free_date_dGammaConf(Patch_T *const patch);


#endif






