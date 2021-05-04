#ifndef obs_LIB_H
#define obs_LIB_H

#include "core_lib.h"
#include "physics_observe_lib.h"
#include "utilities_lib.h"
#include "manifold_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "fields_lib.h"
#include "physics_lib.h"

/* prefix of parameters for this project */
#define P_ "observe_"


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
  /* volume integration type */
  const char *vol_integration_type;/* ["Integral{f(x)dV}[i,f],Spectral",
                                      "Integral{f(x)dV},Spectral"] */
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
  /* interval for volume integration */
  Uint Ii,If;
  Uint Ji,Jf;
  Uint Ki,Kf;
};

double obs_ADM_mass_SV_isotropic(Observe_T *const obs);
double obs_ADM_mass_S2(Observe_T *const obs);
double obs_Komar_mass(Observe_T *const obs);
double obs_ADM_mass(Observe_T *const obs);
void obs_populate_spin_integrands_akv(Patch_T *const patch,const double *const normal[3]);
void obs_populate_spin_integrands_Campanelli(Patch_T *const patch,const double xc[3],const double *const normal[3]);
void obs_Rc_NS(Observe_T *const obs);
void obs_BH_irreducible_mass_CS(Observe_T *const obs);
void obs_ADM_P_Stokes_SV_Ossokine(Observe_T *const obs);
void obs_ADM_J_Stokes_SV_Ossokine(Observe_T *const obs);
void obs_ADM_P_Stokes_SV_Rashti(Observe_T *const obs);
void obs_ADM_P_S_default(Observe_T *const obs);
void obs_ADM_J_S_default(Observe_T *const obs);
void obs_ADM_P_Stokes_SV_constraint(Observe_T *const obs);
void obs_ADM_J_Stokes_SV_constraint(Observe_T *const obs);
double obs_ADM_mass_SV_conformal(Observe_T *const obs);
double obs_integral_SV(Observe_T *const obs,
                        const char *const sS/* integrand for S */,
                        const char *const sV/* intergrand for V */,
                        const char sign_sS/* [+/-] integral of S */,
                        const char sign_sV/* [+/-] integral of V */);
#endif

