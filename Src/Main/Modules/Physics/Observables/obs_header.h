#ifndef obs_LIB_H
#define obs_LIB_H

#include "core_lib.h"
#include "physics_observables_lib.h"
#include "utilities_lib.h"
#include "manifold_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "fields_lib.h"
#include "managers_lib.h"

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

double obs_ADM_mass(Observable_T *const obs);
double obs_BH_ADM_mass(Observable_T *const obs);
double obs_Kommar_mass(Observable_T *const obs);
double obs_ADM_mass(Observable_T *const obs);
void obs_populate_spin_integrands_akv(Patch_T *const patch,const double *const normal[3]);
void obs_populate_spin_integrands_Campanelli(Patch_T *const patch,const double xc[3],const double *const normal[3]);
void obs_Rc_NS(double Rc[3],Obj_Man_T *const obj);
void obs_populate_ADM_integrand_PdS_GdV_binary(const Observable_T *const obs);
void obs_populate_ADM_integrand_PdS_GdV_single(const Observable_T *const obs);

#endif

