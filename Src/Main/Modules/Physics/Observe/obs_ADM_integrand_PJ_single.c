/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "core_lib.h"
#include "physics_observe_lib.h"
#include "utilities_lib.h"
#include "manifold_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "fields_lib.h"
#include "obs_header.h"

#define declare_and_alloc_xi(name) \
  double *name = alloc_double(nn);

#define add_and_get_field(name) \
  if (_Ind(#name) >= 0)\
  {DECLARE_FIELD(name);REMOVE_FIELD(name);}\
  ADD_FIELD(name);REALLOC_v_WRITE_v(name);


void obs_populate_ADM_integrand_PdS_GdV_single(const Observe_T *const obs);
void obs_populate_ADM_integrand_PdS_GdV_single(const Observe_T *const obs)
{
  Physics_T *const phys = obs->phys;
  struct items_S **adm = obs->items;
  const unsigned N = obs->Nitems;
  const double x_cm = Getd("x_CM");
  const double y_cm = Getd("y_CM");
  const double z_cm = Getd("z_CM");
  unsigned p;

  for(p = 0; p < N; ++p)
  {
  Patch_T *patch = adm[p]->patch;
  unsigned nn = patch->nn;
  unsigned ijk;

  /* declaring: */
  READ_v(_A_UiUj_U2U2)
  READ_v(_A_UiUj_U1U2)
  READ_v(_A_UiUj_U1U1)
  READ_v(_A_UiUj_U0U2)
  READ_v(_A_UiUj_U0U1)
  READ_v(_A_UiUj_U0U0)
  READ_v(gConf_D0D2)
  READ_v(gConf_D0D0)
  READ_v(gConf_D0D1)
  READ_v(gConf_D1D2)
  READ_v(gConf_D1D1)
  READ_v(gConf_D2D2)
  READ_v(igConf_U2U2)
  READ_v(igConf_U1U2)
  READ_v(igConf_U1U1)
  READ_v(igConf_U0U2)
  READ_v(igConf_U0U0)
  READ_v(igConf_U0U1)
  READ_v(psi)
  READ_v(K)
  add_and_get_field(ADM_integrand_P_U0)
  add_and_get_field(ADM_integrand_P_U1)
  add_and_get_field(ADM_integrand_P_U2)
  add_and_get_field(ADM_integrand_xiP_U2)
  add_and_get_field(ADM_integrand_xiP_U1)
  add_and_get_field(ADM_integrand_xiP_U0)
  declare_and_alloc_xi(xi_U2)
  declare_and_alloc_xi(xi_U0)
  declare_and_alloc_xi(xi_U1)



   for(ijk = 0; ijk < nn; ++ijk)
   {
   double x    = patch->node[ijk]->x[0];
   double y    = patch->node[ijk]->x[1];
   double z    = patch->node[ijk]->x[2];
   xi_U0[ijk] = x-x_cm;
   xi_U1[ijk] = y-y_cm;
   xi_U2[ijk] = z-z_cm;
   }
      const double *n_U0 = adm[p]->n_U0;
      const double *n_U1 = adm[p]->n_U1;
      const double *n_U2 = adm[p]->n_U2;
      for (ijk = 0; ijk < nn; ++ijk)
      {
      double psi4 = 
pow(psi[ijk], 4);

      double psi6 = 
pow(psi[ijk], 6);

      double P_U0U1 = 
-0.66666666666666663*K[ijk]*igConf_U0U1[ijk] + _A_UiUj_U0U1[ijk]/
psi6;

      double P_U0U0 = 
-0.66666666666666663*K[ijk]*igConf_U0U0[ijk] + _A_UiUj_U0U0[ijk]/
psi6;

      double P_U0U2 = 
-0.66666666666666663*K[ijk]*igConf_U0U2[ijk] + _A_UiUj_U0U2[ijk]/
psi6;

      double P_U2U2 = 
-0.66666666666666663*K[ijk]*igConf_U2U2[ijk] + _A_UiUj_U2U2[ijk]/
psi6;

      double P_U1U1 = 
-0.66666666666666663*K[ijk]*igConf_U1U1[ijk] + _A_UiUj_U1U1[ijk]/
psi6;

      double P_U1U2 = 
-0.66666666666666663*K[ijk]*igConf_U1U2[ijk] + _A_UiUj_U1U2[ijk]/
psi6;

      double Pn_U2 = 
psi4*(P_U0U0*gConf_D0D0[ijk]*gConf_D0D2[ijk]*n_U0[ijk] + P_U0U0*
gConf_D0D1[ijk]*gConf_D0D2[ijk]*n_U1[ijk] + P_U0U0*
pow(gConf_D0D2[ijk], 2)*n_U2[ijk] + P_U0U1*gConf_D0D0[ijk]*
gConf_D1D2[ijk]*n_U0[ijk] + P_U0U1*gConf_D0D1[ijk]*gConf_D0D2[ijk]*
n_U0[ijk] + P_U0U1*gConf_D0D1[ijk]*gConf_D1D2[ijk]*n_U1[ijk] + P_U0U1*
gConf_D0D2[ijk]*gConf_D1D1[ijk]*n_U1[ijk] + 2.0*P_U0U1*gConf_D0D2[ijk]*
gConf_D1D2[ijk]*n_U2[ijk] + P_U0U2*gConf_D0D0[ijk]*gConf_D2D2[ijk]*
n_U0[ijk] + P_U0U2*gConf_D0D1[ijk]*gConf_D2D2[ijk]*n_U1[ijk] + P_U0U2*
pow(gConf_D0D2[ijk], 2)*n_U0[ijk] + P_U0U2*gConf_D0D2[ijk]*
gConf_D1D2[ijk]*n_U1[ijk] + 2.0*P_U0U2*gConf_D0D2[ijk]*gConf_D2D2[ijk]*
n_U2[ijk] + P_U1U1*gConf_D0D1[ijk]*gConf_D1D2[ijk]*n_U0[ijk] + P_U1U1*
gConf_D1D1[ijk]*gConf_D1D2[ijk]*n_U1[ijk] + P_U1U1*
pow(gConf_D1D2[ijk], 2)*n_U2[ijk] + P_U1U2*gConf_D0D1[ijk]*
gConf_D2D2[ijk]*n_U0[ijk] + P_U1U2*gConf_D0D2[ijk]*gConf_D1D2[ijk]*
n_U0[ijk] + P_U1U2*gConf_D1D1[ijk]*gConf_D2D2[ijk]*n_U1[ijk] + P_U1U2*
pow(gConf_D1D2[ijk], 2)*n_U1[ijk] + 2.0*P_U1U2*gConf_D1D2[ijk]*
gConf_D2D2[ijk]*n_U2[ijk] + P_U2U2*gConf_D0D2[ijk]*gConf_D2D2[ijk]*
n_U0[ijk] + P_U2U2*gConf_D1D2[ijk]*gConf_D2D2[ijk]*n_U1[ijk] + P_U2U2*
pow(gConf_D2D2[ijk], 2)*n_U2[ijk]);

      double Pn_U1 = 
psi4*(P_U0U0*gConf_D0D0[ijk]*gConf_D0D1[ijk]*n_U0[ijk] + P_U0U0*
pow(gConf_D0D1[ijk], 2)*n_U1[ijk] + P_U0U0*gConf_D0D1[ijk]*
gConf_D0D2[ijk]*n_U2[ijk] + P_U0U1*gConf_D0D0[ijk]*gConf_D1D1[ijk]*
n_U0[ijk] + P_U0U1*pow(gConf_D0D1[ijk], 2)*n_U0[ijk] + 2.0*P_U0U1*
gConf_D0D1[ijk]*gConf_D1D1[ijk]*n_U1[ijk] + P_U0U1*gConf_D0D1[ijk]*
gConf_D1D2[ijk]*n_U2[ijk] + P_U0U1*gConf_D0D2[ijk]*gConf_D1D1[ijk]*
n_U2[ijk] + P_U0U2*gConf_D0D0[ijk]*gConf_D1D2[ijk]*n_U0[ijk] + P_U0U2*
gConf_D0D1[ijk]*gConf_D0D2[ijk]*n_U0[ijk] + 2.0*P_U0U2*gConf_D0D1[ijk]*
gConf_D1D2[ijk]*n_U1[ijk] + P_U0U2*gConf_D0D1[ijk]*gConf_D2D2[ijk]*
n_U2[ijk] + P_U0U2*gConf_D0D2[ijk]*gConf_D1D2[ijk]*n_U2[ijk] + P_U1U1*
gConf_D0D1[ijk]*gConf_D1D1[ijk]*n_U0[ijk] + P_U1U1*
pow(gConf_D1D1[ijk], 2)*n_U1[ijk] + P_U1U1*gConf_D1D1[ijk]*
gConf_D1D2[ijk]*n_U2[ijk] + P_U1U2*gConf_D0D1[ijk]*gConf_D1D2[ijk]*
n_U0[ijk] + P_U1U2*gConf_D0D2[ijk]*gConf_D1D1[ijk]*n_U0[ijk] + 2.0*
P_U1U2*gConf_D1D1[ijk]*gConf_D1D2[ijk]*n_U1[ijk] + P_U1U2*
gConf_D1D1[ijk]*gConf_D2D2[ijk]*n_U2[ijk] + P_U1U2*
pow(gConf_D1D2[ijk], 2)*n_U2[ijk] + P_U2U2*gConf_D0D2[ijk]*
gConf_D1D2[ijk]*n_U0[ijk] + P_U2U2*pow(gConf_D1D2[ijk], 2)*n_U1[ijk] + 
P_U2U2*gConf_D1D2[ijk]*gConf_D2D2[ijk]*n_U2[ijk]);

      double Pn_U0 = 
psi4*(P_U0U0*pow(gConf_D0D0[ijk], 2)*n_U0[ijk] + P_U0U0*
gConf_D0D0[ijk]*gConf_D0D1[ijk]*n_U1[ijk] + P_U0U0*gConf_D0D0[ijk]*
gConf_D0D2[ijk]*n_U2[ijk] + 2.0*P_U0U1*gConf_D0D0[ijk]*gConf_D0D1[ijk]*
n_U0[ijk] + P_U0U1*gConf_D0D0[ijk]*gConf_D1D1[ijk]*n_U1[ijk] + P_U0U1*
gConf_D0D0[ijk]*gConf_D1D2[ijk]*n_U2[ijk] + P_U0U1*
pow(gConf_D0D1[ijk], 2)*n_U1[ijk] + P_U0U1*gConf_D0D1[ijk]*
gConf_D0D2[ijk]*n_U2[ijk] + 2.0*P_U0U2*gConf_D0D0[ijk]*gConf_D0D2[ijk]*
n_U0[ijk] + P_U0U2*gConf_D0D0[ijk]*gConf_D1D2[ijk]*n_U1[ijk] + P_U0U2*
gConf_D0D0[ijk]*gConf_D2D2[ijk]*n_U2[ijk] + P_U0U2*gConf_D0D1[ijk]*
gConf_D0D2[ijk]*n_U1[ijk] + P_U0U2*pow(gConf_D0D2[ijk], 2)*n_U2[ijk] + 
P_U1U1*pow(gConf_D0D1[ijk], 2)*n_U0[ijk] + P_U1U1*gConf_D0D1[ijk]*
gConf_D1D1[ijk]*n_U1[ijk] + P_U1U1*gConf_D0D1[ijk]*gConf_D1D2[ijk]*
n_U2[ijk] + 2.0*P_U1U2*gConf_D0D1[ijk]*gConf_D0D2[ijk]*n_U0[ijk] + 
P_U1U2*gConf_D0D1[ijk]*gConf_D1D2[ijk]*n_U1[ijk] + P_U1U2*
gConf_D0D1[ijk]*gConf_D2D2[ijk]*n_U2[ijk] + P_U1U2*gConf_D0D2[ijk]*
gConf_D1D1[ijk]*n_U1[ijk] + P_U1U2*gConf_D0D2[ijk]*gConf_D1D2[ijk]*
n_U2[ijk] + P_U2U2*pow(gConf_D0D2[ijk], 2)*n_U0[ijk] + P_U2U2*
gConf_D0D2[ijk]*gConf_D1D2[ijk]*n_U1[ijk] + P_U2U2*gConf_D0D2[ijk]*
gConf_D2D2[ijk]*n_U2[ijk]);

      double xiP_U2 = 
-Pn_U0*xi_U1[ijk] + Pn_U1*xi_U0[ijk];

      double xiP_U1 = 
Pn_U0*xi_U2[ijk] - Pn_U2*xi_U0[ijk];

      double xiP_U0 = 
-Pn_U1*xi_U2[ijk] + Pn_U2*xi_U1[ijk];


      /* populating: */
      ADM_integrand_P_U0[ijk] = Pn_U0;
      ADM_integrand_P_U1[ijk] = Pn_U1;
      ADM_integrand_P_U2[ijk] = Pn_U2;

      /* populating: */
      ADM_integrand_xiP_U2[ijk] = xiP_U2;
      ADM_integrand_xiP_U1[ijk] = xiP_U1;
      ADM_integrand_xiP_U0[ijk] = xiP_U0;
      }

  free(xi_U0);
  free(xi_U1);
  free(xi_U2);
  }
}