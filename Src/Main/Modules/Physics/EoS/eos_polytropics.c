/*
// Alireza Rashti
// June 2019
*/

/* used article:
//  @inproceedings
//  {Smith2013TolmanOppenheimerVolkoffS,
//  title={Tolman–Oppenheimer–Volkoff (TOV) Stars},
//  author={A. Smith},
//  year={2013}
//  } */

#include "eos_polytropics.h"

/* calculate pressure in terms of h for pwp
// ->return value: p(h) */
double EoS_p_h_pwp(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const Uint i = find_threshold_number_h(eos);
  const double K = eos->K[i];
  const double h = eos->h;
  const double a = eos->a[i];
  const double n = eos->n[i];
  
  return pow((h-1-a)/(n+1),n+1)*pow(K,-n);
}

/* calculate rest mass density in terms of h for pwp
// ->return value: rho(h) */
double EoS_rho0_h_pwp(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const Uint i = find_threshold_number_h(eos);    
  const double K = eos->K[i];
  const double h = eos->h;
  const double a = eos->a[i];
  const double n = eos->n[i];
  
  return pow((h-1-a)/(n+1),n)*pow(K,-n);
}

/* calculate d(rest mass density)/dh in terms of h for pwp
// ->return value: d(rho(h))/dh */
double EoS_drho0_dh_h_pwp(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const Uint i = find_threshold_number_h(eos);    
  const double K = eos->K[i];
  const double h = eos->h;
  const double a = eos->a[i];
  const double n = eos->n[i];
  
  return pow((h-1-a)/(n+1),n-1)*pow(K,-n)*n/(n+1);
}

/* calculate d(rest mass density)/dh in terms of h for polytrop
// ->return value: d(rho(h))/dh */
double EoS_drho0_dh_h_p(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const double K = eos->K[0];
  const double h = eos->h;
  const double n = eos->n[0];
  
  return pow((h-1)/(n+1),n-1)*pow(K,-n)*n/(n+1);
}

/* calculate the total energy density in terms of h for pwp
// ->return value: e(h) */
double EoS_e_h_pwp(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const double e0 = EoS_e0_h_pwp(eos);
  
  return EoS_rho0_h_pwp(eos)*(1. + e0);
}

/* calculate the specific internal energy in terms of h for pwp
// ->return value: e0(h) */
double EoS_e0_h_pwp(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const Uint i = find_threshold_number_h(eos);
  const double h = eos->h;
  const double a = eos->a[i];
  const double n = eos->n[i];
  
  return (a+n*(h-1))/(n+1);
}

/* calculate d(total energy density)/dh in terms of h for pwp
// ->return value: de(h)/dh */
double EoS_de_dh_h_pwp(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const Uint i = find_threshold_number_h(eos);
  const double h = eos->h;
  const double a = eos->a[i];
  const double n = eos->n[i];
  
  return EoS_drho0_dh_h_pwp(eos)*(1+(a+n*(h-1))/(n+1))+EoS_rho0_h_pwp(eos)*n/(n+1);
}

/* calculate d(total energy density)/dh in terms of h for polytropic
// ->return value: de(h)/dh */
double EoS_de_dh_h_p(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const double h = eos->h;
  const double n = eos->n[0];
  
  return EoS_drho0_dh_h_p(eos)*(1+n*(h-1)/(n+1))+EoS_rho0_h_p(eos)*n/(n+1);
}

/* calculate pressure in terms of h for polytropic
// ->return value: p(h) */
double EoS_p_h_p(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const double K = eos->K[0];
  const double h = eos->h;
  const double n = eos->n[0];
  
  return pow((h-1)/(n+1),n+1)*pow(K,-n);
}

/* calculate rest mass density in terms of h for polytropic
// ->return value: rho(h) */
double EoS_rho0_h_p(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const double K = eos->K[0];
  const double h = eos->h;
  const double n = eos->n[0];
  
  return pow((h-1)/(n+1),n)*pow(K,-n);
}


/* calculate the total energy density in terms of h for polytropic
// ->return value: e(h) */
double EoS_e_h_p(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const double e0 = EoS_e0_h_p(eos);
  
  return EoS_rho0_h_p(eos)*(1. + e0);
}

/* calculate the specific internal energy in terms of h for polytropic
// ->return value: e0(h) */
double EoS_e0_h_p(EoS_T *const eos)
{
  if (EQL(eos->h,1)) eos->h = 1;
  
  const double h = eos->h;
  const double n = eos->n[0];
  
  return n*(h-1)/(n+1);
}

/* given eos->h it finds this enthalpy takes place in which intervale.
// if it could not find this it gives error.
// note: since it is required that rho0_th be written increasingly
// then h is also comes out increasingly thus:
// 1 < h[0] < h[1] < h[2] < ... < h[n]
// ->return value: the threshold number which h is falls in. */
static Uint find_threshold_number_h(const EoS_T *const eos)
{
  Uint i;
  char str[1000];
  Flag_T flg = NONE;
  
  if (eos->enthalpy_fatal && LSS(eos->h,1))
    Error0("The value of the enthalpy is not set correctly.\n");
  
  if (LSSEQL(eos->h,eos->h_th[1]))
    return 0;
  else if (GRTEQL(eos->h,eos->h_th[eos->N-1]))
    return eos->N-1;
  
  for (i = 1; i < eos->N-1; ++i)
  {
    if (GRTEQL(eos->h,eos->h_th[i]) &&
        LSSEQL(eos->h,eos->h_th[i+1]) )
        {
          flg = FOUND;
          break;
        }
  }
  
  if (flg == NONE)
  {
    sprintf(str,"%e",eos->h);
    Errors("Threshold number could not be found for enthalpy = %s.\n",str);
  }
  
  return i;
}

/* calculate pressure in terms of h for pwp and natural cubic spline
// ->return value: p(h) */
double EoS_p_h_pwp_ncs(EoS_T *const eos)
{
  if (LSSEQL(eos->h, eos->cubic_spline->h_floor))
    return 0.;

  double p;  
  Interpolation_T *const interp_s = eos->cubic_spline->interp_p;
  
  interp_s->N_cubic_spline_1d->h  = eos->h;
  p = execute_interpolation(interp_s);
  
  return (LSSEQL(p,0.) || p == DBL_MAX ? 0. : p);
}

/* calculate rest mass density in terms of h for pwp and natural cubic spline
// ->return value: rho(h) */
double EoS_rho0_h_pwp_ncs(EoS_T *const eos)
{
  if (LSSEQL(eos->h, eos->cubic_spline->h_floor))
    return 0.;
  
  double rho0;  
  Interpolation_T *const interp_s = eos->cubic_spline->interp_rho0;
  
  interp_s->N_cubic_spline_1d->h  = eos->h;
  rho0 = execute_interpolation(interp_s);
  
  return (LSSEQL(rho0,0.) || rho0 == DBL_MAX ? 0. : rho0);
}

/* calculate the total energy density in terms of h for pwp and natural cubic spline
// ->return value: e(h) */
double EoS_e_h_pwp_ncs(EoS_T *const eos)
{
  if (LSSEQL(eos->h, eos->cubic_spline->h_floor))
    return 0.;

  double e;
  Interpolation_T *const interp_s = eos->cubic_spline->interp_e;
  
  interp_s->N_cubic_spline_1d->h  = eos->h;
  e = execute_interpolation(interp_s);
  
  return (LSSEQL(e,0.) || e == DBL_MAX ? 0. : e);
}
