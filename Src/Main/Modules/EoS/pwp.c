/*
// Alireza Rashti
// June 2019
*/

#include "pwp.h"

/* calculate pressure in terms of h for pwp
// ->return value: p(h) */
double EoS_p_h_pwp(EoS_T *const eos)
{
  const unsigned i = find_threshold_number_h(eos);
  const double K = eos->K[i];
  const double h = eos->h;
  const double a = eos->a[i];
  const double n = eos->n[i];
  
  return K*pow((h-1-a)/(K*(n+1)),n+1);
}

/* calculate rest mass density in terms of h for pwp
// ->return value: rho(h) */
double EoS_rho_h_pwp(EoS_T *const eos)
{
  const unsigned i = find_threshold_number_h(eos);    
  const double K = eos->K[i];
  const double h = eos->h;
  const double a = eos->a[i];
  const double n = eos->n[i];
  
  return pow((h-1-a)/(K*(n+1)),n);
}

/* calculate the total energy density in terms of h for pwd
// ->return value: e(h) */
double EoS_e_h_pwp(EoS_T *const eos)
{
  const unsigned i = find_threshold_number_h(eos);
  const double h = eos->h;
  const double a = eos->a[i];
  const double n = eos->n[i];
  
  return EoS_rho_h_pwp(eos)*(1+(a+n*(h-1))/(n+1));
}

/* given eos->h it finds this enthalpy takes place in which intervale.
// if it could not find this it gives error.
// note: since it is required that rho_th be written increasingly
// then h is also comes out increasingly thus:
// 1 < h[0] < h[1] < h[2] < ... < h[n]
// ->return value: the threshold number which h is falls in. */
static unsigned find_threshold_number_h(const EoS_T *const eos)
{
  unsigned i;
  Flag_T flg = NONE;
  
  if (LSS(eos->h,1))
    abortEr("The value of the enthalpy is not set correctly.\n");
  
  if (LSSEQL(eos->h,eos->h_th[0]))
    return 0;
  else if (GRTEQL(eos->h,eos->h_th[eos->N-1]))
    return eos->N;
  
  for (i = 1; i < eos->N-1;++i)
  {
    if (GRTEQL(eos->h,eos->h_th[i-1]) &&
        LSSEQL(eos->h,eos->h_th[i]) )
        {
          flg = FOUND;
          break;
        }
  }
  
  if (flg == NONE)
    abortEr("Threshold number could not be found.\n");
  
  return i;
}
