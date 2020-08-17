/*
// Alireza Rashti
// August 2020
*/

#include "spherical_harmonics.h"

/* ->: d(P^m_l(x))/dtheta. */
double 
d_associated_legendre_dtheta
  (
  const int l/* l in P^l_m(x) */,
  const int m/* m in P^l_m(x), 0 <= m <= l*/,
  const double x/* x in P^l_m(x), -1 <= x=cos(theta) <= 1 */
  )
{
  double dp = 0;
  
  if (!l)
  {
    return 0;
  }
  else if (m == l)
  {
    dp = (l+m)*(l-m+1)*associated_legendre(l,m-1,x);
  }
  else if (m == 0)
  {
    dp = -(1+(l+m)*(l-m+1)/(l*(l+1)))*
          associated_legendre(l,1,x);
  }
  else
  {
    dp = (l+m)*(l-m+1)*associated_legendre(l,m-1,x)
          -associated_legendre(l,m+1,x);
  }
  
  return -dp*0.5;
}

/* ->: associated Legendre polynomial P^m_l(x). 
// this recurrence relation claimed (numerical recipe) to be stable */
double 
associated_legendre
  (
  const int l/* l in P^l_m(x) */,
  const int m/* m in P^l_m(x), 0 <= m <= l*/,
  const double x/* x in P^l_m(x), -1 <= x=cos(theta) <= 1 */
  )
{
  double fact, pll = 0, pmm, pmmp1, somx2;
  int i, ll;

  if (x > 1. || x < -1.)
    Error0("x exceeds from [-1,1].\n");
    
  if(m < 0 || m > l) 
    Error0("m exceeds from [-l,l].\n");

  /* compute P^m_m */
  pmm = 1.0;
  if(m >= 1)
  {
    somx2 = sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for(i=1; i<=m; ++i)
    {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }/* if(m > 0)  */
  
  if(l == m) 
  {
    return pmm;
  }/* if(l == m)  */
  else
  {
    pmmp1 = x*(2*m+1)*pmm;
    
    if(l == (m+1))
    { 
      return pmmp1;
    }
    else
    {
      for(ll=m+2; ll<=l; ++ll)
      {
        pll    = (x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm)/(ll-m);
        pmm    = pmmp1;
        pmmp1 = pll;
      } /* now k=l */
      return pll;
    }
  }
  return DBL_MAX;/* must not reach here */
}

/* ->: Y_l^m(theta,phi) */
double complex 
Ylm
  (
  const int l/* l in Y_l^m(theta,phi) */, 
  const int m/* m in Y_l^m(theta,phi), -l <= m <= l */, 
  const double theta/* theta in Y_l^m(theta,phi) */,
  const double phi/*  phi in Y_l^m(theta,phi) */
  )
{
  if (theta > M_PI || theta < 0)
    Error0("theta exceeds from [0,pi] interval.\n");
  if (phi > 2*M_PI || phi < 0)
    Error0("phi exceeds from [0,2*pi] interval.\n");
  
  const int pm       = m >= 0 ? m : -m;/* positive m */
  const int two_pm   = 2*pm;
  const int lpm      = l+pm;
  double factor      = 1;
  double complex ylm = 0;
  int i;
  
  /* calculate the normalization factor */
  for (i = 0; i < two_pm ; ++i)
  {
    factor /= sqrt(lpm-i);
  }
  factor *= 0.5*sqrt((2*l+1)/M_PI);
  
  if (m >= 0)
  {
    ylm = factor*
          associated_legendre(l,m,cos(theta))*
          cexp(I*(double)m*phi);
  }
  else
  {
    double SIGN[2] = {1,-1};
    int mp = -m;
    ylm = SIGN[mp%2]*factor*
          associated_legendre(l,mp,cos(theta))*
          cexp(I*(double)m*phi);
  }
  
  return ylm;
}

/* ->: dY_l^m(theta,phi)/dphi */
double complex 
dYlm_dphi
  (
  const int l/* l in Y_l^m(theta,phi) */, 
  const int m/* m in Y_l^m(theta,phi), -l <= m <= l */, 
  const double theta/* theta in Y_l^m(theta,phi) */,
  const double phi/*  phi in Y_l^m(theta,phi) */
  )
{
  return (I*(double)m*Ylm(l,m,theta,phi));
}

/* ->: dY_l^m(theta,phi)/dtheta */
double complex 
dYlm_dtheta
  (
  const int l/* l in Y_l^m(theta,phi) */, 
  const int m/* m in Y_l^m(theta,phi), -l <= m <= l */, 
  const double theta/* theta in Y_l^m(theta,phi) */,
  const double phi/*  phi in Y_l^m(theta,phi) */
  )
{
  if (theta > M_PI || theta < 0)
    Error0("theta exceeds from [0,pi] interval.\n");
  if (phi > 2*M_PI || phi < 0)
    Error0("phi exceeds from [0,2*pi] interval.\n");
  
  const int pm        = m >= 0 ? m : -m;/* positive m */
  const int two_pm    = 2*pm;
  const int lpm       = l+pm;
  double factor       = 1;
  double complex dylm = 0;
  const double x      = cos(theta);
  int i;
  
  /* calculate the normalization factor */
  for (i = 0; i < two_pm ; ++i)
  {
    factor /= sqrt(lpm-i);
  }
  factor *= 0.5*sqrt((2*l+1)/M_PI);
  
  if (m >= 0)
  {
    dylm = factor*
          d_associated_legendre_dtheta(l,m,x)*
          cexp(I*(double)m*phi);
  }
  else
  {
    double SIGN[2] = {1,-1};
    int mp = -m;
    dylm = SIGN[mp%2]*factor*
          d_associated_legendre_dtheta(l,mp,x)*
          cexp(I*(double)m*phi);
  }
  
  return dylm;
}

