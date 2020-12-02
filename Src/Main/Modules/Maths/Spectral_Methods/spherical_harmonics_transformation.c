/*
// Alireza Rashti
//  August 2019
*/

#include "spherical_harmonics_transformation.h"

/* calculatinf coefficents of spherical harmonic Clm for:
// f(theta,phi) = \sum_{l=0}^{l=lmax} \sum_{m = -l}^{l} C_{lm}*Y_{lm}(thata,phi)
// Y_{lm}=\left( -1\right) ^{m}\sqrt {\frac {2l+1} {4\pi }\frac {\left( l-m\right) !} 
//         {\left( l+m\right) !}}P_{1}^{m}\left( \cos \theta \right) e^{im\varphi }
// 
// some notes:
// ===========
// o. l is inclusive
// o. f(theta, phi) is given in a 1d array with this mapping:
//    f(i,j) = f[j+i*Nphi] = f[IJ_Ylm(i,j,Nphi)]
// o. It's the best to set Ntheta = Ntheta_Ylm(lmax) and
//    Nphi = Nphi_Ylm(lmax).
// o. Clm = realC[lm2n(l,m,lmax)] + imagC[lm2n(l,m,lmax)]*I;
// o. Clm is only calculated for m >= 0, since C(l,-m) = (-1)^m C*(l,m)
// o. since m >= 0 the total memory is N = (Lmax+1)*Lmax/2 + Lmax+1
// o. an interesting observation is that even some simple functions like:
//    sin(theta)*(cos(phi)-*sin(2*phi)) requires very large Lmax, bigger than 20!
*/
void 
get_Ylm_coeffs(
double *const realClm/* Re(C_{lm}) coeffs */,
double *const imagClm/* Im(C_{lm}) coeffs */,
const double *const f/* f(theta,phi) */,
const Uint Ntheta/* Ntheta points in theta direction */,
const Uint Nphi/* Nphi points in phi direction*/,
const Uint Lmax/* maximum l (inclusive) for the expansion */
)
{
  const char *collocation = "GaussLegendre_EquiSpaced";
  /* Note: later on you can add parameter to direct the collocation */
  if (strcmp_i(collocation,"GaussLegendre_EquiSpaced"))
    get_Ylm_coeffs_GaussLegendre_EquiSpaced(realClm,imagClm,f,Ntheta,Nphi,Lmax);
  else
    Error0(NO_JOB);

}

/* calculating complex coeffs of Ylm in expansion of function f(theta,phi)
// in calculation it assumes:
//    phi direction collocation points are Equispaced and
//    theta direction collocation points are Chebyshev extrema. */
static void 
get_Ylm_coeffs_GaussLegendre_EquiSpaced(
double *const realClm/* Re(C_{lm}) coeffs */,
double *const imagClm/* Im(C_{lm}) coeffs */,
const double *const f/* f(theta,phi) */,
const Uint Ntheta/* Ntheta points in theta direction */,
const Uint Nphi/* Nphi points in phi direction*/,
const Uint Lmax/* maximum l (inclusive) for the expansion */)
{
  double **Ftheta    = calloc(Ntheta,sizeof(*Ftheta));IsNull(Ftheta);
  double *real_fYlm  = alloc_double(Ntheta);/* Re( \int f(theta,phi) Ylm^* sin(theta) exp^{-i m phi} dphi) */
  double *imag_fYlm  = alloc_double(Ntheta);/* Im( \int f(theta,phi) Ylm^* sin(theta) exp^{-i m phi} dphi) */
  double *theta      = alloc_double(Ntheta);
  double complex **v = calloc(Lmax+1,sizeof(*v));IsNull(v);/* v_{m}{theta} = \int f(theta,phi) exp^{-i m phi} dphi */
  Integration_T *I2  = init_integration();
  Uint i,j,l,m,lm;
  
  /* initialize tables */
  init_Legendre_root_function();
  
  /* making the integrands, note: array must be in order */
  for (i = 0; i < Ntheta; ++i)
  {
    Ftheta[i] = alloc_double(Nphi);
    theta[i]  = acos(-Legendre_root_function(i,Ntheta));/* note x starts from -1 and ends to 1 so in order so x = -cos(theta)*/
  }
  for (l = 0; l <= Lmax; ++l)
    for (m = 0; m <= l; ++m)
    {
      lm    = lm2n(l,m);
      for (i = 0; i < Ntheta; ++i)
        for (j = 0; j < Nphi; ++j)
          Ftheta[i][j] = f[j+i*Nphi];/* for each theta */
    }/* end of for (m = 0; m <= l; ++m) */
    
  /* carrying out the phi part of the integral, note: v(-m) = v*(m) */
  for (m = 0; m <= Lmax; ++m)
  { 
    v[m] = alloc_double_complex(Ntheta);
    int mp = (int)m;
    for (i = 0; i < Ntheta; ++i)
      /* for each theta and m calculate integral_{0}^{2pi} f(theta,phi)*exp(-i*m*phi) dphi */
      v[m][i] = integrate_expImphi(Ftheta[i],Nphi,-mp);
  }
  free_2d_mem(Ftheta,Ntheta);
  
  /* now carry out theta part of function: */
  /* use Gauss Legendre quadrature */
  I2->type           = "Gaussian Quadrature Legendre";
  I2->GQ_Legendre->n = Ntheta;
  plan_integration(I2);
  
  for (l = 0; l <= Lmax; ++l)
    for (m = 0; m <= l; ++m)/* m >= 0 */
    {
      lm = lm2n(l,m);
      for (i = 0; i < Ntheta; ++i)
      { 
        real_fYlm[i] = creal(Ylm((int)l,(int)m,theta[i],0))*creal(v[m][i]);
        imag_fYlm[i] = creal(Ylm((int)l,(int)m,theta[i],0))*cimag(v[m][i]);
      }
      I2->GQ_Legendre->f = real_fYlm;
      realClm[lm]        = execute_integration(I2);
      I2->GQ_Legendre->f = imag_fYlm;
      imagClm[lm]        = execute_integration(I2);
    }/* end of for (m = 0; m <= l; ++m) */
  
  /* free: */  
  free(real_fYlm);
  free(imag_fYlm);
  free(theta);
  free_2d_mem(v,Lmax+1);
  free_integration(I2);
}

/* ->return value: given point (theta,phi) and Ylm coeffs, 
// it interpolates to (theta,phi) */
double interpolation_Ylm(const double *const realClm,const double *const imagClm,const Uint Lmax, const double theta, const double phi)
{
  const double sign[2] = {1.,-1.};
  double complex sum = 0;
  Uint l,m,lm;
  
  for (l = 0; l <= Lmax; ++l)
  {
    for (m = 1; m <= l; ++m)
    {
      int mp = (int)m;
      lm   = lm2n(l,m);
      
      sum += (realClm[lm]+imagI*imagClm[lm])*Ylm((int)l,mp,theta,phi);/* m >= 1 */
      sum += sign[m%2]*(realClm[lm]-imagI*imagClm[lm])*Ylm((int)l,-mp,theta,phi);/* m < 0 */
    }
    lm   = lm2n(l,0);
    sum += (realClm[lm]+imagI*imagClm[lm])*Ylm((int)l,0,theta,phi);/* m == 0 */
  }

  return creal(sum);
}

/* ->return value: d(f(theta,phi))/dphi using Ylm expansion, 
// assuming Legendre root in theta direction and EquiSpaced in phi direction. */
double *df_dphi_Ylm(const double *const realClm,const double *const imagClm,const Uint Ntheta, const Uint Nphi,const Uint Lmax)
{
  double *df_dphi = alloc_double(Ntheta*Nphi);
  const double sign[2] = {1.,-1.};
  double theta,phi;
  double complex sum = 0;
  Uint i,j,l,m,lm;
  
  /* initialize tables */
  init_Legendre_root_function();
  
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      sum = 0;
      for (l = 0; l <= Lmax; ++l)
      {
        for (m = 1; m <= l; ++m)
        {
          int mp = (int)m;
          lm   = lm2n(l,m);
          sum += (realClm[lm]+imagI*imagClm[lm])*dYlm_dphi((int)l,mp,theta,phi);/* m >= 0 */
          sum += sign[m%2]*(realClm[lm]-imagI*imagClm[lm])*dYlm_dphi((int)l,-mp,theta,phi);/* m < 0 */
        }
        //lm   = lm2n(l,0);
        //sum += (realClm[lm]+imagI*imagClm[lm])*dYlm_dphi(l,0,theta,phi);/* m == 0 */
      }
      df_dphi[j+i*Nphi] = creal(sum);
    }/* end of for (j = 0; j < Nphi; ++j) */
  }
  
  return df_dphi;
}

/* ->return value: d(f(theta,phi))/dtheta using Ylm expansion, 
// assuming Legendre root in theta direction and EquiSpaced in phi direction. */
double *df_dtheta_Ylm(const double *const realClm,const double *const imagClm,const Uint Ntheta, const Uint Nphi,const Uint Lmax)
{
  double *df_dtheta    = alloc_double(Ntheta*Nphi);
  const double sign[2] = {1.,-1.};
  double theta,phi;
  double complex sum = 0;
  Uint i,j,l,m,lm;
  
  /* initialize tables */
  init_Legendre_root_function();
  
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      sum = 0;
      for (l = 0; l <= Lmax; ++l)
      {
        for (m = 1; m <= l; ++m)
        {
          int mp = (int)m;
          lm   = lm2n(l,m);
          sum += (realClm[lm]+imagI*imagClm[lm])*dYlm_dtheta((int)l,mp,theta,phi);/* m >= 0 */
          sum += sign[m%2]*(realClm[lm]-imagI*imagClm[lm])*dYlm_dtheta((int)l,-mp,theta,phi);/* m < 0 */
        }
        lm   = lm2n(l,0);
        sum += (realClm[lm]+imagI*imagClm[lm])*dYlm_dtheta((int)l,0,theta,phi);/* m == 0 */
      }
      df_dtheta[j+i*Nphi] = creal(sum);
    }/* end of for (j = 0; j < Nphi; ++j) */
  }
  
  return df_dtheta;
}

/* map: (l,m) -> n  for mapping 2-d array to 1-d array for -l <= m <= l */
int lm2n_Ylm(const int l,const int m, const int lmax)
{
  assert(abs(m) <= l);
  int i = m >= 0 ? lmax-l : lmax-l-m;
  int j = m <= 0 ? lmax-l : lmax-l+m;
  return j+(lmax+1)*i;
}

/* map: n = j+i*(lmax+1) -> (l,m) for mapping of 1-d array to 2-d array or -l <= m <= l */
void n2lm_Ylm(const int n, int *const l, int *const m,const int lmax)
{
  int j = n%(lmax+1);
  int i = (n-j)/(lmax+1);
  
  *m = j-i;
  *l = *m >= 0 ? lmax-i : lmax-j;
}

/* map: (l,m) -> n  for mapping 2-d array to 1-d array for 0 <= m <= l, 
// the inverse of this map won't be needed.
// note: the order is C_{0}^{0},C_{1}^{0},C_{1}^{1},C_{2}^{0},C_{2}^{1},C_{2}^{2}, ...  */
Uint lm2n(const Uint l, const Uint m)
{
  return l*(l+1)/2+m;
}

/* ->return value \integral_{0}^{2pi} f(phi)*exp(imagI*m*phi)dphi (trapezoidal rule) */
static double complex integrate_expImphi(const double *const f, const Uint n/* f array dimension */,const int m)
{
  double complex i0 = 0;
  const double complex phi0 = 2.*imagI*M_PI/n;
  Uint i;
  
  for (i = 0; i < n; ++i)
    i0 += f[i]*cexp(m*phi0*i);
  i0 *= 2.*M_PI/n;
 
  return i0;
}

/* ->return value: allocating memory for Clm coeffs of Ylm expansion 
// for given Lmax, which is (Lmax+1)*Lmax/2 + Lmax+1. */
double *alloc_ClmYlm(Uint Lmax)
{
  return alloc_double((Lmax+1)*Lmax/2 + Lmax+1);
}
