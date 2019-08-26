/*
// Alireza Rashti
//  August 2019
*/

#include "spherical_harmonics_transformation.h"

/* calculatinf coefficents of spherical harmonic Clm for:
// f(theta,phi) = \sum_{l=0}^{l=lmax} \sum_{m = -l}^{l} C_{lm}*Y_{lm}(thata,phi)
// Y_{lm}=\left( -1\right) ^{m}\sqrt {\frac {2l+1} {4\pi }\frac {\left( l-m\right) !} {\left( l+m\right) !}}P_{1}^{m}\left( \cos \theta \right) e^{im\varphi }
// 
// some notes:
// ===========
// o. l is inclusive
// o. f(theta, phi) is given in a 1d array with this mapping:
//    f(i,j) = f[j+i*Nphi]
// o. Clm = realC[lm2n(l,m,lmax)] + imagC[lm2n(l,m,lmax)]*I;
// o. Clm is only calculated for m >= 0, since C(l,-m) = (-1)^m C*(l,m)
// o. since m >= 0 the total memory is N = (Lmax+1)*Lmax/2 + Lmax+1
*/
void 
get_Ylm_coeffs(
double *const realClm/* Re(C_{lm}) coeffs */,
double *const imagClm/* Im(C_{lm}) coeffs */,
const double *const f/* f(theta,phi) */,
const unsigned Ntheta/* Ntheta points in theta direction */,
const unsigned Nphi/* Nphi points in phi direction*/,
const unsigned Lmax/* maximum l (inclusive) for the expansion */
)
{
  const char *collocation = "GaussLegendre_EquiSpaced";
  /* Note: later on you can add parameter to direct the collocation */
  if (strcmp_i(collocation,"GaussLegendre_EquiSpaced"))
    get_Ylm_coeffs_GaussLegendre_EquiSpaced(realClm,imagClm,f,Ntheta,Nphi,Lmax);
  else
    abortEr(NO_JOB);

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
const unsigned Ntheta/* Ntheta points in theta direction */,
const unsigned Nphi/* Nphi points in phi direction*/,
const unsigned Lmax/* maximum l (inclusive) for the expansion */)
{
  const unsigned N = (Lmax+1)*Lmax/2 + Lmax+1;
  double **Ftheta  = calloc(Ntheta,sizeof(*Ftheta));pointerEr(Ftheta);
  double *g_r      = alloc_double(Ntheta);
  double *g_i      = alloc_double(Ntheta);
  double **g       = calloc(N,sizeof(*g));pointerEr(g);
  double *theta    = alloc_double(Ntheta);
  double complex **v = calloc(Lmax+1,sizeof(*v));pointerEr(v);
  Integration_T *I2  = init_integration();
  unsigned i,j,l,m,lm;
  
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
      g[lm] = alloc_double(Ntheta);
      for (i = 0; i < Ntheta; ++i)
      {
        g[lm][i]  = creal(Ylm(l,(int)m,theta[i],0));/* for phi = 0 Ylm is real */
        for (j = 0; j < Nphi; ++j)
          Ftheta[i][j] = f[j+i*Nphi];/* for each theta */
      }
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
        g_r[i] = g[lm][i]*creal(v[m][i]);
        g_i[i] = g[lm][i]*cimag(v[m][i]);
      }
      I2->GQ_Legendre->f = g_r;
      realClm[lm]        = execute_integration(I2);
      I2->GQ_Legendre->f = g_i;
      imagClm[lm]        = execute_integration(I2);
    }/* end of for (m = 0; m <= l; ++m) */
  
  /* free: */  
  free(g_r);
  free(g_i);
  free(theta);
  free_2d_mem(g,N);
  free_2d_mem(v,Lmax+1);
  free_integration(I2);
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
unsigned lm2n(const unsigned l, const unsigned m)
{
  return l*(l+1)/2+m;
}

/* ->return value \integral_{0}^{2pi} f(phi)*exp(I*m*phi)dphi (trapezoidal rule) */
static double complex integrate_expImphi(const double *const f, const unsigned n/* f array dimension */,const int m)
{
  double complex i0 = 0;
  const double complex phi0 = 2.*I*M_PI/n;
  unsigned i;
  
  for (i = 0; i < n; ++i)
    i0 += f[i]*cexp(m*phi0*i);
  i0 *= 2.*M_PI/n;
 
  return i0;
}

/* ->return value: allocating memory for Clm coeffs of Ylm expansion 
// for given Lmax, which is (Lmax+1)*Lmax/2 + Lmax+1. */
double *alloc_ClmYlm(unsigned Lmax)
{
  return alloc_double((Lmax+1)*Lmax/2 + Lmax+1);
}
