/*
// Alireza Rashti
// July 2018
*/

#include "fourier_transformation.h"


/* 1-D regular Fourier transformation -as opposed to 
// fast Fourier transformation.
// getting values and number of values,
// it finds NORMALIZED coefficients of a function with Chebyshev basis of 
// first kind expansion with Chebyshev Extrema collocation.
//
// what it calculates are c's which stand for coeffs as follows :
// v(i) = c(0) + (-1)^i*c(n-1) + 2 \sum_{j=1}^{n-2} c(j)*cos(i*j*PI/(n-1))
// the algorithm is simple : one needs to write the above c's in 
// terms of v's as the following:
// c(i)*scale = v(0) + (-1)^i*v(n-1) + 2 \sum_{j=1}^{n-2} v(j)*cos(i*j*PI/(n-1))
// in wich scale is 2*(n-1).
// note: the coeffs are NORMALIZED.
*/
void rft_1d_ChebyshevExtrema_coeffs(double *const values ,double *const coeffs,const Uint n)
{
  Uint i,j;
  const double SIGN[2] = {1.0,-1.0};
  const double th0 = M_PI/(n-1);
  const double scale = 2*(n-1);
  
  for (i = 0; i < n; ++i)
  {
    double th = i*th0;
    double c = 0;/* coeffs */
    for (j = 1; j < n-1; ++j)
      c += values[j]*cos(j*th);
    
    c *= 2;
    c += values[0]+SIGN[i%2]*values[n-1];
    coeffs[i] = c/scale;
  }
}


/* 1-D regular Fourier transformation -as opposed to 
// fast Fourier transformation.
// getting values and number of values,
// it finds NORMALIZED coefficients of a function with Chebyshev basis of 
// first kind expansion with Chebyshev Nodes collocation
//
// what it calculates are c's which stand for coeffs as follows :
// v(i) = c(0)+2*\sum_{j=1}^{n-1} c(j)*cos(j*PI(i+1/2)/n)
// the algorithm is simple : one needs to write the above c's in 
// terms of v's as the following:
// c(i)*scale = 2*\sum_{j=0}^{n-1} v(j)*cos(i*PI(j+1/2)/n)
// in wich scale is 2n.
// note: the coeffs are NORMALIZED.
*/
void rft_1d_ChebyshevNodes_coeffs(double *const values ,double *const coeffs,const Uint n)
{
  Uint i,j;
  const double th0 = M_PI/n;
  const double scale = n;
  
  for (i = 0; i < n; ++i)
  {
    double th = i*th0;
    double c = 0;/* coeffs */
    for (j = 0; j < n; ++j)
      c += values[j]*cos(th*(j+0.5));
      
    coeffs[i] = c/scale;
  }
}

/* fourier transformation from real value to complex coeffs:
// f(x) = \sum_{m=-l+1}^{l-1} c(m)*exp(imagI*m*x), where x is in [0,2*pi],
// and l = n/2+1 (if n is odd, it is rounded down) thus: we have:
// c(m) = 1/(2*pi)*\integral_{0}^{2*pi} f(x)*exp(-imagI*m*x) dx.
//
// some notes:
// ============
// note: l = n/2 +1 if n is odd it is rounded down.
// note: since f(x) is real we have: c(-m) = c*(m) 
//       and we only calculate c(j), j = 0,...,l-1
// note: argument n is the number of grid points in x direction.
// note: of course, works best when f function is periodic in [0,2*pi]
// note: it allocates memory so at some point you need to free it.
// note: we use trapezoidal rule to carry out the intergral. 
// ->return value: c(m) */
void *r2cft_1d_EquiSpaced_coeffs(const double *const value,const Uint n)
{
  if (n == 0)
    Error0("Fourier Transformation: No points!\n");
  
  const Uint l = n/2+1;/* note: if n is not even, it is rounded down */
  double complex *const coeffs = alloc_double_complex(l);
  const double complex x0 = -2.*imagI*M_PI/n;/* - included */
  Uint m;
  
  for (m = 0; m < l; ++m)/* note: l is excluded, otherwise one have aliasing and then error */
  {
    Uint i;
    
    coeffs[m] = 0;
    for (i = 0; i < n; ++i)
      coeffs[m] += value[i]*cexp(m*i*x0);
    
    coeffs[m] /= n;
  }
  
  return coeffs;
}

/* fourier transformation from complex coeffs to real values:
// f(x) = \sum_{m=-l+1}^{l-1} c(m)*exp(imagI*m*x), where x is in [0,2*pi],
// and l = n/2+1 (if n is odd, it is rounded down) thus: we have:
// c(m) = 1/(2*pi)*\integral_{0}^{2*pi} f(x)*exp(-imagI*m*x) dx. 
// ->return value : f(x) */
double *c2rft_1d_EquiSpaced_values(void *const coeffs,const Uint N)
{
  if (N == 0)
    Error0("Fourier Transformation: No points!\n");
  
  const double complex *const c = coeffs;
  double *f = alloc_double(N);
  const Uint l = N/2+1;
  const double complex x0 = 2.*imagI*M_PI/N;
  Uint i,j;
  
  for (i = 0; i < N; ++i)
  {
    double complex x = i*x0;
    f[i] = creal(c[0]);
    for (j = 1; j < l; ++j)
      f[i] += 2.*creal(c[j]*cexp((double)j*x));
  }

  return f;
}




/* fourier transformation from real value to complex coeffs for 2d on S2
// this is double fourier transformation on S2, since the given function f
// is extended on the whole sphere:
//                 | f(theta,phi),                         theta in [0,pi] and phi in [0,2pi)
// f`(theta,phi) = |
//                 | f(theta+pi,phi) = f(pi-theta,phi+pi), theta+pi in [pi,2pi), phi in [0,2pi)
//                 | or 
//                 | f(theta+pi,phi) = f(pi-theta,phi), theta+pi in [pi,2pi), phi in [0,2pi) if this is better.
//
//
// notes:
// o. Nphi MUST be even
// o. values of f at theta = pi MUST be given.
// o. it is cursed to slowly converge!
// o. f expansion => f(theta,phi) = F(phi0,phi1) =
//    \sum_{m0=-l0,m1=-l1}^{m0=l0,m1=l1}{Cm0m1 exp(I.m0.phi0) exp(I.m1.phi1)}.
//    =>  Cm0m1 = 1/(2*pi)^2 *\integral_{0}^{2*pi}\integral_{0}^{2*pi} 
//              f(phi0,phi1) exp(-imagI*m0*phi0) exp(-imagI*m1*phi1) dphi0 dphi1.
// o. theta is in [0,pi] and phi is in [0,2pi]
// o. theta = phi0/2 and phi = phi1
// o. phi1 and phi2 are in [0,2 pi]
// o. collocation poinst are EquiSpaced for theta and phi
// o. f(theta(i),phi(j)) = f[i][j] = f[IJ(i,j,Nphi)], where IJ is the macro in the header
// o. Cm0m1's are composed of two parts, Cr[IJ(m0,m1,l1)] and Ci[l0*l1+IJ(m0,m1,l1)]
// o. syntax:
// =========
// double f = data;
// double *realC,*imagC;
//
// r2cft_2d_coeffs_S2(f,Ntheta,Nphi,&realC,&imagC);
// ...
// free(realC);
// free(imagC);
//
// ->: Cm0m1  if improve==1 it returns the max magnitude of 
// the last coeffs otherwise DBL_MAX. */
double
r2cft_2d_coeffs_S2
(
  const double *const f/* field values given on theta and phi coords. */,
  Uint Ntheta/* number of point in theta direction */, 
  const Uint Nphi/* number of point in phi direction */,
  double **const realC/* real part of coeffs, allocates memory */,
  double **const imagC/* imag part of coeffs, allocates memory*/,
  const int improve/* if 1, it tries to improve the expansion, otherwise no. */
)
{
  if (!f)
    Error0("Bad argument: no value\n!");
  
  if (Nphi%2)
    Error0("Number of points in phi direction must be even.\n!");
  
  Ntheta -= 1;/* adjust Ntheta */
  const double COEFF_THRESHOLD = 1E-6;/* if coeffs error is bigger than this 
                                      // change the continuation method */
  const Uint TwiceNtheta = 2*Ntheta;
  double *const F = alloc_double(TwiceNtheta*Nphi); IsNull(F);
  double ret = DBL_MAX;
  Uint ij,i,j,k,l;
  
  /* f(theta,phi), theta in [0,pi] and phi in [0,2pi) */
  for (i = 0; i < Ntheta; ++i)
  {
    for (j = 0; j < Nphi; ++j)
    {
      ij        = IJ(i,j,Nphi);
      F[ij]     = f[ij];
    }
  }
  
  /* f(theta+pi,phi) = f(pi-theta,phi+pi), theta+pi in [pi,2pi), phi in [0,2pi) */
  for (i = Ntheta; i < TwiceNtheta; ++i)
  {
    k = TwiceNtheta-i;
    for (j = 0; j < Nphi; ++j)
    {
      l               = (j+Nphi/2)%Nphi;
      F[IJ(i,j,Nphi)] = f[IJ(k,l,Nphi)];
    }
  }
  r2cft_2d_coeffs(F,TwiceNtheta,Nphi,realC,imagC);
  
  /* check the coeffs and if possible improve the expansion
  // by using different continuation */
  if (improve)
  {
    double *realC2 = 0,*imagC2 = 0;
    double max1,max2;
    
    max1 = r2cft_2d_last_coeffs_max_mag_S2(TwiceNtheta,Nphi,*realC,*imagC);
    ret  = max1;
    /* if passes the threshold */
    if (max1 > COEFF_THRESHOLD)
    {
      /* f(theta+pi,phi) = f(pi-theta,phi), theta+pi in [pi,2pi), phi in [0,2pi) */
      for (i = Ntheta; i < TwiceNtheta; ++i)
      {
        k = TwiceNtheta-i;
        for (j = 0; j < Nphi; ++j)
        {
          F[IJ(i,j,Nphi)] = f[IJ(k,j,Nphi)];
        }
      }
      
      r2cft_2d_coeffs(F,TwiceNtheta,Nphi,&realC2,&imagC2);
      max2 = r2cft_2d_last_coeffs_max_mag_S2(TwiceNtheta,Nphi,realC2,imagC2);
      
      if (max2 > max1)
      {
        free(realC2);
        free(imagC2);
        ret = max1;
      }
      else
      {
        free(*realC);
        free(*imagC);
        *realC = realC2;
        *imagC = imagC2;
        ret = max2;
      }
    }/* if (max1 > COEFF_THRESHOLD) */
  }/* if (improve) */
  
  free(F);
  return ret;
}

/* -> the max magnitude of the last few coeffs.
// find the max magnitude of the last few coeffs to estimate the error 
// in the expansion for r2cft_2d_coeffs_S2 */
static double 
r2cft_2d_last_coeffs_max_mag_S2
(
  Uint Ntheta/* number of point in theta direction */, 
  const Uint Nphi/* number of point in phi direction */,
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs,*/
)
{
  const Uint FEW  = 2;/* the last few 'FEW' coeffs */
  const Uint l0   = Ntheta/2+1;/* note: there are TwiceNtheta coeffs */
  const Uint l1   = Nphi/2+1;
  const Uint l0l1 = l0*l1;
  double max_r,max_i,abs_coeff;
  Uint m0,m1,m0m1;
  
  assert(l0>=FEW);
  assert(l1>=FEW);
  
  /* check the last few coeffs and make sure the last coeffs are small */
  max_r = max_i = 0;
  for (m0 = l0-FEW; m0 < l0; ++m0)
  {
    for (m1 = 0; m1 < l1; ++m1)
    {
      m0m1 = IJ(m0,m1,l1);
      
      abs_coeff = fabs(realC[m0m1]);
      if (abs_coeff > max_r)
        max_r = abs_coeff;
        
      abs_coeff = fabs(imagC[m0m1]);
      if (abs_coeff > max_i)
        max_i = abs_coeff;
    }
  }
  for (m0 = 0; m0 < l0; ++m0)
  {
    for (m1 = l1-FEW; m1 < l1; ++m1)
    {
      m0m1 = IJ(m0,m1,l1);
      
      abs_coeff = fabs(realC[m0m1]);
      if (abs_coeff > max_r)
        max_r = abs_coeff;
        
      abs_coeff = fabs(imagC[m0m1]);
      if (abs_coeff > max_i)
        max_i = abs_coeff;
    }
  }
  for (m0 = l0-FEW; m0 < l0; ++m0)
  {
    for (m1 = 0; m1 < l1; ++m1)
    {
      m0m1 = IJ(m0,m1,l1)+l0l1;
      
      abs_coeff = fabs(realC[m0m1]);
      if (abs_coeff > max_r)
        max_r = abs_coeff;
        
      abs_coeff = fabs(imagC[m0m1]);
      if (abs_coeff > max_i)
        max_i = abs_coeff;
    }
  }
  for (m0 = 0; m0 < l0; ++m0)
  {
    for (m1 = l1-FEW; m1 < l1; ++m1)
    {
      m0m1 = IJ(m0,m1,l1)+l0l1;
      
      abs_coeff = fabs(realC[m0m1]);
      if (abs_coeff > max_r)
        max_r = abs_coeff;
        
      abs_coeff = fabs(imagC[m0m1]);
      if (abs_coeff > max_i)
        max_i = abs_coeff;
    }
  }
  return MaxMag_d(max_r,max_i);
}

/* fourier transformation from real value to complex coeffs for 2d.
// notes:
// o. f expansion => f(phi0,phi1) = 
//    \sum_{m0=-l0,m1=-l1}^{m0=l0,m1=l1}{Cm0m1 exp(I.m0.phi0) exp(I.m1.phi1)}.
//    =>  Cm0m1 = 1/(2*pi)^2 *\integral_{0}^{2*pi}\integral_{0}^{2*pi} 
//              f(phi0,phi1) exp(-imagI*m0*phi0) exp(-imagI*m1*phi1) dphi0 dphi1.
// o. phi1 and phi2 are in [0,2 pi]
// o. collocation poinst are EquiSpaced
// o. f(phi0(i),phi1(j)) = f[i][j] = f[IJ(i,j,Nphi1)], where IJ is the macro in the header
// o. Cm0m1's are composed of two parts, Cr[IJ(m0,m1,l1)] and Ci[l0*l1+IJ(m0,m1,l1)]
// o. l0 = Nphi0/2+1 since f is real, note: if Nphi0 is odd it is rounded down
// o. l1 = Nphi1/2+1 since f is real, note: if Nphi1 is odd it is rounded down
// o. syntax:
// =========
// double f = data;
// double *realC,*imagC;
//
// r2cft_2d_coeffs(f,Nphi0,Nphi1,&realC,&imagC);
// ...
// free(realC);
// free(imagC);
//
// ->: Cm0m1  */
void
r2cft_2d_coeffs
(
  const double *const f/* field values */,
  const Uint Nphi0/* number of point in phi0 direction */, 
  const Uint Nphi1/* number of point in phi1 direction */,
  double **const realC/* real part of coeffs, allocates memory */,
  double **const imagC/* imag part of coeffs, allocates memory*/
)
{
  if (!f)
    Error0("Bad argument: no value\n!");
    
  const Uint l0   = Nphi0/2+1;
  const Uint l1   = Nphi1/2+1;
  const Uint l0l1 = l0*l1;
  const double complex x0 = -2.*imagI*M_PI/Nphi0;/* - included */
  const double complex x1 = -2.*imagI*M_PI/Nphi1;/* - included */
  double *const Rc = alloc_double(2*l0l1);
  double *const Ic = alloc_double(2*l0l1);
  double *crr      = alloc_double(l0l1);
  double *cri      = alloc_double(l0l1);
  double *cir      = alloc_double(l0l1);
  double *cii      = alloc_double(l0l1);
  double **cfr     = calloc(Nphi0,sizeof(*cfr));IsNull(cfr);
  double **cfi     = calloc(Nphi0,sizeof(*cfi));IsNull(cfi);
  Uint i,j,m0,m1;
  
  /* FT in 2nd index */
  for (i = 0; i < Nphi0; ++i)
  {
    cfr[i] = alloc_double_complex(l1);
    cfi[i] = alloc_double_complex(l1);
    for (m1 = 0; m1 < l1; ++m1)
    {
      double complex m1x1 = m1*x1;
      double complex cf   = 0;
      for (j = 0; j < Nphi1; ++j)
          cf += f[IJ(i,j,Nphi1)]*cexp(j*m1x1);
      cf /= Nphi1;
      cfr[i][m1] = creal(cf);
      cfi[i][m1] = cimag(cf);
    }
    cfr[i][0] /= 2;
    cfi[i][0] /= 2;
  }
  
  /* FT for each component in 1st index */
  for (m1 = 0; m1 < l1; ++m1)
  {
    for (m0 = 0; m0 < l0; ++m0)
    {
      double complex m0x0 = m0*x0;
      double complex cr = 0, ci = 0;
      Uint m0m1 = IJ(m0,m1,l1);
      for (i = 0; i < Nphi0; ++i)
      {
        cr += cfr[i][m1]*cexp(i*m0x0);
        ci += cfi[i][m1]*cexp(i*m0x0);
      }
      cr /= Nphi0;
      ci /= Nphi0;
      crr[m0m1] = creal(cr);
      cri[m0m1] = cimag(cr);
      cir[m0m1] = creal(ci);
      cii[m0m1] = cimag(ci);
    }
  }
  m0 = 0;
  for (m1 = 0; m1 < l1; ++m1)
  {
    Uint m0m1 = IJ(m0,m1,l1);
    crr[m0m1] /= 2;
    cri[m0m1] /= 2;
    cir[m0m1] /= 2;
    cii[m0m1] /= 2;
  }
  /* decompose real and imag parts */
  for (m1 = 0; m1 < l1; ++m1)
    for (m0 = 0; m0 < l0; ++m0)
    {
      Uint m0m1 = IJ(m0,m1,l1);
      Rc[m0m1]= crr[m0m1];
      Ic[m0m1]= cri[m0m1];
    }
  
  /* decompose real and imag parts */
  for (m1 = 0; m1 < l1; ++m1)
    for (m0 = 0; m0 < l0; ++m0)
    {
      Uint m0m1 = IJ(m0,m1,l1);
      Rc[m0m1+l0l1]= cir[m0m1];
      Ic[m0m1+l0l1]= cii[m0m1];
    }
  
  free_2d_mem(cfr,Nphi0);
  free_2d_mem(cfi,Nphi0);
  free(crr);
  free(cri);
  free(cir);
  free(cii);
  *realC = Rc;
  *imagC = Ic;
}

/* -> interpolation at (theta,phi) using 2-d Fourier transformation on S2 */
double 
r2cft_2d_interpolation_S2
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Ntheta/* number of point in theta direction */,
  const Uint Nphi/* number of point in phi direction */,
  const double theta/* point of interest at theta dir */,
  const double phi/* point of interest at phi dir */
)
{
  return r2cft_2d_interpolation(realC,imagC,2*(Ntheta-1),Nphi,theta,phi);
}

/* -> interpolation at (phi0,phi1) using 2-d Fourier transformation 
// r2cft_2d.
// note: phi0 and phi1 must be in radian. */
double 
r2cft_2d_interpolation
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Nphi0/* number of point in phi0 direction */,
  const Uint Nphi1/* number of point in phi1 direction */,
  const double phi0/* point of interest at phi0 dir */,
  const double phi1/* point of interest at phi0 dir */
)
{
  if(!realC || !imagC)
    Error0("Bad argument: no coefficients!\n");
    
  const Uint l0 = Nphi0/2+1;
  const Uint l1 = Nphi1/2+1;
  const Uint l0l1 = l0*l1;
  double complex interp = 0;
  Uint m0,m1;
  
  /* sum */
  for (m0 = 0; m0 < l0; ++m0)
  {
    for (m1 = 0; m1 < l1; ++m1)
    {
      Uint m0m1 = IJ(m0,m1,l1);
      interp += (realC[m0m1]    - imagC[l0l1+m0m1]+
                 imagI*(imagC[m0m1] + realC[l0l1+m0m1]))*
                 cexp(imagI*((double)m0*phi0+(double)m1*phi1));
      interp += (realC[m0m1]    + imagC[l0l1+m0m1] +
                 imagI*(imagC[m0m1] - realC[l0l1+m0m1]))*
                 cexp(imagI*((double)m0*phi0-(double)m1*phi1));
    }
  }
  return 2*creal(interp);
}

/* -> taking derivative : df(theta,phi)/dtheta on S2. */
double *
r2cft_2d_df_dtheta_S2
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Ntheta/* number of point in theta direction */,
  const Uint Nphi/* number of point in phi direction */
)
{
  return r2cft_2d_df_dphi0(realC,imagC,2*(Ntheta-1),Nphi);
}

/* -> taking derivative : df(theta,phi)/dphi on S2. */
double *
r2cft_2d_df_dphi_S2
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Ntheta/* number of point in theta direction */,
  const Uint Nphi/* number of point in phi direction */
)
{
  return r2cft_2d_df_dphi1(realC,imagC,2*(Ntheta-1),Nphi);
}

/* -> taking derivative : df(phi0,phi1)/dphi0. */
double *
r2cft_2d_df_dphi0
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Nphi0/* number of point in phi0 direction */,
  const Uint Nphi1/* number of point in phi1 direction */
)
{
  if(!realC || !imagC)
    Error0("Bad argument: no coefficients!\n");
    
  const Uint l0 = Nphi0/2+1;
  const Uint l1 = Nphi1/2+1;
  const Uint l0l1 = l0*l1;
  const double x0 = 2.*M_PI/Nphi0;
  const double x1 = 2.*M_PI/Nphi1;
  double *df        = alloc_double(Nphi0*Nphi1);
  Uint i,j,m0,m1,ij,m0m1;
  
  for (i = 0; i < Nphi0; ++i)
  {
    double phi0 = i*x0;
    for (j = 0; j < Nphi1; ++j)
    {
      double phi1        = j*x1;
      double complex dfc = 0;
      ij = IJ(i,j,Nphi1);
      for (m0 = 0; m0 < l0; ++m0)
      {
        double complex Im0 = imagI*(double)m0;
        double complex expIm0phi0 = cexp(Im0*phi0);
        for (m1 = 0; m1 < l1; ++m1)
        {
          m0m1 = IJ(m0,m1,l1);
          dfc += Im0*(realC[m0m1]   - imagC[l0l1+m0m1]+
                     imagI*(imagC[m0m1] + realC[l0l1+m0m1]))*
                     expIm0phi0*cexp(imagI*(double)m1*phi1);
          dfc += Im0*(realC[m0m1]   + imagC[l0l1+m0m1] +
                     imagI*(imagC[m0m1] - realC[l0l1+m0m1]))*
                     expIm0phi0*cexp(-imagI*(double)m1*phi1);
        }
      }
      df[ij] = 2*creal(dfc);
    }
  }
  
  return df;
}

/* -> taking derivative : df(phi0,phi1)/dphi1. */
double *
r2cft_2d_df_dphi1
(
  const double *const realC/* real part of coeffs */,
  const double *const imagC/* imag part of coeffs */,
  const Uint Nphi0/* number of point in phi0 direction */,
  const Uint Nphi1/* number of point in phi1 direction */
)
{
  if(!realC || !imagC)
    Error0("Bad argument: no coefficients!\n");
    
  const Uint l0 = Nphi0/2+1;
  const Uint l1 = Nphi1/2+1;
  const Uint l0l1 = l0*l1;
  const double x0 = 2.*M_PI/Nphi0;
  const double x1 = 2.*M_PI/Nphi1;
  double *df        = alloc_double(Nphi0*Nphi1);
  Uint i,j,m0,m1,ij,m0m1;
  
  for (i = 0; i < Nphi0; ++i)
  {
    double phi0 = i*x0;
    for (j = 0; j < Nphi1; ++j)
    {
      double phi1        = j*x1;
      double complex dfc = 0;
      ij = IJ(i,j,Nphi1);
      for (m1 = 0; m1 < l1; ++m1)
      {
        double complex Im1 = imagI*(double)m1;
        double complex expIm1phi1 = cexp(Im1*phi1);
        for (m0 = 0; m0 < l0; ++m0)
        {
          m0m1 = IJ(m0,m1,l1);
          dfc += Im1*(realC[m0m1]   - imagC[l0l1+m0m1]+
                     imagI*(imagC[m0m1] + realC[l0l1+m0m1]))*
                     cexp(imagI*(double)m0*phi0)*expIm1phi1;
          dfc += -Im1*(realC[m0m1]  + imagC[l0l1+m0m1] +
                     imagI*(imagC[m0m1] - realC[l0l1+m0m1]))*
                     cexp(imagI*(double)m0*phi0)/expIm1phi1;
        }
      }
      df[ij] = 2*creal(dfc);
    }
  }
  
  return df;
}

