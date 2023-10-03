/*
// Alireza Rashti
// July 2018
*/

#include "spectral_tests.h"
//#include "interpolations.h"
//#include "interpolations.c"
#include <time.h>
#define ArgM(a) a,#a/*used for being more accurate in naming and fast */
#define MAXSTR (400)

/* testing Fourier transformation functions.
// ->return value: EXIT_SUCCESS */
int fourier_transformation_tests(Grid_T *const grid)
{
  
  if (DO)
  {
    printf("Fourier transformation test: real to complex Fourier transformation 1-D: \n");
    cft_c2r_r2c_1d_EquiSpaced_test(grid);
  }
  if (DO)
  {
    printf("Fourier transformation test: real to complex Fourier transformation 2-D: \n");
    r2cft_2d_EquiSpaced_test(grid);
  }
  if (DO)
  {
    printf("Fourier transformation test: real to complex Fourier transformation 2-D on S2: \n");
    r2cft_2d_EquiSpaced_S2_test(grid);
  }
  
  return EXIT_SUCCESS;
}

/* testing spherical harmonica transformation functions.
// ->return value: EXIT_SUCCESS */
int Ylm_transformation_tests(Grid_T *const grid)
{
  
  if (DO)
  {
    printf("\nSpherical harmonic transformation test: \n");
    Ylm_trans_test(grid);
  }
  if (DO)
  {
    printf("\nSpherical harmonic derivatives test: \n");
    Ylm_derivatives_test(grid);
  }
  
  return EXIT_SUCCESS;
}

/* testing : get_Ylm_coeffs function.
// ->return value: TEST_SUCCESSFUL */
static int Ylm_trans_test(Grid_T *const grid)
{
  const char *const par = Pgets("Test_Ylm_Transformation");
  const double sign[2] = {1.,-1.};
  Uint Ntheta = 0;
  Uint Nphi   = 0;
  Uint lmax   = 0;
  Uint N;
  double *f;/* f(theta,phi) theta = [0,phi], phi = [0,2pi] */
  double  *realClm;/* f(theta,phi) = Re(Clm)*Ylm(theta,phi) */
  double  *imagClm;/* f(theta,phi) = Im(Clm)*Ylm(theta,phi) */
  double phi,theta,df;
  Uint Smem;
  Uint i,j,l,m,lm;
  
  /* initialize tables */
  init_Legendre_root_function();
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    N = (Uint)atoi(s);
    Free(s);
  }
  else
    N = 10;
  
  /* setting grid size */  
  Ntheta = Nphi = N;
  lmax = (Ntheta-1)/2;
  Smem = (lmax+1)*lmax/2 + lmax+1;
  
  /* allocating memory */
  f = alloc_double(Ntheta*Nphi);
  realClm = alloc_double(Smem);
  imagClm = alloc_double(Smem);
  
  /* populating f(theta,phi) */  
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      //f[j+Nphi*i] = (sin(theta)*cos(phi)-sin(theta)*sin(2*phi)); (this required high Lmax)
      f[j+Nphi*i] = creal(Ylm(4,3,theta,phi));
    }
  }
  /* calculating coeffs */
  get_Ylm_coeffs(realClm,imagClm,f,Ntheta,Nphi,lmax);
  
  df = 0.;
  /* now let's see how Ylm sum works: */
  printf("--> For Ntheta = %u, Nphi = %u, Lmax = %u:\n",Ntheta,Nphi,lmax);
  printf("sum:\n");
  /* for each theta and phi */
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      double fr = f[j+Nphi*i];
      
      /* sum: */
      double complex sum = 0;
      for (l = 0; l <= lmax; ++l)
      {
        for (m = 1; m <= l; ++m)
        {
          int mp = (int)m;
          lm   = lm2n(l,m);
          
          sum += (realClm[lm]+imagI*imagClm[lm])*Ylm((int)l,mp,theta,phi);/* m >= 0 */
          sum += sign[m%2]*(realClm[lm]-imagI*imagClm[lm])*Ylm((int)l,-mp,theta,phi);/* m < 0 */
        }
        lm   = lm2n(l,0);
        sum += (realClm[lm]+imagI*imagClm[lm])*Ylm((int)l,0,theta,phi);/* m == 0 */
      }
      if (df < fabs(creal(sum)-fr))
        df = fabs(creal(sum)-fr);
    }
  }/* for (i = 0; i < Ntheta; ++i) */
  printf("Max Error = %e\n",df);
  
  /* let's do some interpolation: */
  if (1)
  {
    double *ran_theta = make_random_number(Ntheta,0,M_PI);
    double *ran_phi   = make_random_number(Nphi,0,2*M_PI);
    printf("Interpolation:       f(x):                diff:\n");
    df = 0;
    for (i = 0; i < Ntheta; ++i)
    {
      double x = ran_theta[i];
      for (j = 0; j < Nphi; ++j)
      {
        double y = ran_phi[j];
        double fr = creal(Ylm(4,3,x,y));
        double fi = interpolation_Ylm(realClm,imagClm,lmax,x,y);
        printf("%+0.15f   %+0.15f   %+e\n",fi,fr,fi-fr);
        if (df<fabs(fi-fr))
          df = fabs(fi-fr);
      }
    }
    printf("Max Error = %e\n",df);
    Free(ran_theta);
    Free(ran_phi);
  }
  
  free(f);
  free(realClm);
  free(imagClm);
  UNUSED(grid);
  
  return TEST_SUCCESSFUL;
}

/* testing :  df_dphi_Ylm and df_dtheta_Ylm functions.
// ->return value: TEST_SUCCESSFUL */
static int Ylm_derivatives_test(Grid_T *const grid)
{
  const char *const par = Pgets("Test_Ylm_Transformation");
  Uint Ntheta = 0;
  Uint Nphi   = 0;
  Uint lmax   = 0;
  Uint N;
  double *f;/* analytic: f(theta,phi), theta = [0,phi], phi = [0,2pi] */
  double *f_dphi;/* analytic: df(theta,phi)/dphi, theta = [0,phi], phi = [0,2pi] */
  double *f_dtheta;/* analytic: df(theta,phi)/dtheta, theta = [0,phi], phi = [0,2pi] */
  double *df_dphi;/* numeric: df(theta,phi)/dphi */
  double *df_dtheta;/* numeric: df(theta,phi)/dtheta */
  double *realClm;/* f(theta,phi) = Re(Clm)*Ylm(theta,phi) */
  double *imagClm;/* f(theta,phi) = Im(Clm)*Ylm(theta,phi) */
  double phi,theta,df;
  Uint Smem;
  Uint i,j;
  
  /* initialize tables */
  init_Legendre_root_function();
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    N = (Uint)atoi(s);
    Free(s);
  }
  else
    N = 14;
  
  /* setting grid size */  
  Ntheta = Nphi = N;
  lmax = (Ntheta-1)/2;
  Smem = (lmax+1)*lmax/2 + lmax+1;
  
  /* allocating memory */
  f        = alloc_double(Ntheta*Nphi);
  f_dphi   = alloc_double(Ntheta*Nphi);
  f_dtheta = alloc_double(Ntheta*Nphi);
  realClm  = alloc_double(Smem);
  imagClm  = alloc_double(Smem);
  
  /* populating f(theta,phi) */  
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      f[j+Nphi*i]        = creal(Ylm(6,-5,theta,phi));
      f_dphi[j+Nphi*i]   = creal(-15.0/32.0*sqrt(1001)*imagI*cexp(-5*imagI*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI));
      f_dtheta[j+Nphi*i] = creal((3.0/32.0)*sqrt(1001)*(5 - 6*pow(sin(theta), 2))*cexp(-5*imagI*phi)*pow(sin(theta), 4)/sqrt(M_PI));
    }
  }
  /* calculating coeffs */
  get_Ylm_coeffs(realClm,imagClm,f,Ntheta,Nphi,lmax);
  
  df = 0.;
  /* now let's see how Ylm derivatives: works: */
  df_dphi = df_dphi_Ylm(realClm,imagClm,Ntheta,Nphi,lmax);
  printf("df(theta,phi)/dphi:\n");
  printf("--> For Ntheta = %u, Nphi = %u, Lmax = %u:\n",Ntheta,Nphi,lmax);
  /* for each theta and phi */
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      if (df < fabs(f_dphi[j+Nphi*i]-df_dphi[j+Nphi*i]))
        df = fabs(f_dphi[j+Nphi*i]-df_dphi[j+Nphi*i]);
    }
  }/* for (i = 0; i < Ntheta; ++i) */
  printf("Max Error = %e\n",df);
  
  df = 0.;
  /* now let's see how Ylm derivatives: works: */
  df_dtheta = df_dtheta_Ylm(realClm,imagClm,Ntheta,Nphi,lmax);
  printf("df(theta,phi)/dtheta:\n");
  printf("--> For Ntheta = %u, Nphi = %u, Lmax = %u:\n",Ntheta,Nphi,lmax);
  /* for each theta and phi */
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      if (df < fabs(f_dtheta[j+Nphi*i]-df_dtheta[j+Nphi*i]))
        df = fabs(f_dtheta[j+Nphi*i]-df_dtheta[j+Nphi*i]);
    }
  }/* for (i = 0; i < Ntheta; ++i) */
  printf("Max Error = %e\n",df);
  
  free(f);
  free(f_dphi);
  free(f_dtheta);
  free(df_dphi);
  free(df_dtheta);
  free(realClm);
  free(imagClm);
  UNUSED(grid);
  
  return TEST_SUCCESSFUL;
}


/* testing : r2cft_2d functions interpolations, derivatives, etc
// ->return value: TEST_SUCCESSFUL */
static int r2cft_2d_EquiSpaced_test(Grid_T *const grid)
{
  const Uint Nphi0 = 9;
  const Uint Nphi1 = 13;
  const Uint l0 = Nphi0/2+1;
  const Uint l1 = Nphi1/2+1;
  const Uint l0l1 = l0*l1;
  double *f = alloc_double(Nphi0*Nphi1);
  double *df_dx = alloc_double(Nphi0*Nphi1);
  double *df_dy = alloc_double(Nphi0*Nphi1);
  double *df_dphi0 = 0, *df_dphi1 = 0;
  double *realC = 0;
  double *imagC = 0; 
  double x,y;
  Uint i,j,m0,m1,ij;
  
  /* populate the values */
  for (i = 0; i < Nphi0; ++i)
  {
    x = 2.*i*M_PI/Nphi0;
    for (j = 0; j < Nphi1; ++j)
    {
      y = 2.*j*M_PI/Nphi1;
      
      f[IJ(i,j,Nphi1)] = 1 + Cos(x) - Power(Cos(x),2) + 
          Cos(y) + Sin(x) + Sin(4*x)/2. + Sin(y) + 
          Power(Sin(x),2)*Power(Sin(y),2) - Power(Sin(x) + Sin(y),2) - 
          Sin(6*y)/2.;
          
      df_dx[IJ(i,j,Nphi1)] = Cos(x) + 2*Cos(4*x) - 
            Sin(x) - 2*Cos(x)*Sin(y) + Sin(2*x)*Power(Sin(y),2);
      df_dy[IJ(i,j,Nphi1)] = Cos(y) - 3*Cos(6*y) - Sin(y) - 
            2*Cos(y)*(Sin(x) + Sin(y)) + Power(Sin(x),2)*Sin(2*y);


    }
  }
  
  /* calculating coeffs */
  r2cft_2d_coeffs(f,Nphi0,Nphi1,&realC,&imagC);
  
  /* print bases */
  if(1)
  {
    printf("Bases: Crr+Cri*I\n");
    for (m0 = 0; m0 < l0; ++m0)
    {
      for (m1 = 0; m1 < l1; ++m1)
      {
        Uint m0m1 = IJ(m0,m1,l1);
        printf("Crr(%u,%u)+Cri(%u,%u)I = %+.2f %+.2fI\n",
              m0,m1,m0,m1,realC[m0m1],imagC[m0m1]);
      }
    }
    
    printf("Bases: Cir+Cii*I\n");
    for (m0 = 0; m0 < l0; ++m0)
    {
      for (m1 = 0; m1 < l1; ++m1)
      {
        Uint m0m1 = IJ(m0,m1,l1);
        printf("Cir(%u,%u)+Cii(%u,%u)I = %+.2f %+.2fI\n",
              m0,m1,m0,m1,realC[l0l1+m0m1],imagC[l0l1+m0m1]);
      }
    }
  }
  /* let's do some interpolation: */
  double *ran_phi0 = make_random_number(Nphi0,0,2*M_PI);
  double *ran_phi1 = make_random_number(Nphi1,0,2*M_PI);
  if (1)
  {
    printf("Interpolation:       f(x):                diff:\n");
    for (i = 0; i < Nphi0; ++i)
    {
      x = ran_phi0[i];
      for (j = 0; j < Nphi1; ++j)
      {
        y = ran_phi1[j];
        double fr = 1 + Cos(x) - Power(Cos(x),2) + 
          Cos(y) + Sin(x) + Sin(4*x)/2. + Sin(y) + 
          Power(Sin(x),2)*Power(Sin(y),2) - Power(Sin(x) + Sin(y),2) - 
          Sin(6*y)/2.;
                    
        double fi = r2cft_2d_interpolation(realC,imagC,Nphi0,Nphi1,x,y);
        printf("%+0.15f   %+0.15f   %+e\n",fi,fr,fi-fr);
      }
    }
  }

  /* derivative tests: */
  df_dphi0 = r2cft_2d_df_dphi0(realC,imagC,Nphi0,Nphi1);
  df_dphi1 = r2cft_2d_df_dphi1(realC,imagC,Nphi0,Nphi1);
  
  printf("df(x)/dphi0|num:     df(x)/dphi0|anl:     diff:\n");
  for (i = 0; i < Nphi0; ++i)
  {
    for (j = 0; j < Nphi1; ++j)
    {
      ij=IJ(i,j,Nphi1);
      printf("%+0.15f   %+0.15f   %+e\n",
          df_dphi0[ij],df_dx[ij],df_dphi0[ij]-df_dx[ij]);
    }
  }
    
  printf("df(x)/dphi1|num:     df(x)/dphi1|anl:     diff:\n");
  for (i = 0; i < Nphi0; ++i)
  {
    for (j = 0; j < Nphi1; ++j)
    {
      ij=IJ(i,j,Nphi1);
      printf("%+0.15f   %+0.15f   %+e\n",
          df_dphi1[ij],df_dy[ij],df_dphi1[ij]-df_dy[ij]);
    }
  }
    
  free(f);
  free(realC);
  free(imagC);
  free(df_dphi0);
  free(df_dphi1);
  free(df_dx);
  free(df_dy);
  free(ran_phi0);
  free(ran_phi1);
  
  UNUSED(grid);
  return TEST_SUCCESSFUL;
}

/* testing : r2cft_2d functions interpolations, derivatives, etc on S2
// ->return value: TEST_SUCCESSFUL */
static int r2cft_2d_EquiSpaced_S2_test(Grid_T *const grid)
{
  const Uint Ntheta = 10;
  const Uint Nphi = 10;
  const Uint l0   = Ntheta;
  const Uint l1   = Nphi/2+1;
  const Uint l0l1 = l0*l1;
  double *f = alloc_double(Ntheta*Nphi);
  double *df_dx = alloc_double(Ntheta*Nphi);
  double *df_dy = alloc_double(Ntheta*Nphi);
  double *df_dtheta = 0, *df_dphi = 0;
  double *realC = 0;
  double *imagC = 0; 
  double x,y;
  Uint i,j,ij,m0,m1;
 
  /* populate the values */
  for (i = 0; i < Ntheta; ++i)
  {
    x = i*M_PI/(Ntheta-1);
    for (j = 0; j < Nphi; ++j)
    {
      y = 2.*j*M_PI/Nphi;
      
      f[IJ(i,j,Nphi)] =  (2*sin(2*x)*cos(y)+1)*cos(x);
           
      df_dx[IJ(i,j,Nphi)] = 4*Power(Cos(x),3)*Cos(y) - Sin(x) 
                            - 8*Cos(x)*Cos(y)*Power(Sin(x),2);
      df_dy[IJ(i,j,Nphi)] = -4*Power(Cos(x),2)*Sin(x)*Sin(y);


    }
  }
  
  /* calculating coeffs */
  r2cft_2d_coeffs_S2(f,Ntheta,Nphi,&realC,&imagC,1);
  /* print bases */
  if(1)
  {
    const double eps = 1E-6;
    printf("Large Bases: Crr+Cri*I\n");
    for (m0 = 0; m0 < l0; ++m0)
    {
      for (m1 = 0; m1 < l1; ++m1)
      {
        Uint m0m1 = IJ(m0,m1,l1);
        if ( GRT(fabs(realC[m0m1]),eps) ||
             GRT(fabs(imagC[m0m1]),eps) )
        printf("Crr(%u,%u)+Cri(%u,%u)I = %g %gI\n",
              m0,m1,m0,m1,realC[m0m1],imagC[m0m1]);
      }
    }
    
    printf("Bases: Cir+Cii*I\n");
    for (m0 = 0; m0 < l0; ++m0)
    {
      for (m1 = 0; m1 < l1; ++m1)
      {
        Uint m0m1 = IJ(m0,m1,l1);
        if ( GRT(fabs(realC[l0l1+m0m1]),eps) ||
             GRT(fabs(imagC[l0l1+m0m1]),eps))
        printf("Cir(%u,%u)+Cii(%u,%u)I = %+g %+gI\n",
              m0,m1,m0,m1,realC[l0l1+m0m1],imagC[l0l1+m0m1]);
      }
    }
  }
  /* let's do some interpolation: */
  double *ran_theta = make_random_number(Ntheta,0,M_PI);
  double *ran_phi = make_random_number(Nphi,0,2*M_PI);
  if (1)
  {
    printf("Interpolation:       f(x):                diff:\n");
    for (i = 0; i < Ntheta; ++i)
    {
      x = ran_theta[i];
      for (j = 0; j < Nphi; ++j)
      {
        y = ran_phi[j];
        double fr = (2*sin(2*x)*cos(y)+1)*cos(x);
        double fi = r2cft_2d_interpolation_S2(realC,imagC,Ntheta,Nphi,x,y);
        printf("%+0.15f   %+0.15f   %+e\n",fi,fr,fi-fr);
      }
    }
  }

  /* derivative tests: */
  if (1)
  {
    df_dtheta = r2cft_2d_df_dtheta_S2(realC,imagC,Ntheta,Nphi);
    df_dphi   = r2cft_2d_df_dphi_S2(realC,imagC,Ntheta,Nphi);
    
    printf("df(x)/dtheta|num:     df(x)/dtheta|anl:     diff:\n");
    for (i = 0; i < Ntheta; ++i)
    {
      for (j = 0; j < Nphi; ++j)
      {
        ij=IJ(i,j,Nphi);
        printf("%+0.15f   %+0.15f   %+e\n",
            df_dtheta[ij],df_dx[ij],df_dtheta[ij]-df_dx[ij]);
      }
    }
      
    printf("df(x)/dphi|num:     df(x)/dphi|anl:     diff:\n");
    for (i = 0; i < Ntheta; ++i)
    {
      for (j = 0; j < Nphi; ++j)
      {
        ij=IJ(i,j,Nphi);
        printf("%+0.15f   %+0.15f   %+e\n",
            df_dphi[ij],df_dy[ij],df_dphi[ij]-df_dy[ij]);
      }
    }
    free(df_dtheta);
    free(df_dphi);
  }  
  free(f);
  free(realC);
  free(imagC);
  free(df_dx);
  free(df_dy);
  free(ran_theta);
  free(ran_phi);
  
  UNUSED(grid);
  return TEST_SUCCESSFUL;
}

/* testing : r2cft_1d_EquiSpaced_coeffs function.
// ->return value: TEST_SUCCESSFUL */
static int cft_c2r_r2c_1d_EquiSpaced_test(Grid_T *const grid)
{
  Uint N = 0;
  const char *const par = Pgets("Test_FourierTransformation");
  double *f;/* f(x) x = [0,2pi] */
  double complex *c;/* f(x) = c_i*exp(imagI*i*x) */
  double x;
  Uint i,j;
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    N = (Uint)atoi(s);
    Free(s);
  }
  else
    N = 20;
 
  f = alloc_double_complex(N);
  for (i = 0; i < N; ++i)
  {
    x    = 2.*i*M_PI/N;
    f[i] = sin(2*x)*cos(x)+Pow2(sin(4*x));
  }
  /* calculating coeffs */
  c = r2cft_1d_EquiSpaced_coeffs(f,N);
  
  /* now let's see how Fourier sum works: */
  printf("%-*s%-*s diff:\n",40,"Fourier sum:",20,"f(x):");
  for (i = 0; i < N; ++i)
  {
    x = 2.*i*M_PI/N;
    double complex fc = 0;
    fc = c[0];
    for (j = 1; j < N/2+1; ++j)
      fc += c[j]*cexp(imagI*(double)j*x)+conj(c[j])*cexp(-imagI*(double)j*x);
    printf("%+0.15f%+0.15fI   %+0.15f   %+e\n",creal(fc),cimag(fc),f[i],creal(fc)-f[i]);
  }
  
  /* let's do some interpolation too: */
  double *rand = make_random_number(N,0,2*M_PI);
  printf("\nFor n = %u:\n",N);
  printf("%-*s%-*s diff:\n",40,"Interpolation:",20,"f(x):");
  for (i = 0; i < N; ++i)
  {
    x = rand[i];
    double complex fi = 0;
    double fr = sin(2*x)*cos(x)+Pow2(sin(4*x));
    fi = c[0];
    for (j = 1; j < N/2+1; ++j)
      fi += c[j]*cexp(imagI*(double)j*x)+conj(c[j])*cexp(-imagI*(double)j*x);
    printf("%+0.15f%+0.15fI   %+0.15f   %+e\n",creal(fi),cimag(fi),fr,creal(fi)-fr);
  }
  
  /* let's also check the inverse transformation: */
  double *f_inv = c2rft_1d_EquiSpaced_values(c,N);
  printf("\nFor n = %u:\n",N);
  printf("Checking inverse Fourier function:\n");
  printf("%-*s%-*s diff:\n",21,"Fourier sum:",20,"f(x):");
  for (i = 0; i < N; ++i)
    printf("%+0.15f   %+0.15f   %+e\n",f_inv[i],f[i],f_inv[i]-f[i]);
  
/* undefining I complex */   
#ifdef I
#undef I
#endif

  free(f);
  free(c);
  free(rand);
  free(f_inv);
  UNUSED(grid);
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation functions.
// takes bunch of random points inside each patch and compares
// the value of interpolation and exact value for a field.
// ->return value: EXIT_SUCCESS
*/
int interpolation_tests(Grid_T *const grid)
{
  FUNC_TIC
  Uint p;
  double *X,*Y,*Z;
  Field_T *field;
  int status;
  
  printf("Note: to test spectral interpolation, a fifth order \n"
         "      polynomial is used, which might not get well resolved in \n"
         "      outermost patches which use compactification.\n\n");
  
  FOR_ALL_PATCHES(p,grid) //Un-comment for full tests.
  {
  if (DO_NOT) { /////////////////Remove if(DO_NOT) to test functions other than
                /////////////////1-D interpolation
    Patch_T *patch = grid->patch[p];
    const Uint *n = patch->n;
    const double *min = patch->min;
    const double *max = patch->max;
    
    /* making a field which have analytic properties */
    field = add_field("interpolant","(3dim)",patch,NO);
    field->v = poly5_f(patch);/* a polynomial */
    
    /* making 2n random number in each direction */
    X = make_random_number(n[0],min[0],max[0]);
    Y = make_random_number(n[1],min[1],max[1]);
    Z = make_random_number(n[2],min[2],max[2]);
    
    if (DO)
    {
      printf("Interpolation test:      X direction, patch %10s:\n",patch->name);
      status = interpolation_tests_X(field,X,n[0]);
      check_test_result(status);
    }
    
    if (DO)
    {
      printf("Interpolation test:      Y direction, patch %10s:\n",patch->name);
      status = interpolation_tests_Y(field,Y,n[1]);
      check_test_result(status);
    }
    if (DO)
    {
      printf("Interpolation test:      Z direction, patch %10s:\n",patch->name);
      status = interpolation_tests_Z(field,Z,n[2]);
      check_test_result(status);
    }
    if (DO)
    {
      printf("Interpolation test: X & Y directions, patch %10s:\n",patch->name);
      status = interpolation_tests_XY(field,X,Y,n[0],n[1]);
      check_test_result(status);
    }
    if (DO)
    {
      printf("Interpolation test: X & Z directions, patch %10s:\n",patch->name);
      status = interpolation_tests_XZ(field,X,Z,n[0],n[2]);
      check_test_result(status);
    } 
    if (DO)
    {
      printf("Interpolation test: Y & Z directions, patch %10s:\n",patch->name);
      status = interpolation_tests_YZ(field,Y,Z,n[1],n[2]);
      check_test_result(status);
    }
    if (DO)
    {
      printf("Interpolation test:              3-D, patch %10s:\n",patch->name);
      status = interpolation_tests_XYZ(field,X,Y,Z,n[0],n[1],n[2]);
      check_test_result(status);
    }
    
    
    /* freeing */
    free(X);
    free(Y);
    free(Z);
    remove_field(field);
  }}////////////////////////////
  
  if (DO_NOT)
  {
      printf("Interpolation test:            Neville Iterative Method =>");
      status = interpolation_tests_Neville_1d();
      check_test_result(status);
  }
  if (DO_NOT)
  {
      printf("Interpolation test:            Natural Cubic Spline Method =>");
      status = interpolation_tests_N_cubic_spline_1d();
      check_test_result(status);
  }
  if (DO_NOT)
  {
      printf("Interpolation test:            1D Spline =>\n");
      status = interpolation_tests_spline_1d();
      check_test_result(status);
  }
  if (DO_NOT)
  {
      printf("Interpolation test:           Single Interpolant =>\n");
      status = interpolation_tests_single_interpolant();
      check_test_result(status);
  }
  if (DO_NOT)
  {
      printf("Interpolation test:            Fixed Point =>\n");
      status = interpolation_tests_fixed();
      check_test_result(status);
  }
  // Convergence tests for 1D interpolation
  // Edit parameters for number of trials and number of
  // spline knots in each trial here.
  if (DO_NOT)
  {
      printf("Interpolation test:             Convergence Test=>\n");
      printf("Interpolation method: %s\n", Pgets("Interpolation_Method"));
      status = interpolation_tests_Convergence(500, 30, 1000);
      check_test_result(status);
  }
  // Convergence test for Hermite spline order
  // (i.e. error vs order of interpolating polynomial).
  if (DO_NOT)
  {
      printf("Interpolation test:             Hermite Order Test=>\n");
      status = interpolation_tests_Hermite_Order();
      check_test_result(status);
  }
  
  if (DO)
  {
      printf("Interpolation test:             Hermite derivative test=>\n");
      status = interpolation_tests_Hermite_Derivative();
      check_test_result(status);
  }
  
  ///////////////////////////////Finite difference tests////////////////////////////
  if (DO_NOT)
  {
      srand((Uint)time(NULL));
      printf("Finite Difference test:         Fornberg and Uniform methods =>\n");
      status = interpolation_tests_FDM();
      check_test_result(status);
  }
  if (DO_NOT)
  {
      srand((Uint)time(NULL));
      printf("Finite Difference test:         Fornberg method (random grid) =>\n");
      status = interpolation_tests_Fornberg();
      check_test_result(status);
  }
  if (DO_NOT)
  {
      srand((Uint)time(NULL));
      printf("Finite difference test:         FDM Convergence =>\n");
      status = interpolation_tests_FDM_convergence();
      check_test_result(status);
  }
  if (DO_NOT)
  {
      printf("Finite difference test:         Fixed Point =>\n");
      status = interpolation_tests_FDM_fixed();
      check_test_result(status);
  }
  
  FUNC_TOC
  return EXIT_SUCCESS;
}
   
// Prints f(x) to text file for
// visual confirmation of interpolation tests.
static void print_arrays(const char *const fileName,
                         double* x, double* f, double N)
{
  const char *const path_par = Pgets("top_directory");
  char outName[1000]; //*path;
  //path = make_directory(path_par,"Test_Data");
  sprintf(outName, "%s/%s.txt", path_par, fileName);
  FILE* outFile = fopen(outName, "w");
  
  if (outFile == NULL)
  { printf("Error printing interpolation test data: could not open file.\n"); }
  
  for (Uint j = 0; j < N; j++)
  { fprintf(outFile, "%E\t%E\n", x[j], f[j]); }
  
  fclose(outFile);
}

// Test function used in interpolation and FDM tests.
static double test_fxn(double x)
{
  //return x;
  //return 3.0*x*x*x + 2.0*x*x - 4.0*x;
  //return x*x*x*x*x*x - 10.0*x*x*x*x*x + 20.0*x*x*x;
  return log(x) * cos(x) + x;
}

// First derivative
static double test_fxn_p(double x)
{
  //UNUSED(x);
  //return 1.0;
  //return 9.0*x*x + 4.0*x - 4.0;
  //return 6.0*x*x*x*x*x - 50.0*x*x*x*x + 60.0*x*x;
  return - log(x)*sin(x) + cos(x)/x + 1;
}

// Third derivative
static double test_fxn_ppp(double x)
{
  //UNUSED(x);
  //return 0.0;
  //return 18.0;
  //return 120.0*x*x*x - 600.0*x*x + 120.0;
  return (2*cos(x)/Pow3(x) + 3*sin(x)/Pow2(x)
        - 3*cos(x)/x + log(x)*sin(x));
}

// Wrapper which selects the desired derivative of test_fxn
static double test_fxn_derivative(double x, Uint derivative)
{
  double (*fxns[4])(double a);
  fxns[0] = test_fxn;
  fxns[1] = test_fxn_p;
  fxns[2] = NULL;
  fxns[3] = test_fxn_ppp;
  
  return fxns[derivative](x);
}

// Returns random integer between a and b.
static int randRange(int a, int b)
{ return a + rand() % (b - a + 1); }

// Returns random double between a and b.
static double rand_double(double a, double b)
{ return a + (rand() * (b-a) / RAND_MAX); }

/* test Natural Cubic Spline method for 1-d arrays.
// Generates spline knots from analytical test function,
// then attempts to interpolate points on random grid.
// Prints results to files.
// ->return value: result of test. */
static int interpolation_tests_N_cubic_spline_1d(void)
{
  Interpolation_T *interp_s = init_interpolation();
  const Uint N                    = (Uint)Pgeti("n_interp");
  double *f                       = alloc_double(N);
  double *f_derivative            = alloc_double(N);
  double *x                       = alloc_double(N);
  double *error                   = alloc_double(N);
  double *error_derivative        = alloc_double(N);
  double *interp_vals             = alloc_double(N);
  double *interp_vals_derivative  = alloc_double(N);
  double *x_vals                  = alloc_double(N);
  const double a = -M_PI, b = 3/4*M_PI;/* an arbitrary interval  */
  double *hs = make_random_number(N,a,b);
  double s = (b-a)/(N-1);
  double t,interp,interp_derivative;
  double time;
  double error_total = 0;
  double error_total_derivative = 0;
  clock_t time1 = clock();
  Flag_T flg = NONE;
  Uint i;
  
  for (i = 0; i < N; ++i)
  {
    t = x[i] = a+i*s;
    f[i] = cos(t)+t*t*t;/* arbitrary function */
    f_derivative[i] = - sin(t) + 3*t*t; // Derivative
  }
  
  interp_s->method         = "Natural_Cubic_Spline_1D";
  interp_s->N_cubic_spline_1d->f   = f;
  interp_s->N_cubic_spline_1d->x   = x;
  interp_s->N_cubic_spline_1d->N   = N;
  assign_interpolation_ptrs(interp_s);
  *interp_s->f   = f;
  *interp_s->x   = x;
  *interp_s->N   = N;
  plan_interpolation(interp_s);
  
  for (i = 0; i < N; ++i)
  {
    double diff;
    double diff_derivative;
    t = hs[i];
    x_vals[i] = hs[i];
    interp_s->N_cubic_spline_1d->h = t;
    interp = execute_interpolation(interp_s);
    interp_derivative = execute_derivative_interpolation(interp_s);
    diff = interp-(cos(t)+t*t*t);
    error[i] = Pow2(diff);
    diff_derivative = interp_derivative - (-sin(t) + 3*t*t);
    error_derivative[i] = Pow2(diff_derivative);
    
    interp_vals[i] = interp;
    interp_vals_derivative[i] = interp_derivative;
    
    if (GRT(fabs(diff),s))
    {
      fprintf(stderr,"diff = %g\n",diff);
      flg = FOUND;
      break;
    }
  }
  
  time = (double)(clock() - time1)/CLOCKS_PER_SEC;
  printf("\nNCS test time: %E\n", time);
  for (Uint j; j < N; j++)
  { 
    error_total += error[j];
    error_total_derivative += error_derivative[j];
  }
  printf("Total sum-of-squares error: %E\n", error_total);
  printf("Total sum-of-squares error for first derivative: %E\n",
        error_total_derivative);
  
  if (strstr_i(PgetsEZ("print_interpolation_tests"), "yes"))
  { 
    print_arrays("NCS_Function", x, f, N);
    print_arrays("NCS_Interp", x_vals, interp_vals, N);
    print_arrays("NCS_Error", x_vals, error, N);
    print_arrays("NCS_Derivative", x, f_derivative, N);
    print_arrays("NCS_Interp_Derivative", x_vals, interp_vals_derivative, N);
    print_arrays("NCS_Error_Derivative", x_vals, error_derivative, N);
  }
  
  free_interpolation(interp_s);
  /* let's test the reveres order for x's */
  s = -(b-a)/(N-1);
  interp_s = init_interpolation();
  for (i = 0; i < N; ++i)
  {
    t = x[i] = b+i*s;
    f[i] = cos(t)+t*t*t;/* arbitrary function */
  }
    
  interp_s->method         = "Natural_Cubic_Spline_1D";
  interp_s->N_cubic_spline_1d->f   = f;
  interp_s->N_cubic_spline_1d->x   = x;
  interp_s->N_cubic_spline_1d->N   = N;
  plan_interpolation(interp_s);
  
  for (i = 0; i < N; ++i)
  {
    double diff;
    t = hs[i];
    interp_s->N_cubic_spline_1d->h = t;
    interp = execute_interpolation(interp_s);
    diff = interp-(cos(t)+t*t*t);
    
    if (GRT(fabs(diff),-s))
    {
      fprintf(stderr,"diff = %g\n",diff);
      flg = FOUND;
      break;
    }
  }
  free_interpolation(interp_s);
  free(f);
  free(f_derivative);
  free(x);
  free(hs);
  free(x_vals);
  free(interp_vals);
  free(interp_vals_derivative);
  free(error);
  free(error_derivative);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  return TEST_SUCCESSFUL;
}

// Tests spline for 1d arrays
// Available for multiple interpolation methods.
// Generates spline knots from analytical test function,
// then interpolates over random grid.
// Outputs analytical function, interpolated function,
// and error to files.
// Returns: result of test.
static int interpolation_tests_spline_1d(void)
{
  Interpolation_T *interp_s = init_interpolation();
  const Uint N                        = (Uint)Pgeti("Interpolation_grid_points");
  const Uint test_pts                 = (Uint)Pgeti("Interpolation_test_points");
  double *f                           = alloc_double(N);        // Analytical spline knots
  double *x                           = alloc_double(N);        // Grid points
  double *error                       = alloc_double(test_pts); // Error vs x points
  double *interp_vals                 = alloc_double(test_pts); // Interpolated values
  double *x_vals                      = alloc_double(test_pts); // For test grid points.
  const double a = 1, b = 4; // an arbitrary interval
  double *hs    = make_random_number(test_pts,a,b);
  double s      = (b-a)/(N-1);
  double t,interp;
  double time;
  double error_total = 0; // Cumulative error
  clock_t time1 = clock();
  Flag_T flg = NONE;
  Uint i;
  
  // Generate grid and spline knots
  for (i = 0; i < N; ++i)
  {
    t = x[i] = a+i*s;
    f[i] = test_fxn(t);   // Analytical function
  }
  
  // Set up interpolation
  interp_s->method = Pgets("Interpolation_method");
  assign_interpolation_ptrs(interp_s);
  *interp_s->f   = f;
  *interp_s->x   = x;
  *interp_s->N   = N;
  plan_interpolation(interp_s);
  
  // Interpolate over random grid
  for (i = 0; i < test_pts; ++i)
  {
    double diff;
    t                   = hs[i];
    x_vals[i]           = hs[i];
    *interp_s->h        = t;
    interp              = execute_interpolation(interp_s);
    diff                = interp - test_fxn(t);
    error[i]            = ABSd(diff);
    interp_vals[i]      = interp;
    
    // Failure in case of really high error
    if (GRT(fabs(diff),s))
    {
      printf("Maximum error exceeded at t = %E.\n",t);
      printf("Function value: %E\n", test_fxn(t));
      printf("Interpolated value: %E\n", interp);
      printf("Difference: %E\n", diff);
      fprintf(stderr,"diff = %g\n",diff);
      flg = FOUND;
    }
  }
  
  time = (double)(clock() - time1)/CLOCKS_PER_SEC;
  printf("\nSpline test time: %E\n", time);
  
  for (Uint j=0; j < test_pts; j++)
  { error_total += error[j]; }
  printf("Cumulative error: %E\n", error_total);
  
  // Output arrays to files
  if (strstr_i(PgetsEZ("print_interpolation_tests"), "yes"))
  { 
    print_arrays("Analytical", x, f, N);
    print_arrays("Interpolated", x_vals, interp_vals, test_pts);
    print_arrays("Error", x_vals, error, test_pts);
  }
  
  free_interpolation(interp_s);
  free(hs);
  free(x_vals);
  free(interp_vals);
  free(error);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  return TEST_SUCCESSFUL;
}

static int interpolation_tests_Hermite_Derivative(void)
{
  printf("Check 1\n");////////
  Interpolation_T *interp_s = init_interpolation();
  const Uint N                        = (Uint)Pgeti("Interpolation_grid_points");
  const Uint test_pts                 = (Uint)Pgeti("Interpolation_test_points");
  double *f                           = alloc_double(N);        // Analytical spline knots
  double *f_p                         = alloc_double(N);        // Analytical derivative values
  double *x                           = alloc_double(N);        // Grid points
  double *error                       = alloc_double(test_pts); // Error vs x points
  double *interp_vals                 = alloc_double(test_pts); // Interpolated values
  double *x_vals                      = alloc_double(test_pts); // For test grid points.
  const double a = 1, b = 4; // an arbitrary interval
  double *hs    = make_random_number(test_pts,a,b);
  double s      = (b-a)/(N-1);
  double t,interp;
  double time;
  double error_total = 0; // Cumulative error
  clock_t time1 = clock();
  Flag_T flg = NONE;
  Uint i;
  
  printf("Check 2\n");////////
  // Generate grid and spline knots
  for (i = 0; i < N; ++i)
  {
    t = x[i] = a+i*s;
    f[i] = test_fxn(t);   // Analytical function
    f_p[i] = test_fxn_p(t); // Analytical derivative
  }
  
  printf("Check 3\n");////////
  // Set up interpolation
  interp_s->method = "Hermite";
  assign_interpolation_ptrs(interp_s);
  *interp_s->f   = f;
  *interp_s->x   = x;
  *interp_s->N   = N;
  printf("Check 4\n");////////
  plan_interpolation(interp_s);
  printf("Check 5\n");////////
  
  // Interpolate over random grid
  for (i = 0; i < test_pts; ++i)
  {
    double diff;
    t                   = hs[i];
    x_vals[i]           = hs[i];
    *interp_s->h        = t;
    interp              = execute_derivative_interpolation(interp_s);
    diff                = interp - test_fxn_p(t);
    error[i]            = ABSd(diff);
    interp_vals[i]      = interp;
    
    // Failure in case of really high error
    if (GRT(fabs(diff),s))
    {
      printf("Maximum error exceeded at t = %E.\n",t);
      printf("Function value: %E\n", test_fxn(t));
      printf("Interpolated value: %E\n", interp);
      printf("Difference: %E\n", diff);
      fprintf(stderr,"diff = %g\n",diff);
      flg = FOUND;
    }
  }
  
  time = (double)(clock() - time1)/CLOCKS_PER_SEC;
  printf("\nSpline test time: %E\n", time);
  
  for (Uint j=0; j < test_pts; j++)
  { error_total += error[j]; }
  printf("Cumulative error: %E\n", error_total);
  
  // Output arrays to files
  if (strstr_i(PgetsEZ("print_interpolation_tests"), "yes"))
  { 
    print_arrays("Analytical", x, f_p, N);
    print_arrays("Interpolated", x_vals, interp_vals, test_pts);
    print_arrays("Error", x_vals, error, test_pts);
  }
  
  free_interpolation(interp_s);
  free(f_p);
  free(hs);
  free(x_vals);
  free(interp_vals);
  free(error);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  return TEST_SUCCESSFUL;
}

// Test single spline interval.
// Useful for checking exact agreement between interpolant
// and low-degree polynomials.
static int interpolation_tests_single_interpolant(void)
{
  Interpolation_T* interp = init_interpolation();
  Uint N = 4;
  Uint test_pts = 100;
  double a = -2, b = 1;
  double* x_full = alloc_double(N);
  double* f = alloc_double(N);
  double* x_vals = alloc_double(test_pts);
  double* analytical_vals = alloc_double(test_pts);
  double* interp_vals = alloc_double(test_pts);
  double* errors = alloc_double(test_pts);
  Uint j;
  
  double dx = (b-a)/(double)(N-1);
  // Create knots
  for (j = 0; j < N; j++)
  { 
    x_full[j]  = a + j*dx;
    f[j] = test_fxn(x_full[j]);
  }
  
  interp->method = Pgets("Interpolation_method");
  assign_interpolation_ptrs(interp);
  *interp->f   = f;
  *interp->x   = x_full;
  *interp->N   = N;
  plan_interpolation(interp);
  
  for (j = 0; j < test_pts; j++)
  {
    x_vals[j] = rand_double(a, b);
    *interp->h = x_vals[j];
    interp_vals[j] = execute_interpolation(interp);
    analytical_vals[j] = test_fxn(x_vals[j]);
    errors[j] = interp_vals[j] - analytical_vals[j];
  }
  
  print_arrays("Analytical", x_vals, analytical_vals, test_pts);
  print_arrays("Interpolation", x_vals, interp_vals, test_pts);
  print_arrays("Error", x_vals, errors, test_pts);
  
  free_interpolation(interp);
  free(x_vals);
  free(interp_vals);
  free(analytical_vals);
  free(errors);
  return TEST_SUCCESSFUL;
}

// Tests interpolation error vs grid refinement.
// Each iteration, double the number of grid points.
// Track error on fixed middle points vs number of grid points
// to test convergence.
static int interpolation_tests_fixed(void)
{
  const Uint MaxIter = 9;
  double* K_vals = alloc_double(MaxIter-3);
  Uint* N_array = malloc((MaxIter-3)*sizeof(*N_array));
  //const Uint Nf = (Uint)pow(2, MaxIter) * N0;
  const Uint N0 = 9;
  const double a = 0.5, b = 1.5;
  const double c = a + 0.5*(b-a); // Fixed central point
  
  // Each iteration:
  // Generate uniformly spaced grid from a to b,
  // with the halfway point fixed, and double the
  // points of the last iteration/half the grid spacing.
  Uint N = N0;
  double interp_1 = 0;
  double interp_2 = 0;
  double interp_3 = 0;
  for (Uint j = 0; j < MaxIter; j++)
  {
    N = 2*N-1;
    double dx = (b-a)/(N-1);
    double* x = alloc_double(N);
    double* f = alloc_double(N);
    
    // Generate uniform grid
    for (Uint k = 0; k < N; k++)
    { 
      x[k] = a + k*dx;
      f[k] = test_fxn(x[k]);
    }
    
    Interpolation_T* interp = init_interpolation();
    interp->method = Pgets("interpolation_method");
    assign_interpolation_ptrs(interp);
    *interp->x = x;
    *interp->f = f;
    *interp->N = N;
    plan_interpolation(interp);
    
    // Compute difference from previous approximation
    // on more refined grid.
    if (j > 2)
    {
      K_vals[j-3] = (log(ABSd((interp_2 - interp_1) / 
                    (interp_3 - interp_2))) / log(2.0));
      N_array[j-3] = N;
    }
    
    *interp->h = c;
    interp_1 = interp_2;
    interp_2 = interp_3;
    interp_3 = execute_interpolation(interp);
    
    free_interpolation(interp);
  }
  
  // Print results
  printf("i \t | N \t | K\n");
  for (Uint j = 0; j < MaxIter-3; j++)
  { printf("%i \t | %i \t | %E\n", j, N_array[j], K_vals[j]); }
         
  free(K_vals);
  free(N_array);
  return TEST_SUCCESSFUL;
}
    
// Uniform-grid finite difference methods, only used for testing
// (e.g. in comparison with Fornberg or semi-fixed FDM).
// Use fixed coefficients and central finite difference.

// 1st derivative, 6th-order accuracy.
// Parameters:
//    x: Array of evenly-spaced grid points.
//    f: Array of analytical function values f[j] = f(x[j])
//    h: Interpolation point x[0] < h < x[last]
//    N: length of x and f arrays, i.e. number of total grid points.
static double FDM_Uniform_1_6(double* x, double* f, double h, Uint N)
{
  Uint i = 0;
 
  // Find interval 
  for (i = 0; i < N-1; i++)
  {
    if (GRTEQL(h,x[i]) && LSS(h,x[i+1]))
    { break; }
  }
  
  // Error if out of bounds
  if (i < 3 || i > N - 4)
  {
    printf("Uniform_FDM_1_6: Interval: %u\n", i);
    printf("Uniform_FDM_1_6: Total points: %u\n", N); 
    printf("First derivative at 6-th order accuracy requires " 
           "3 grid points on either side of target point.\n");
    Error0("Uniform FDM: out-of-bounds.");
  }
  
  double ret = 0;
  ret += -(1.0/60.0) * f[i-3];
  ret += (3.0/20.0) * f[i-2];
  ret += -(3.0/4.0) * f[i-1];
  ret += (3.0/4.0) * f[i+1];
  ret += -(3.0/20.0) * f[i+2];
  ret += (1.0/60.0) * f[i+3];
  ret /= (ABSd(x[i] - x[i-1]));
  
  return ret;
}

// 3rd derivative, 6th-order accuracy.
// Parameters:
//    x: Array of evenly-spaced grid points.
//    f: Array of analytical function values f[j] = f(x[j])
//    h: Interpolation point x[4] < h < x[last-4]
//    N: length of x and f arrays, i.e. number of total grid points.
// Only valid for h = x[j] for some 
static double FDM_Uniform_3_6(double* x, double* f, double h, Uint N)
{
  Uint i = 0;
 
  // Find spline interval 
  for (i = 0; i < N-1; ++i)
  {
    if (GRTEQL(h,x[i]) && LSS(h,x[i+1]))
    {
      //flg = FOUND;
      break;
    }
  }
  
  // Error if out of bounds
  if (i < 4 || i > N - 5)
  {
    printf("Uniform_FDM_3_6: Interval: %u\n", i);
    printf("Uniform_FDM_3_6: Total points: %u\n", N); 
    printf("Third derivative at 6-th order accuracy requires " 
           "4 grid points on either side of target point.\n");
    Error0("Uniform FDM: out-of-bounds.");
  }
  
  double ret = 0;
  ret += -(7.0/240.0)     * f[i-4];
  ret += (3.0/10.0)       * f[i-3];
  ret += -(169.0/120.0)   * f[i-2];
  ret += (61.0/30.0)      * f[i-1];
  ret += -(61.0/30.0)     * f[i+1];
  ret += (169.0/120.0)    * f[i+2];
  ret += -(3.0/10.0)      * f[i+3];
  ret += (7.0/240.0)      * f[i+4];
  ret /= Pow3(ABSd(x[i+1] - x[i]));
  
  return ret;
}

// Tests uniform-grid and Fornberg finite difference methods
// against an analytical function by randomly sampling points 
// from a uniform grid.
// Plots approximated first and third derivatives of function,
// and associated error vs analytical value.
static int interpolation_tests_FDM(void)
{
  Uint N = (Uint)Pgeti("n_interp");
  Uint samples = N-2*4;
  double* x                 = alloc_double(N);        // Uniform grid
  double* f                 = alloc_double(N);        // Analytical values
  double* f_p               = alloc_double(N);        // Analytical values for 1st derivative
  double* f_ppp             = alloc_double(N);        // Analytical values for 3rd derivative
  double* x_samples         = alloc_double(samples);  // Test points
  double* f_p_uniform       = alloc_double(samples);  // Estimated 1st derivative points via uniform method
  double* f_p_Fornberg      = alloc_double(samples);  // Estimated 1st derivative points via Fornberg method
  double* f_ppp_uniform     = alloc_double(samples);  // Estimated 3rd derivative points via uniform method
  double* f_ppp_Fornberg    = alloc_double(samples);  // Estimated 3rd derivative points via Fornberg method
  double* err_p_uniform     = alloc_double(samples);  // Error for 1st derivative via uniform method
  double* err_p_Fornberg    = alloc_double(samples);  // Error for 1st derivative via Fornberg method
  double* err_ppp_uniform   = alloc_double(samples);  // Error for 3rd derivative via uniform method
  double* err_ppp_Fornberg  = alloc_double(samples);  // Error for 3rd derivative via Fornberg method
  double* difference        = alloc_double(samples);  // 1st-derivative differences between uniform and Fornberg methods
  double* difference_p      = alloc_double(samples);  // 3rd-derivative differences between uniform and Fornberg methods
  double a = 1, b = 4;                                // Grid end points
  Uint j;
  
  // Generate uniform grid from analytical function.
  double dx = (b-a)/N;
  for (j = 0; j < N; j++)
  {
    x[j]       = a + j*dx;
    f[j]       = test_fxn(x[j]);
    f_p[j]     = test_fxn_p(x[j]);
    f_ppp[j]   = test_fxn_ppp(x[j]);
  }
  
  // Randomly sample grid points
  // Note: The Fornberg method works anywhere on the x domain,
  // but is restricted here to the sample points for comparison
  // against the uniform-grid method. 
  for (j = 0; j < samples; j++)
  {
    //Uint r = (Uint)randRange(5, (int)N-5);
    UNUSED(randRange);
    x_samples[j]        = x[j+4];
    f_p_uniform[j]      = FDM_Uniform_1_6(x, f, x_samples[j], N);
    f_p_Fornberg[j]     = finite_difference_Fornberg(x, f, x_samples[j], 1, 6, N);
    f_ppp_uniform[j]    = FDM_Uniform_3_6(x, f, x_samples[j], N);
    f_ppp_Fornberg[j]   = finite_difference_Fornberg(x, f, x_samples[j], 3, 6, N);
    err_p_uniform[j]    = ABSd(f_p_uniform[j] - test_fxn_p(x_samples[j]));
    err_p_Fornberg[j]   = ABSd(f_p_Fornberg[j] - test_fxn_p(x_samples[j]));
    err_ppp_uniform[j]  = ABSd(f_ppp_uniform[j] - test_fxn_ppp(x_samples[j]));
    err_ppp_Fornberg[j] = ABSd(f_ppp_Fornberg[j] - test_fxn_ppp(x_samples[j]));
    difference[j]       = f_p_Fornberg[j] - f_p_uniform[j];
    difference_p[j]     = f_ppp_Fornberg[j] - f_ppp_uniform[j];
  }
  
  // Calculates cumulative error scaled by number of sample points.
  double total_err_p_uniform = 0;
  double total_err_p_Fornberg = 0;
  double total_err_ppp_uniform = 0;
  double total_err_ppp_Fornberg = 0;
  for (j = 0; j < samples; j++)
  {
    total_err_p_uniform += err_p_uniform[j];
    total_err_p_Fornberg += err_p_Fornberg[j];
    total_err_ppp_uniform += err_ppp_uniform[j];
    total_err_ppp_Fornberg += err_ppp_Fornberg[j];
  }
  total_err_p_uniform /= samples;
  total_err_p_Fornberg /= samples;
  total_err_ppp_uniform /= samples;
  total_err_ppp_Fornberg /= samples;
  
  // Write results to files.
  print_arrays("Analytical_f", x, f, N);
  print_arrays("Analytical_f_p", x, f_p, N);
  print_arrays("Analytical_f_ppp", x, f_ppp, N);
  print_arrays("Uniform_f_p", x_samples, f_p_uniform, samples);
  print_arrays("Fornberg_f_p", x_samples, f_p_Fornberg, samples);
  print_arrays("Uniform_f_ppp", x_samples, f_ppp_uniform, samples);
  print_arrays("Fornberg_f_ppp", x_samples, f_ppp_Fornberg, samples);
  print_arrays("Error_uniform_f_p", x_samples, err_p_uniform, samples);
  print_arrays("Error_Fornberg_f_p", x_samples, err_p_Fornberg, samples);
  print_arrays("Error_uniform_f_ppp", x_samples, err_ppp_uniform, samples);
  print_arrays("Error_Fornberg_f_ppp", x_samples, err_ppp_Fornberg, samples);
  print_arrays("Difference", x_samples, difference, samples);
  print_arrays("Difference_p", x_samples, difference_p, samples);
  
  printf("First derivative uniform method error: %E\n", total_err_p_uniform);
  printf("First derivative Fornberg method error: %E\n", total_err_p_Fornberg);
  printf("Third derivative uniform method error: %E\n", total_err_ppp_uniform);
  printf("Third derivative Fornberg method error: %E\n", total_err_ppp_Fornberg);
  
  free(x);
  free(f);
  free(f_p);
  free(f_ppp);
  free(x_samples);
  free(f_p_uniform);
  free(f_p_Fornberg);
  free(f_ppp_uniform);
  free(f_ppp_Fornberg);
  free(err_p_uniform);
  free(err_p_Fornberg);
  free(err_ppp_uniform);
  free(err_ppp_Fornberg);
  free(difference);
  free(difference_p);
  
  return TEST_SUCCESSFUL;
}

// Needed solely for qsort
static int qsort_compare(const void* a, const void* b)
{ return (*(const double*)a == *(const double*)b ? 0 :
         (*(const double*)a > *(const double*)b ? 1 : -1)); }

// Tests Fornberg method on non-uniform grid, off grid points.
// Similar to interpolation, the Fornberg method can estimate derivatives
// for points in between grid points, i.e. points h such that x[j] < h < x[j+1],
// in contrast to other finite difference methods that are only valid on
// grid points x[j] = h.
// The Fornberg method also works on non-uniform grid spacings, which is also
// tested here.
static int interpolation_tests_Fornberg(void)
{
  Uint N            = (Uint)Pgeti("FDM_grid_points");
  Uint test_pts     = (Uint)Pgeti("FDM_test_points");
  Uint n            = (Uint)Pgeti("fd_accuracy_order");
  Uint derivative   = 1; // Degree of derivative
  double a = 1, b = 4;
  double* x         = alloc_double(N);
  double* f         = alloc_double(N);
  double* f_p       = alloc_double(N);
  double* x_vals    = alloc_double(test_pts);
  double* FDM_vals  = alloc_double(test_pts);
  double* error     = alloc_double(test_pts);
  double avg_error  = 0;
  Uint j;
  
  // Fill non-uniform analytical grid.
  for (j = 0; j < N; j++)
  { x[j] = rand_double(a, b); }
  qsort(x, (size_t)N, sizeof(double), qsort_compare);
  x[0] = a;
  x[N-1] = b;
  for (j = 0; j < N; j++)
  {
    f[j] = test_fxn(x[j]);
    f_p[j] = test_fxn_derivative(x[j], derivative);
  }
  
  // Test points at random and calculate error vs analytical function.
  for (j = 0; j < test_pts; j++)
  {
    x_vals[j] = rand_double(a, b);
    FDM_vals[j] = finite_difference_Fornberg(x, f, x_vals[j], derivative, n, N);
    error[j] = FDM_vals[j] - test_fxn_derivative(x_vals[j], derivative);
    avg_error += ABSd(error[j]);
  }
  avg_error /= (double)test_pts;
  
  // Write arrays to files.
  print_arrays("Analytical_Function", x, f, N);
  print_arrays("Analytical_Derivative", x, f_p, N);
  print_arrays("FDM_Values", x_vals, FDM_vals, test_pts);
  print_arrays("Error", x_vals, error, test_pts);
  printf("\t Grid points: %i\n", N);
  printf("\t Test points: %i\n", test_pts);
  printf("\t Average error: %E\n", avg_error);

  free(x);
  free(f);
  free(f_p);
  free(x_vals);
  free(FDM_vals);
  free(error);
    
  return TEST_SUCCESSFUL;
} 

// Tests error convergence (error vs number of grid points)
// for Fornberg and uniform-grid finite difference methods.
// For each trial, generate a uniform grid of
// analytical function values.
// Then calculate derivative values using finite difference
// methods on a random sample of grid points, and calculate
// error from analytical value.
// Find scaled cumulative error from each trial, and plot
// error vs number of grid points. 
static int interpolation_tests_FDM_convergence(void)
{
  Uint trials = 100;
  Uint samples;
  Uint N0 = 20, Nf = 1000;
  Uint dN = (Uint)ceil(((double)Nf-(double)N0)/(double)trials);
  double a = 1, b = 4;
  double* N_array = alloc_double(trials);
  double* errors_p_uniform = alloc_double(trials);    // 1st-derivative error for uniform method vs # grid points
  double* errors_p_Fornberg = alloc_double(trials);   // 1st-derivative error for Fornberg method vs # grid points
  double* errors_ppp_uniform = alloc_double(trials);  // 3rd-derivative error for uniform method vs # gird points
  double* errors_ppp_Fornberg = alloc_double(trials); // 3rd-derivative error for Fornberg method vs # grid points

  Uint N;
  double dx;
  double* x;                // Grid values
  double* f;                // Analytical values f[j] = f(x[j])
  double* f_p;              // Analytical values for 1st derivative
  double* f_ppp;            // Analytical values for 3rd derivative
  double* x_samples;        // Grid samples for testing
  double* f_p_uniform;      // 1st-derivative estimates via uniform method
  double* f_p_Fornberg;     // 1st-derivative estimates via Fornberg method
  double* f_ppp_uniform;    // 3rd-derivative estimates via uniform method
  double* f_ppp_Fornberg;   // 3rd-derivative 
  double err_p_uniform;     // Cumulative error for 1st derivative via uniform method
  double err_p_Fornberg;    // Cumulative error for 1st derivative via Fornberg method
  double err_ppp_uniform;   // Cumulative error for 3rd derivative via uniform method
  double err_ppp_Fornberg;  // Cumulative error for 3rd derivative via Fornberg method
  
  Uint j, k, trial = 0;
  for (N = N0; N <= Nf; N += dN)
  {
    samples = N;
    x                 = alloc_double(N);
    f                 = alloc_double(N);
    f_p               = alloc_double(N);
    f_ppp             = alloc_double(N);
    x_samples         = alloc_double(samples);
    f_p_uniform       = alloc_double(samples);
    f_p_Fornberg      = alloc_double(samples);
    f_ppp_uniform     = alloc_double(samples);
    f_ppp_Fornberg    = alloc_double(samples);
    err_p_uniform     = 0;
    err_p_Fornberg    = 0;
    err_ppp_uniform   = 0;
    err_ppp_Fornberg  = 0;
    
    // Generate uniform grid
    dx = (b - a) / N;
    for (j = 0; j < N; j++)
    {
      x[j]      = a + j*dx;
      f[j]      = test_fxn(x[j]);
      f_p[j]    = test_fxn_p(x[j]);
      f_ppp[j]  = test_fxn_ppp(x[j]);
    }
    
    // Sample grid and compute error vs analytical
    for (j = 5; j < samples-5; j++)
    {
      //k = randRange(5, N-5);
      k = j;
      x_samples[j]        = x[k];
      f_p_uniform[j]      = FDM_Uniform_1_6(x, f, x[k], N);
      f_p_Fornberg[j]     = finite_difference_Fornberg(x, f, x[k], 1, 6, N);
      f_ppp_uniform[j]    = FDM_Uniform_3_6(x, f, x[k], N);
      f_ppp_Fornberg[j]   = finite_difference_Fornberg(x, f, x[k], 3, 6, N);
      err_p_uniform       += ABSd(f_p_uniform[j] - test_fxn_p(x[k]));
      err_p_Fornberg      += ABSd(f_p_Fornberg[j] - test_fxn_p(x[k]));
      err_ppp_uniform     += ABSd(f_ppp_uniform[j] - test_fxn_ppp(x[k]));
      err_ppp_Fornberg    += ABSd(f_ppp_Fornberg[j] - test_fxn_ppp(x[k]));
    }
    err_p_uniform     /= samples;
    err_p_Fornberg    /= samples;
    err_ppp_uniform   /= samples;
    err_ppp_Fornberg  /= samples;
    
    N_array[trial] = (double)N;
    errors_p_uniform[trial] = err_p_uniform;
    errors_p_Fornberg[trial] = err_p_Fornberg;
    errors_ppp_uniform[trial] = err_ppp_uniform;
    errors_ppp_Fornberg[trial] = err_ppp_Fornberg;
    
    free(x);
    free(f);
    free(f_p);
    free(f_ppp);
    free(x_samples);
    free(f_p_uniform);
    free(f_p_Fornberg);
    free(f_ppp_uniform);
    free(f_ppp_Fornberg);
    
    trial++;
  }
  
  // Remove the junk end point
  N_array[trials-1] = N_array[trials-2];
  errors_p_uniform[trials-1] = errors_p_uniform[trials-2];
  errors_p_Fornberg[trials-1] = errors_p_Fornberg[trials-2];
  errors_ppp_uniform[trials-1] = errors_ppp_uniform[trials-2];
  errors_ppp_Fornberg[trials-1] = errors_ppp_Fornberg[trials-2];
  
  // Write arrays to data files
  print_arrays("Error_Uniform_p", N_array, errors_p_uniform, trials);
  print_arrays("Error_Fornberg_p", N_array, errors_p_Fornberg, trials);
  print_arrays("Error_Uniform_ppp", N_array, errors_ppp_uniform, trials);
  print_arrays("Error_Fornberg_ppp", N_array, errors_ppp_Fornberg, trials);
  
  free(N_array);
  free(errors_p_uniform);
  free(errors_p_Fornberg);
  free(errors_ppp_uniform);
  free(errors_ppp_Fornberg);
  
  return TEST_SUCCESSFUL;
}

// Tests error vs grid refinement on fixed point
// for uniform and Fornberg finite difference methods.
// Evenly spaced grid.
static int interpolation_tests_FDM_fixed(void)
{
  const Uint MaxIter = 9;
  double* K_uniform = alloc_double(MaxIter-3);
  double* K_Fornberg = alloc_double(MaxIter-3);
  Uint* N_array = malloc((MaxIter-3)*sizeof(*N_array));
  const Uint N0 = 9;
  const double a = 0, b = 2;
  const double c = a + 0.5*(b-a); // Fixed central point
  Uint order = 3; //Order of approximation
  
  // Each iteration:
  // Generate uniformly spaced grid from a to b,
  // with the halfway point fixed, and double the
  // points of the last iteration/half the grid spacing.
  Uint N = N0;
  double f_p_uniform_1 = 0;
  double f_p_uniform_2 = 0;
  double f_p_uniform_3 = 0;
  double f_p_Fornberg_1 = 0;
  double f_p_Fornberg_2 = 0;
  double f_p_Fornberg_3 = 0;
  for (Uint j = 0; j < MaxIter; j++)
  {
    N = 2*N-1;
    double dx = (b-a)/(N-1);
    double* x = alloc_double(N);
    double* f = alloc_double(N);
    
    // Generate uniform grid
    for (Uint k = 0; k < N; k++)
    { 
      x[k] = a + k*dx;
      f[k] = test_fxn(x[k]);
    }
    
    // Compute difference from previous approximation
    // on more refined grid.
    if (j > 2)
    {
      K_uniform[j-3] = (log(ABSd((f_p_uniform_2 - f_p_uniform_1) / 
                          (f_p_uniform_3 - f_p_uniform_2))) / log(2.0));
      K_Fornberg[j-3] = (log(ABSd((f_p_Fornberg_2 - f_p_Fornberg_1) /
                           (f_p_Fornberg_3 - f_p_Fornberg_2))) / log(2.0));
      N_array[j-3] = N;
    }
    f_p_uniform_1 = f_p_uniform_2;
    f_p_Fornberg_1 = f_p_Fornberg_2;
    f_p_uniform_2 = f_p_uniform_3;
    f_p_Fornberg_2 = f_p_Fornberg_3;
    f_p_uniform_3 = FDM_Uniform_1_6(x, f, c, N);
    f_p_Fornberg_3 = finite_difference_Fornberg(x, f, c, 1, order, N);
    
    free(x);
    free(f);
  }
  
  // Print results
  printf("FDM order of accuracy: %i\n", order);
  for (Uint j = 0; j < MaxIter-3; j++)
  { printf("N[%i] == %i\n", j, N_array[j]); }
  for (Uint j = 0; j < MaxIter-3; j++)
  { printf("Fornberg K[%i] == %E\n", j, K_Fornberg[j]); }
  for (Uint j = 0; j < MaxIter-3; j++)
  { printf("Uniform K[%i] == %E\n", j, K_uniform[j]); }
  printf("For finite difference methods with order of accuracy n,\n"
         "K should converge to the constant n.\n");
         
  free(K_uniform);
  free(K_Fornberg);
  free(N_array);
  return TEST_SUCCESSFUL;
}
    
/* test Neville iterative method for 1-d arrays.
// ->return value: result of test. */
static int interpolation_tests_Neville_1d(void)
{
  Interpolation_T *interp_s = init_interpolation();
  const Uint N = (Uint)Pgeti("n_a");
  double *f = alloc_double(N);
  double *x = alloc_double(N);
  const double a = -M_PI, b = 3/4*M_PI;/* an arbitrary interval  */
  double *hs = make_random_number(N,a,b);
  const double s = (b-a)/(N-1);
  double t,interp;
  Flag_T flg = NONE;
  Uint i;
  
  for (i = 0; i < N; ++i)
  {
    t = x[i] = a+i*s;
    f[i] = cos(t)+t*t*t;/* arbitrary function */
  }
    
  interp_s->method         = "Neville_1D";
  interp_s->Neville_1d->f   = f;
  interp_s->Neville_1d->x   = x;
  interp_s->Neville_1d->N   = N;
  interp_s->Neville_1d->max = N/2;
  plan_interpolation(interp_s);
  
  for (i = 0; i < N; ++i)
  {
    double diff;
    t = hs[i];
    interp_s->Neville_1d->h = t;
    interp = execute_interpolation(interp_s);
    diff = interp-(cos(t)+t*t*t);
    
    if (GRT(fabs(diff),s))
    {
      fprintf(stderr,"diff = %g\n",diff);
      flg = FOUND;
      break;
    }
  }
  free_interpolation(interp_s);
  free(f);
  free(x);
  free(hs);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
    
  return TEST_SUCCESSFUL;
}

// Tests error versus number of interpolated points
// for 1D interpolation methods.
// Interpolation methods must be specified in the parameter file,
// as well as interpolation parameters (e.g. derivative type).
static int interpolation_tests_Convergence(const Uint tests, const Uint N0, const Uint Nf)
{
  if (Nf <= N0 || N0 == 0 || Nf == 0 || tests == 0)
  { Error0("Error with interpolation convergence test parameters.\n"); }
  
  Flag_T flg = NONE;
  double* error_vals            = alloc_double(tests);
  double* error_derivative_vals = alloc_double(tests);
  double* times                 = alloc_double(tests);
  double* N_array               = alloc_double(tests);
  double time_total = 0;
  Uint N = N0;
  Uint delta_N = (Uint)floor((Nf-N0) / tests);
  for (Uint test = 0; test < tests; test++)
  {
    Interpolation_T *interp_s = init_interpolation();
    double *f                      = alloc_double(N);
    double *f_derivative           = alloc_double(N);
    double *x                      = alloc_double(N);
    double *error                  = alloc_double(N);
    double *error_derivative       = alloc_double(N);
    double *interp_vals            = alloc_double(N);
    double *interp_vals_derivative = alloc_double(N);
    double *x_vals                 = alloc_double(N);
    double *hs                     = alloc_double(N);
    const double a = 1, b = 4;// an arbitrary interval
    //double *hs = make_random_number(N,a,b);
    double s = (b-a)/(N-1);
    double t,interp,interp_derivative;
    double error_RMS = 0;
    double error_RMS_derivative = 0;
    clock_t time1 = clock();
    Uint i;
    
    for (i = 0; i < N; ++i)
    {
      t = x[i] = a+i*s;
      f[i] = cos(t)+t*t*t;                // arbitrary analytical function
      f_derivative[i] = - sin(t) + 3*t*t; // Derivative
      hs[i] = rand_double(a,b);           // Make random grid
    }
    
    interp_s->method = Pgets("Interpolation_Method");
    assign_interpolation_ptrs(interp_s);
    *interp_s->f   = f;
    *interp_s->x   = x;
    *interp_s->N   = N;
    plan_interpolation(interp_s);
    
    Uint cutoff = 6; // # knots we do not attempt to interpolate
    for (i = cutoff; i < N-cutoff; ++i)
    {
      double diff;
      double diff_derivative;
      t = hs[i];
      x_vals[i] = hs[i];
      *interp_s->h = t;
      interp = execute_interpolation(interp_s);
      interp_derivative = execute_derivative_interpolation(interp_s);
      diff = interp-(cos(t)+t*t*t);
      error[i] = Pow2(diff);
      diff_derivative = interp_derivative - (-sin(t) + 3*t*t);
      error_derivative[i] = Pow2(diff_derivative);
      
      interp_vals[i] = interp;
      interp_vals_derivative[i] = interp_derivative;
      
      // Cutoff due to high error is removed because trials with low
      // # knots will typically have high error.
      /*
      if (GRT(fabs(diff),s))
      {
        flg = FOUND;
        break;
      }
      */
    }
    
    //printf("\nTest time: %E\n", time);
    for (Uint j = cutoff; j < N-cutoff; j++)
    { 
      error_RMS += Pow2(error[j]);
      error_RMS_derivative += Pow2(error_derivative[j]);
    }
    error_RMS /= N-2*cutoff;
    error_RMS_derivative /= N-2*cutoff;
    
    error_vals[test] = error_RMS;
    error_derivative_vals[test] = error_RMS_derivative;
    times[test] = (double)(clock() - time1)/CLOCKS_PER_SEC;
    N_array[test] = (double)N;
    N += delta_N;
    
    free(f_derivative);
    free(hs);
    free(x_vals);
    free(interp_vals);
    free(interp_vals_derivative);
    free(error);
    free(error_derivative);
    free(f);
    free(x);
    set_interp_alloc_mem_flag(interp_s, 0);
    free_interpolation(interp_s);
  }
  
  for (Uint test = 0; test < tests; test++)
  { time_total += times[test]; }
  printf("Total time: %E\n", time_total);
  
  print_arrays("Error_Values", N_array, error_vals, tests);
  print_arrays("Derivative_Error_Values",
                N_array, error_derivative_vals, tests);
  print_arrays("Times", N_array, times, tests);
  
  for (Uint test = 0; test < tests; test++)
  { time_total += times[test]; }
  printf("Total time: %E\n", time_total);
  
  free(error_vals);
  free(error_derivative_vals);
  free(times);
  free(N_array);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
    
  return TEST_SUCCESSFUL;
}

// Tests convergence vs order of Hermite spline approximation.
// To run this test, "interpolations.h" and "interpolations.c"
// must be included in spectral_tests.c.
// When not running, comment out the interpolation planning block.
static int interpolation_tests_Hermite_Order(void)
{
  Uint min_order = 3;
  Uint max_order = 10;
  Uint tests = max_order - min_order + 1;
  Uint N = 200;
  Uint test_pts = 100;
  double* scaled_errors = alloc_double(tests);
  double* times = alloc_double(tests);
  double* test_vals = alloc_double(tests);
  double error;
  clock_t time1;
  double time_total = 0;
  Uint test;
  Uint j, k;
  Interpolation_T* interp;
  double* x_vals = alloc_double(test_pts);
  double* f_vals = alloc_double(test_pts);
  
  double* f = alloc_double(N);
  double* x = alloc_double(N);
  double x0 = 0;
  double xf = 3.14159;
  double dx = (xf - x0)/N;
  double* hs = make_random_number(N,x0 + max_order*dx, xf - max_order*dx);
  for (j = 0; j < N; j++)
  {
    x[j] = x0 + j*dx;
    f[j] = sin(x[j])*x[j];
  }
  
  k = 0;
  for (test = min_order; test <= max_order; test++)
  {
    time1 = clock();
    interp = init_interpolation();
    interp->method = "Hermite";
    assign_interpolation_ptrs(interp);
    *interp->x = x;
    *interp->f = f;
    *interp->N = N;
    
    /*
    ////////////////////COMMENT OUT WHEN NOT IN USE///////////////////////////
    // We have to manually plan the interpolation to override the parameter file.
    interp->spline_order = test;
    interp->Hermite_1d->spline_order = test;
    order_arrays_Hermite_1d(interp);
    interp->fd_accuracy_order = (Uint)Pgeti("Interpolation_finite_diff_order");
    find_coeffs_Hermite_1d(interp);
    interp->interpolation_func = interpolation_Hermite_1d;
    interp->interpolation_derivative_func = derivative_Hermite_1d;
    ///////////////////////////////////////////////////////////////////////////////
    */
    
    error = 0;
    for (j = 0; j < test_pts; j++)
    {
      interp->Hermite_1d->h = hs[j];
      *interp->h = hs[j];
      
      if (test == 5)
      {
        f_vals[j] = execute_interpolation(interp);
        x_vals[j] = hs[j];
      }
      error += fabs(execute_interpolation(interp) - sin(hs[j])*hs[j]);
    }
    scaled_errors[k] = ABSd(error / (double)test_pts);
    test_vals[k] = (double)test;
    times[k] = (double)(clock() - time1)/CLOCKS_PER_SEC;
    time_total += times[k];
    
    free(interp->Hermite_1d->a);
    free(interp->Hermite_1d->fp);
    free(interp);
    k++;
  }
  
  printf("Hermite Order Test:\n");
  for (test = 0; test < max_order - min_order; test++)
  {
    printf("\nTest %i:\n", test);
    printf("Spline order: %i\n", test + min_order);
    printf("Scaled error: %E\n", scaled_errors[test]);
    printf("Time: %E\n", times[test]);
  }
  
  print_arrays("5th_Order",x_vals,f_vals,test_pts);
  print_arrays("Analytical", x, f, N);
  print_arrays("Scaled_Error", test_vals, scaled_errors, max_order-min_order+1);
  printf("\nTotal time: %E\n", time_total);
  free(scaled_errors);
  free(test_vals);
  free(times);
  free(f);
  free(x);
  free(hs);
  free(x_vals);
  free(f_vals);
  return TEST_SUCCESSFUL;
}

/* testing interpolation in X direction and comparing
// to analytical value and returning the result.
// ->return value: result of test. */
static int interpolation_tests_X(Field_T *const field,const double *const X,const Uint N)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  const double max_f = L_inf(patch->nn,field->v);
  const double tol = n[0]*n[0]*1e-14*(max_f > 1 ? max_f : 1.);
  double Y,Z;
  double xc[3];/* Cartesian coord */
  double Xc[3];/* Curvilinear coord */
  double max = 0;
  Flag_T flg;
  Uint i,j;
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->X_dir_flag = 1;
  plan_interpolation(interp_s);
  
  //////////
  //printf("\ninterpolation_tests_x:\n");
  //printf("\tpatch total nodes: %i\n", patch->nn);
  //Uint l;
  //FOR_ALL_POINTS(l, patch)
  //{ printf("\t\tfield->v[%i] == %E\n", l, field->v[l]); }
  //////////
  
  
  flg = NONE;
  for (j = 0; j < N; ++j)/* -> choose different slices for testing */
  {
    //printf("\nj == %i\n", j);//////////
    interp_s->J = (Uint) floor(random_double(0,n[1],j));
    interp_s->K = (Uint) floor(random_double(0,n[2],1));
    
    Y = node[i_j_k_to_ijk(n,0,interp_s->J,0)]->X[1];
    Z = node[i_j_k_to_ijk(n,0,0,interp_s->K)]->X[2];
    Xc[1] = Y;
    Xc[2] = Z;
    
    for (i = 0; i < N; ++i)
    {
      interp_s->X = X[i];
      Xc[0] = X[i];
      x_of_X(xc,Xc,patch);
      
      //diff = poly5_f_point(xc[0],xc[1],xc[2])-execute_interpolation(interp_s);
      double a = execute_interpolation(interp_s);//////////
      diff = poly5_f_point(xc[0],xc[1],xc[2])-a;//////////
      //////////
      //printf("\ni == %i\n", i);
      //printf("\tCoords: (%.4E, %.4E, %.4E)\n", xc[0], xc[1], xc[2]);
      //printf("\tInterpolation: %E\n", a);
      //printf("\tDiff: %E\n", diff);
      //printf("----------------------\n");
      if (GRT(fabs(diff),tol))
      {
        flg = FOUND; max = (max<fabs(diff)?fabs(diff):max);
      }
    }
  }
  
  free_interpolation(interp_s);
  
  if (flg == FOUND)
  {printf(Pretty0"Max difference = %e\n",max); return TEST_UNSUCCESSFUL;}
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in Y direction and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_Y(Field_T *const field,const double *const Y,const Uint N)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  const double max_f = L_inf(patch->nn,field->v);
  const double tol = n[1]*n[1]*1e-14*(max_f > 1 ? max_f : 1.);
  double X,Z;
  double xc[3];/* Cartesian coord */
  double Xc[3];/* Curvilinear coord */
  double max = 0;
  Flag_T flg;
  Uint i,j;
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->Y_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (j = 0; j < N; ++j)/* -> choose different slices for testing */
  {
    interp_s->I = (Uint) floor(random_double(0,n[0],j));
    interp_s->K = (Uint) floor(random_double(0,n[2],1));
    
    X = node[i_j_k_to_ijk(n,interp_s->I,0,0)]->X[0];
    Z = node[i_j_k_to_ijk(n,0,0,interp_s->K)]->X[2];
    Xc[0] = X;
    Xc[2] = Z;
    for (i = 0; i < N; ++i)
    {
      interp_s->Y = Y[i];
      
      Xc[1] = Y[i];
      x_of_X(xc,Xc,patch);
      diff = poly5_f_point(xc[0],xc[1],xc[2])-execute_interpolation(interp_s);
      
      if (GRT(fabs(diff),tol))
      {
        flg = FOUND; max = (max<fabs(diff)?fabs(diff):max);
      }
    }
    
  }
  
  free_interpolation(interp_s);
  
  if (flg == FOUND)
  {printf(Pretty0"Max difference = %e\n",max); return TEST_UNSUCCESSFUL;}
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in Z direction and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_Z(Field_T *const field,const double *const Z,const Uint N)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  const double max_f = L_inf(patch->nn,field->v);
  const double tol = n[2]*n[2]*1e-14*(max_f > 1 ? max_f : 1.);
  double Y,X;
  double xc[3];/* Cartesian coord */
  double Xc[3];/* Curvilinear coord */
  double max = 0;
  Flag_T flg;
  Uint i,j;
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->Z_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (j = 0; j < N; ++j)/* -> choose different slices for testing */
  {
    interp_s->J = (Uint) floor(random_double(0,n[1],j));
    interp_s->I = (Uint) floor(random_double(0,n[0],1));
    
    Y = node[i_j_k_to_ijk(n,0,interp_s->J,0)]->X[1];
    X = node[i_j_k_to_ijk(n,interp_s->I,0,0)]->X[0];
    Xc[0] = X;
    Xc[1] = Y;
    for (i = 0; i < N; ++i)
    {
      interp_s->Z = Z[i];
      
      Xc[2] = Z[i];
      x_of_X(xc,Xc,patch);
      diff = poly5_f_point(xc[0],xc[1],xc[2])-execute_interpolation(interp_s);
      
      if (GRT(fabs(diff),tol))
      {
        flg = FOUND; max = (max<fabs(diff)?fabs(diff):max);
      }
    }
  }
  
  free_interpolation(interp_s);
  
  if (flg == FOUND)
  {printf(Pretty0"Max difference = %e\n",max); return TEST_UNSUCCESSFUL;}
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in X&Y directions and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_XY(Field_T *const field,const double *const X,const double *const Y,const Uint Nx,const Uint Ny)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  const double max_f = L_inf(patch->nn,field->v);
  const double tol = n[0]*n[0]*n[1]*n[1]*1e-14*(max_f > 1 ? max_f : 1.);
  double Z;
  double xc[3];/* Cartesian coord */
  double Xc[3];/* Curvilinear coord */
  double max = 0;
  Flag_T flg;
  Uint a,b,c;
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->XY_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (a = 0; a < Nx; ++a)/* -> choose different slices for testing */
  {
    interp_s->K = (Uint) floor(random_double(0,n[2],a));
    
    Z = node[i_j_k_to_ijk(n,0,0,interp_s->K)]->X[2];
    Xc[2] = Z;
    for (b = 0; b < Nx; ++b)
    {
      interp_s->X = X[b];
      Xc[0] = X[b];
      for (c = 0; c < Ny; ++c)
      {
        interp_s->Y = Y[c];
        Xc[1] = Y[c];
        
        x_of_X(xc,Xc,patch);
        diff = poly5_f_point(xc[0],xc[1],xc[2])-execute_interpolation(interp_s);
        if (GRT(fabs(diff),tol))
        {
          flg = FOUND; max = (max<fabs(diff)?fabs(diff):max);
        }
      }
    
    }
  }
  
  free_interpolation(interp_s);
  
  if (flg == FOUND)
  {printf(Pretty0"Max difference = %e\n",max); return TEST_UNSUCCESSFUL;}
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in X&Z directions and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_XZ(Field_T *const field,const double *const X,const double *const Z,const Uint Nx,const Uint Nz)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  const double max_f = L_inf(patch->nn,field->v);
  const double tol = n[0]*n[0]*n[2]*n[2]*1e-14*(max_f > 1 ? max_f : 1.);
  double Y;
  double xc[3];/* Cartesian coord */
  double Xc[3];/* Curvilinear coord */
  double max = 0;
  Flag_T flg;
  Uint a,b,c;
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->XZ_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (a = 0; a < Nx; ++a)/* -> choose different slices for testing */
  {
    interp_s->J = (Uint) floor(random_double(0,n[1],a));
    
    Y = node[i_j_k_to_ijk(n,0,interp_s->J,0)]->X[1];
    Xc[1] = Y;

    for (b = 0; b < Nx; ++b)
    {
      interp_s->X = X[b];
      Xc[0] = X[b];
      for (c = 0; c < Nz; ++c)
      {
        interp_s->Z = Z[c];
        Xc[2] = Z[c];
        x_of_X(xc,Xc,patch);
        diff = poly5_f_point(xc[0],xc[1],xc[2])-execute_interpolation(interp_s);
        
        if (GRT(fabs(diff),tol))
        {
          flg = FOUND; max = (max<fabs(diff)?fabs(diff):max);
        }
      }
    
    }
  }
  
  free_interpolation(interp_s);
  
  if (flg == FOUND)
  {printf(Pretty0"Max difference = %e\n",max); return TEST_UNSUCCESSFUL;}
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in Y&Z directions and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_YZ(Field_T *const field,const double *const Y,const double *const Z,const Uint Ny,const Uint Nz)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  const double max_f = L_inf(patch->nn,field->v);
  const double tol = n[2]*n[2]*n[1]*n[1]*1e-14*(max_f > 1 ? max_f : 1.);
  double X;
  double xc[3];/* Cartesian coord */
  double Xc[3];/* Curvilinear coord */
  double max = 0;
  Flag_T flg;
  Uint a,b,c;
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->YZ_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (a = 0; a < Ny; ++a)/* -> choose different slices for testing */
  {
    interp_s->I = (Uint) floor(random_double(0,n[0],a));
    
    X = node[i_j_k_to_ijk(n,interp_s->I,0,0)]->X[0];
    Xc[0] = X;
    
    for (b = 0; b < Ny; ++b)
    {
      interp_s->Y = Y[b];
      Xc[1] = Y[b];
      for (c = 0; c < Nz; ++c)
      {
        interp_s->Z = Z[c];
        Xc[2] = Z[c];
        x_of_X(xc,Xc,patch);
        
        diff = poly5_f_point(xc[0],xc[1],xc[2])-execute_interpolation(interp_s);
        if (GRT(fabs(diff),tol))
        {
          flg = FOUND; max = (max<fabs(diff)?fabs(diff):max);
        }
      }
    
    }
  }
  
  free_interpolation(interp_s);
  
  if (flg == FOUND)
  {printf(Pretty0"Max difference = %e\n",max); return TEST_UNSUCCESSFUL;}
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in X&Y&Z directions and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_XYZ(Field_T *const field,const double *const X,const double *const Y,const double *const Z,const Uint Nx,const Uint Ny,const Uint Nz)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const Uint *const n = patch->n;
  double diff = 0;
  double max_f = L_inf(patch->nn,field->v);
  const double tol = n[0]*n[0]*n[1]*n[1]*n[2]*n[2]*1e-14*(max_f > 1 ? max_f : 1.);
  double xc[3];/* Cartesian coord */
  double Xc[3];/* Curvilinear coord */
  double max = 0;
  Flag_T flg;
  Uint a,b,c;
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (a = 0; a < Nx; ++a)/* -> choose different slices for testing */
  {
    interp_s->X = X[a];
    Xc[0] = X[a];
    for (b = 0; b < Ny; ++b)
    {
      interp_s->Y = Y[b];
      Xc[1] = Y[b];
      for (c = 0; c < Nz; ++c)
      {
        interp_s->Z = Z[c];
        Xc[2] = Z[c];
        
        x_of_X(xc,Xc,patch);
        diff = poly5_f_point(xc[0],xc[1],xc[2])-execute_interpolation(interp_s);
        
        if (GRT(fabs(diff),tol))
        {
          flg = FOUND; max = (max<fabs(diff)?fabs(diff):max);
        }
      }
    
    }
  }
  
  free_interpolation(interp_s);
  
  if (flg == FOUND)
  {printf(Pretty0"Max difference = %e\n",max); return TEST_UNSUCCESSFUL;}
  
  return TEST_SUCCESSFUL;
}

/* testing:
// ========
//
// given a set of analytic functions, it expands them in these basis and
// take its various derivatives ans finally 
// compares the analytic results with numeric results.
// note: the results will be printed accordingly in 
// "Derivative_Tests" folder.
// note: only those patches that use basis will be compared.
// ->return value: EXIT_SUCCESS;
*/
int derivative_tests(Grid_T *const grid)
{
  sFunc_Patch2Pdouble_T **DataBase_func;
  
  const char *path_par;
  char *path;
  char der_s[MAXSTR];
  Uint fi;
  Flag_T flg;
  
  path_par = Pgets("top_directory");
  path = make_directory(path_par,"Derivative_Tests");

  
  init_func_Patch2Pdouble(&DataBase_func);
  /* below is all the available analytic functions with their derivatives. 
  // one can add more function, with the same style below.
  // note: derivative must be shown by underline '_' and direction of
  // derivative as it has been shown. it is "recommended" to use the same
  // notation for naming of functions.
  // note: functions are defined in Analytic folder in Maths.
  */
  
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f));/* note: although, this function is of order of 5, when you change coordinate system
                                                       // in the new coordinate, it is not necessarily of order of 5, or it might not even
                                                       // be polynomial at all. */
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_xyz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f));/* this function has branch point and pole at +i and -i
                                                         // in complex plane, thus, at r close to zero one needs more coeffs to resolve. So
                                                         // if in testing of derivative you notice error, be aware of this.
                                                         // this function is to be used mostly for large distances from the origin. */
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xyz));

  /*add_func_Patch2Pdouble(&DataBase_func,ArgM(c_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(x_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(y_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(z_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(c_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(c_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(c_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(x_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(x_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(x_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(x_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(y_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(y_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(y_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(y_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(z_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(z_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(z_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(z_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_xyz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_xyz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_xyz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_xyz));
  */
  
  FOR_ALL(fi,DataBase_func)
  {
    /* avoid counting some of funcs which has already been counted */
    if (DataBase_func[fi]->flg == 1) continue;
    
    sFunc_Patch2Pdouble_T *F[N_FUNC];
    
    double *anac[N_FUNC];/* analytic */
    double *numc[N_FUNC];/* numeric */
    enum FUNC_E e;
    struct Error_S error[1] = {0};
    Uint p;
    
    /* initializing, make them to point to 0 */
    for (e = FUNC; e < N_FUNC; ++e)
    {
      anac[e] = 0;
      numc[e] = 0;
    }
    
    /* read F from data base */
    flg = read_F(F,DataBase_func,(enum FUNC_E)fi);
    
    /* if the is no derivative for this function it continues */
    if (flg == NO) continue;
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      
      Field_T *df_num = add_field("$Numerica_derivative_TEST_FUNCTION$","(3dim)",patch,NO);
      
      /* compute anac if any */
      for (e = FUNC; e < N_FUNC; ++e)
        if (F[e])  anac[e] = F[e]->func(patch);
      
      for (e = FUNC_x; e < N_FUNC; ++e)
      {
        /* if anac is define and so not null, compare with numeric */
        if (anac[e])
        {
          printf("Testing Derivative: patch=%s, function:%15s\t",patch->name,F[e]->name);
          df_num->v = anac[FUNC];
          free_v2(df_num);
          enum2str(e,der_s);
          numc[e] = Partial_Derivative(df_num,der_s);
          flg = compare_derivative(F[e]->name,numc[e],anac[e],df_num,e,patch,path,error);
          free(anac[e]);
          free(numc[e]);
          
          if (flg == YES)
            printf("[+].(E_an,E_nu)=(%e,%e)\n",error->E_an,error->E_nu);
          else
            printf("[-].(E_an,E_nu)=(%e,%e)\n",error->E_an,error->E_nu);
        }
      }/* end of for (e = FUNC_x; e < N_FUNC; ++e) */
      free(anac[FUNC]);
      df_num->v = 0;
      remove_field(df_num);
  
    }/* end of FOR_ALL_PATCHES(pa,grid) */
    
  }/* end of FOR_ALL(fi,DataBase_func) */
  
  free_func_Patch2Pdouble(DataBase_func);
  free(path);
  
  return EXIT_SUCCESS;
}

/* get a function and based on its name, find all of its 
// derivative in data base.
// ->return value: if related derivative found YES, NO otherwise. 
*/
static Flag_T read_F(sFunc_Patch2Pdouble_T **const F,sFunc_Patch2Pdouble_T **const DataBase_func,const enum FUNC_E fn)
{
  const char *fname = DataBase_func[fn]->name;/* function name */
  char fname_derivative[MAXSTR];/* e.g. poly5_f_xy */
  enum FUNC_E e;
  Flag_T flg = NO;
  
  F[FUNC] = get_func_Patch2Pdouble(fname,DataBase_func);
  
  if (!F[FUNC])
    Errors("There is no function %s .\n",fname);
  else
    F[FUNC]->flg = 1;
    
  for (e = FUNC_x; e < N_FUNC; ++e)
  {
    sprintf(fname_derivative,"%s",fname);
    enum2strcat(e,fname_derivative);
    F[e] = get_func_Patch2Pdouble(fname_derivative,DataBase_func);
    
    /* if there is no such a derivative continue */
    if(!F[e]) continue;
    
    F[e]->flg = 1;
    flg = YES;
  }
  
  return flg;
}

/* appending a string based on the given enum e to the fname_derivatives */
static void enum2strcat(enum FUNC_E e,char *const fname_derivative)
{
  switch (e)
  {
    case FUNC_x:
      strcat(fname_derivative,"_x");
      break;
    case FUNC_y:
      strcat(fname_derivative,"_y");
      break;
    case FUNC_z:
      strcat(fname_derivative,"_z");
      break;
    case FUNC_xx:
      strcat(fname_derivative,"_xx");
      break;
    case FUNC_yy:
      strcat(fname_derivative,"_yy");
      break;
    case FUNC_zz:
      strcat(fname_derivative,"_zz");
      break;
    case FUNC_xy:
      strcat(fname_derivative,"_xy");
      break;
    case FUNC_xz:
      strcat(fname_derivative,"_xz");
      break;
    case FUNC_yz:
      strcat(fname_derivative,"_yz");
      break;
    case FUNC_xyz:
      strcat(fname_derivative,"_xyz");
      break;
    default:
      Error0("There is no such derivative defined.\n"
      "If you added more kind of derivative please add" 
        "to enum FUNC_E and consequently other locations.\n");
  }
}

/* filling a string str based on the given enum e  */
static void enum2str(enum FUNC_E e,char *const str)
{
  switch (e)
  {
    case FUNC_x:
      sprintf(str,"x");
      break;
    case FUNC_y:
      sprintf(str,"y");
      break;
    case FUNC_z:
      sprintf(str,"z");
      break;
    case FUNC_xx:
      sprintf(str,"x,x");
      break;
    case FUNC_yy:
      sprintf(str,"y,y");
      break;
    case FUNC_zz:
      sprintf(str,"z,z");
      break;
    case FUNC_xy:
      sprintf(str,"x,y");
      break;
    case FUNC_xz:
      sprintf(str,"x,z");
      break;
    case FUNC_yz:
      sprintf(str,"y,z");
      break;
    case FUNC_xyz:
      sprintf(str,"x,y,z");
      break;
    default:
      Error0("There is no such derivative defined.\n"
      "If you added more kind of derivative please add" 
        "to enum FUNC_E and consequently other locations.\n");
  }
}

/* comparing the values obtained from numeric and with analytic one */
static Flag_T compare_derivative(const char *const name,const double *const numc,const double *const anac,const Field_T *const func,const enum FUNC_E fn,const Patch_T *const patch,const char *const path,struct Error_S *const error)
{
  char prefix[MAXSTR];
  double E_anac, E_numc;/* numeric and analytic error */
  Flag_T flg;
  
  sprintf(prefix,"%s/%s.DiffByNode",path,name);
  error->E_nu = E_numc = pr_derivatives_DiffByNode(numc,anac,patch,prefix);
  error->E_an = E_anac = calculate_expected_precision_for_derivative(func,fn,patch);
  
  if (LSSEQL(E_numc/E_anac,10.0))
    flg = YES;
  else
    flg = NO;
    
  return flg;
}

/* calculating precision of derivative using the fact 
// that computer is using finite number of digits.
// ->return value: error in calculation = general idea is as follow:
// 1e-14*max(func)*max(Jacobian)^(order of derivative )*n*n^(2*order of derivative)*10  */
static double calculate_expected_precision_for_derivative(const Field_T *const func,const enum FUNC_E fn,const Patch_T *const patch)
{
  Uint o = order_of_derivative(fn);
  
  UNUSED(patch);
  
  return spectral_derivative_max_error(func,o);
}

/* ->retun value: order of derivative according to fn */
static Uint order_of_derivative(const enum FUNC_E fn)
{
  Uint o = 0;
  
  switch (fn)
  {
    case FUNC_x:
      o = 1;
      break;
    case FUNC_y:
      o = 1;
      break;
    case FUNC_z:
      o = 1;
      break;
    case FUNC_xx:
      o = 2;
      break;
    case FUNC_yy:
      o = 2;
      break;
    case FUNC_zz:
      o = 2;
      break;
    case FUNC_xy:
      o = 2;
      break;
    case FUNC_xz:
      o = 2;
      break;
    case FUNC_yz:
      o = 2;
      break;
    case FUNC_xyz:
      o = 3;
      break;
    default:
      Error0("There is no such derivative defined.\n"
      "If you added more kind of derivative please add" 
        "to enum FUNC_E and consequently other locations.\n");
  }
  
  return o;
}

/* making N random number in the within min and max.
// ->return value: random number(s)
*/
static double *make_random_number(const Uint N,const double min,const double max)
{
  double *x = calloc(N,sizeof(*x));
  Uint i;
  
  for (i = 0; i < N; ++i)
  {
    x[i] = random_double(min,max,i);
  }
  
  return x;
}

/* free function structure form grid to pointer to double */
static void free_func_Patch2Pdouble(sFunc_Patch2Pdouble_T **func)
{
  Uint i;
  
  for (i = 0; func[i] != 0; ++i)
  {
    free(func[i]->name);
    free(func[i]);
  }
  
  free(func);
}

