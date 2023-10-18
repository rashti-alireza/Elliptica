/*
// Alireza Rashti
// July 2018
*/

#include "spectral_tests.h"
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
  
  FOR_ALL_PATCHES(p,grid)
  {
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
  }
  
  if (DO)
  {
      printf("Interpolation test:            Neville Iterative Method =>");
      status = interpolation_tests_Neville_1d();
      check_test_result(status);
  }
  if (DO)
  {
      printf("Interpolation test:            Natural Cubic Spline Method =>");
      status = interpolation_tests_N_cubic_spline_1d();
      check_test_result(status);
  }
  if (DO)
  {
      printf("Interpolation test:            Hermite 1D Method =>");
      status = interpolation_tests_Hermite_1d();
      check_test_result(status);
  }
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* test Natural Cubic Spline method for 1-d arrays.
// ->return value: result of test. */
static int interpolation_tests_N_cubic_spline_1d(void)
{
  Interpolation_T *interp_s = init_interpolation();
  const Uint N = (Uint)Pgeti("n_a");
  double *f = alloc_double(N);
  double *x = alloc_double(N);
  const double a = -M_PI, b = 3./4.*M_PI;/* an arbitrary interval  */
  double *hs = make_random_number(N,a,b);
  double s = (b-a)/(N-1);
  double t,interp;
  Flag_T flg = NONE;
  Uint i;
  
  for (i = 0; i < N; ++i)
  {
    t = x[i] = a+i*s;
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
    
    if (GRT(fabs(diff),s))
    {
      fprintf(stderr,"diff = %g\n",diff);
      flg = FOUND;
      break;
    }
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
  free(x);
  free(hs);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
    
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
  const double a = -M_PI, b = 3./4.*M_PI;/* an arbitrary interval  */
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
  
  flg = NONE;
  for (j = 0; j < N; ++j)/* -> choose different slices for testing */
  {
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

int derivative_tests(Grid_T *const grid)
{
  if (DO)
  {
    derivative_tests_spectral_method(grid);
  }
  if (DO)
  {
    derivative_tests_finite_diff_method(grid);
  }
  
  return EXIT_SUCCESS;
}

// testing finite difference derivatives
static void derivative_tests_finite_diff_method(Grid_T *const grid)
{
  const Uint N = (Uint)Pgeti("n_a");
  const Uint fd_acc = 3; // note: we are using 3rd order polynomial
  double *f   = alloc_double(N);
  double *fp  = alloc_double(N);
  double *fpp = alloc_double(N);
  double *x   = alloc_double(N);
  double (*f_t)(const double x) = f_poly_3deg1;
  double (*fp_t)(const double x) = df_poly_3deg1;
  double (*fpp_t)(const double x) = ddf_poly_3deg1;
  const double a = -M_PI, b = 3./4.*M_PI;/* an arbitrary interval  */
  double s = (b-a)/(N-1);
  double t;
  Uint i;
  
  for (i = 0; i < N; ++i)
  {
    x[i]   = 
    t      = a+i*s;
    f[i]   = f_t(t);
    fp[i]  = fp_t(t);
    fpp[i] = fpp_t(t);
  }
  
  for (i = 0; i < N; ++i)
  {
    t = a+i*s;
    double df_an  = finite_difference_Fornberg(f,x,t,N,1,fd_acc);
    double ddf_an = finite_difference_Fornberg(f,x,t,N,2,fd_acc);
    
    printf("|df_an  - df_nu|  = %0.6e\n", fabs(fp_t(t)  - df_an));
    printf("|ddf_an - ddf_nu| = %0.6e\n", fabs(fpp_t(t) - ddf_an));
  }
  
  Free(x);
  Free(f);
  Free(fp);
  Free(fpp);
  UNUSED(grid);
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
*/
static void derivative_tests_spectral_method(Grid_T *const grid)
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

/* test Hermite 1-d method
// ->return value: result of test. */
static int interpolation_tests_Hermite_1d(void)
{
  Interpolation_T *interp_s = init_interpolation();
  const Uint N = (Uint)Pgeti("n_a");
  const Uint Nh = N/2; assert(N>2);// num tests
  const Uint n_pnt  = 3;
  const Uint fd_acc = 3;
  double *f = alloc_double(N);
  double *fp = alloc_double(N);
  double *x = alloc_double(N);
  double (*f_t)(const double x) = f_poly_3deg1;
  double (*fp_t)(const double x) = df_poly_3deg1;
  const double a = -M_PI, b = 3./4.*M_PI;/* an arbitrary interval  */
  double *hs = make_random_number(Nh,a,b/5.); // b/5 so not too close to the end
  double s = (b-a)/(N-1);
  double t,interp;
  Flag_T flg = NONE;
  Uint i;
  
  for (i = 0; i < N; ++i)
  {
    t = x[i] = a+i*s;
    f[i] = f_t(t);
  }
  
  // with the internal derivatives  
  interp_s->method = "Hermite1D";
  interp_s->Hermite_1d->f = f;
  interp_s->Hermite_1d->x = x;
  interp_s->Hermite_1d->N = N;
  interp_s->Hermite_1d->fd_accuracy_order = fd_acc;
  interp_s->Hermite_1d->num_points = n_pnt;
  plan_interpolation(interp_s);
  
  for (i = 0; i < Nh; ++i)
  {
    double diff;
    t = hs[i];
    interp_s->Hermite_1d->h = t;
    interp = execute_interpolation(interp_s);
    diff = interp-f_t(t);
    
    if (GRT(fabs(diff),s))
    {
      fprintf(stderr,"diff = %g\n",diff);
      flg = FOUND;
      break;
    }
  }
  free_interpolation(interp_s);
  
  /* let's test the reveres order for x's */
  s = -(b-a)/(N-1);
  interp_s = init_interpolation();
  for (i = 0; i < N; ++i)
  {
    t = x[i] = b+i*s;
    f[i] = f_t(t);
    fp[i] = fp_t(t);
  }
  
  // with given derivative
  interp_s->method = "Hermite1D";
  interp_s->Hermite_1d->f = f;
  interp_s->Hermite_1d->fp = fp;
  interp_s->Hermite_1d->x = x;
  interp_s->Hermite_1d->N = N;
  interp_s->Hermite_1d->num_points = n_pnt;
  plan_interpolation(interp_s);
  
  for (i = 0; i < Nh; ++i)
  {
    double diff;
    t = hs[i];
    interp_s->Hermite_1d->h = t;
    interp = execute_interpolation(interp_s);
    diff = interp-f_t(t);
    
    if (GRT(fabs(diff),-s))
    {
      fprintf(stderr,"diff = %g\n",diff);
      flg = FOUND;
      break;
    }
  }
  
  free_interpolation(interp_s);
  free(f);
  free(fp);
  free(x);
  free(hs);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
    
  return TEST_SUCCESSFUL;
}

// an arbitray 3rd degree polynomial
static double f_poly_3deg1(const double x)
{
  return 1.33*pow(x,3.) + 0.33*pow(x,2.) + 5.2;
}

// the 1st derivative of an arbitray 3rd degree polynomial
static double df_poly_3deg1(const double x)
{
  return 3.*1.33*pow(x,2.) + 2.*0.33*x;
}

// the 2nd derivative of an arbitray 3rd degree polynomial
static double ddf_poly_3deg1(const double x)
{
  return 2*3.*1.33*x + 2.*0.33;
}

