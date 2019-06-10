/*
// Alireza Rashti
// June 2019
*/

#include "TOV.h"

/* solving TOV equations. it approximates the properties of a NS.
// it only deals with inside of the star, since the outside is analytic.
// ->return value: solution of TOV equations for inside of a star with give EoS */
TOV_T *TOV_solution(TOV_T *const TOV)
{
  const double Fac = 0.25;
  double h_cent_prev = 1;
  double a = 0,b = 0;
  double m = 0;/* ADM mass of NS for various central enthalpy */
  Flag_T increase  = YES;
  Flag_T bisection = NO;
  EoS_T *eos = 0;
  unsigned i;
 
  pr_line_custom('=');
  printf("Solving TOV equations for %s ...\n",TOV->description);
  
  /* some initialization and preparation */
  TOV->h_cent = 1.5;
  TOV->m = alloc_double(TOV->N);
  TOV->r = alloc_double(TOV->N);
  TOV->h = alloc_double(TOV->N);
  TOV->p = alloc_double(TOV->N);
  
  /* find enthalpy at the center of NS such that 
  // the ADM mass reaches the desired value.
  // the spirit of the approach is bisection method in root finders */
  while (!EQL(m,TOV->ADM_m))
  {
    solve_ODE_enthalpy_approach(TOV);
    m = TOV->m[TOV->N-1];
    
    /* as long as the mass is less than desired mass, increase enthalpy
    // to find the range of enthalpy */
    if (increase == YES)
    {
      if (LSS(m,TOV->ADM_m))
      {
        h_cent_prev = TOV->h_cent;
        /* increase the enthalpy: */
        TOV->h_cent += TOV->h_cent*Fac;
      }
      else
      {
        increase  = NO;
        bisection = YES;/* active bisection */
        a = h_cent_prev;
        b = TOV->h_cent;
      }
    }
    /* now we know the range of enthalpy and we use bisection method
    // to pin down the correct value of the central enthalpy */
    if (bisection == YES)
    {
      if (GRT(m,TOV->ADM_m))
      {
        b = TOV->h_cent;
      }
      else
      {
        a = TOV->h_cent;
      }
      TOV->h_cent = a+0.5*(b-a);
    }
  }
  TOV->ADM_m = TOV->m[TOV->N-1];
  TOV->bar_m = calculate_baryonic_mass(TOV);
  calculate_phi(TOV);
  calculate_ADM_and_Komar_mass(TOV);
  
  /* having known every thing, now populate pressure */
  eos = initialize_EoS();
  for (i = 0; i < TOV->N; ++i)
  {
    eos->h = TOV->h[i];
    TOV->p[i] = eos->pressure(eos);
  }
  
  /* print some informations about TOV */
  printf("TOV properties:\n");
  printf("--> NS radius (Schwarzschild Coords.) = %g\n",TOV->r[TOV->N-1]);
  printf("--> NS radius (Isotropic Coords.)     = %g\n",-10000.);
  printf("--> ADM mass                          = %g\n",TOV->ADM_m);
  printf("--> baryonic mass                     = %g\n",TOV->bar_m);
  printf("--> central pressure                  = %g\n",TOV->p[0]);
  eos->h = TOV->h_cent;
  printf("--> central energy density            = %g\n",eos->energy_density(eos));
  
  free_EoS(eos);
  
  printf("Solving TOV equations for %s ==> Done.\n",TOV->description);
  pr_clock();
  pr_line_custom('=');
  
  return TOV;
}

/* calculate phi in g_00 = - exp[2phi] for the metric of space time */
static void calculate_phi(TOV_T *const TOV)
{
  const double R = TOV->r[TOV->N-1],
               M = TOV->m[TOV->N-1],
               *const h = TOV->h;
  unsigned i;
  TOV->phi = alloc_double(TOV->N);
  
  for (i = 0; i < TOV->N; ++i)
    TOV->phi[i] = -log(h[i])+log(sqrt(1-2*M/R));
  
}

/* calculate Komar mass and ADM mass and check virial theorem */
static void calculate_ADM_and_Komar_mass(TOV_T *const TOV)
{
  const double tol = 1e-4;/* this is due to the integration */
  Integration_T *I = init_integration();
  double *f = 0;
  double Komar_mass, ADM_mass;
  
  I->type = "Composite Simpson's Rule 1D";
  I->Composite_Simpson_1D->b = 1;/* since the variable is enthalpy */
  I->Composite_Simpson_1D->a = TOV->h_cent;
  I->Composite_Simpson_1D->n = TOV->N;
  f = Komar_mass_integrand(TOV);
  I->Composite_Simpson_1D->f = f;
  plan_integration(I);
  Komar_mass = execute_integration(I);
  free(f);
  
  f = ADM_mass_integrand(TOV);
  I->Composite_Simpson_1D->f = f;
  plan_integration(I);
  ADM_mass = execute_integration(I);
  free(f);
  free_integration(I);
  
  /* some test control */
  if (GRT(fabs(Komar_mass-ADM_mass),tol))/* virial theorem */
  {
    fprintf(stderr,"Komar mass = %g, ADM mass = %g\n",Komar_mass,ADM_mass);
    abortEr("Komar mass and ADM mass must be equal!\n");
  }
}

/* populate the integrand needed in calculation of Komar mass.
// ->return value: integrand for Komar mass calculation. */
static double *Komar_mass_integrand(const TOV_T *const TOV)
{
  double *f = alloc_double(TOV->N);
  const double *const r = TOV->r,
               *const m = TOV->m,
               *const h = TOV->h,
               *const phi = TOV->phi;
  EoS_T *eos = initialize_EoS();
  unsigned i;
  
  f[0] = 0;
  for (i = 1; i < TOV->N; ++i)
  {
    eos->h = h[i];
    f[i] = 4*M_PI*(eos->energy_density(eos)+3*eos->pressure(eos))*
            exp(phi[i])/sqrt(1-2*m[i]/r[i])*SQR(r[i])*dr_dh(h[i],r[i],m[i]);
  }
  
  free_EoS(eos);
  
  return f;
}

/* populate the integrand needed in calculation of ADM mass.
// ->return value: integrand for ADM mass calculation. */
static double *ADM_mass_integrand(const TOV_T *const TOV)
{
  double *f = alloc_double(TOV->N);
  const double *const r = TOV->r,
               *const m = TOV->m,
               *const h = TOV->h;
  EoS_T *eos = initialize_EoS();
  unsigned i;
  
  f[0] = 0;
  for (i = 1; i < TOV->N; ++i)
  {
    eos->h = h[i];
    f[i] = 4*M_PI*eos->energy_density(eos)*SQR(r[i])*dr_dh(h[i],r[i],m[i]);
  }
  
  free_EoS(eos);
  
  return f;
}

/* taking the integral of 4\pi \int ^{R}_{0}\rho \left( 1-\frac {2m}{r}\right) ^{-\frac {1}{2}}r^{2}dr.
// ->return value: baryonic rest mass */
static double calculate_baryonic_mass(const TOV_T *const TOV)
{
  Integration_T *I = init_integration();
  double *f = baryonic_mass_integrand(TOV);
  double integral;/* resultant */
  
  I->type = "Composite Simpson's Rule 1D";
  I->Composite_Simpson_1D->b = 1;/* since the variable is enthalpy */
  I->Composite_Simpson_1D->a = TOV->h_cent;
  I->Composite_Simpson_1D->n = TOV->N;
  I->Composite_Simpson_1D->f = f;
  plan_integration(I);
  integral = execute_integration(I);
  
  free(f);
  free_integration(I);
  
  return integral;
}

/* populate the integrand needed in calculation of baryonic rest mass.
// ->return value: integrand for baryonic rest mass calculation. */
static double *baryonic_mass_integrand(const TOV_T *const TOV)
{
  double *f = alloc_double(TOV->N);
  double rho;
  const double *const r = TOV->r,
               *const m = TOV->m,
               *const h = TOV->h;
  EoS_T *eos = initialize_EoS();
  double s;
  unsigned i;
  
  f[0] = 0;
  for (i = 1; i < TOV->N; ++i)
  {
    eos->h = h[i];
    rho = (eos->energy_density(eos)+eos->pressure(eos))/h[i];
    s = 1-2*m[i]/r[i];
    f[i] = 4*M_PI*rho/sqrt(s)*SQR(r[i])*dr_dh(h[i],r[i],m[i]);
  }
  
  free_EoS(eos);
  
  return f;
}

/* given central enthalpy it solves ODE of TOV equations
// using enthalpy as independent variable. this approach is useful
// since it avoids singularity that TOV equations have when r is used
// as independent variable.
// in this approach we write equation for m and r in terms of h and for
// phi we solved it analytically, m and r equations are:
// 
// r:
// \frac {dr}{dh}=-\frac {r\left( r-2m\right) }{\left( m+4\pi r^{3}p\right) h}\\ .
// m:
// \frac {dm}{dh}=4\pi r^{2}e\frac {dr}{dh}\\ .
//
// To solve ODE, we use Runge-Kutta method of 4th order. */
static void solve_ODE_enthalpy_approach(TOV_T *const TOV)
{
  double *const r = TOV->r;
  double *const m = TOV->m;
  double *const h = TOV->h;
  const double b = 1;/* h at the NS surface */
  const double a = TOV->h_cent;/* h at the center of NS */
  const double s = (b-a)/(TOV->N-1);/* step size of h. note:
                               // it's negative since we start 
                               // from the center which has greater value
                               // of h compare to the NS's surface. */
  double t;/* independent variable h, we conventionally called it t */
  unsigned i;
  enum {R = 0,M = 1};
  
  /* initialization */
  h[0] = a;/* h at center */
  h[1] = a+s;
  r[0] = 0;/* r at center is 0 */
  r[1] = r_approx(h[1],a);
  m[0] = 0;/* m at center is 0 */
  m[1] = m_approx(h[1],a);
  t = h[1];
  
  /* for all points */
  for (i = 2; i < TOV->N; ++i)
  {
    double k1[2],k2[2],k3[2],k4[2];/* variables for Runge-Kutta method of 4th order */
    
    /* dr/dh and dm/dh equations: */
    k1[R] = s*dr_dh(t,r[i-1],m[i-1]);
    k1[M] = s*dm_dh(t,r[i-1],m[i-1]);
    
    k2[R] = s*dr_dh(t+s/2,r[i-1]+k1[R]/2,m[i-1]+k1[M]/2);
    k2[M] = s*dm_dh(t+s/2,r[i-1]+k1[R]/2,m[i-1]+k1[M]/2);
    
    k3[R] = s*dr_dh(t+s/2,r[i-1]+k2[R]/2,m[i-1]+k2[M]/2);
    k3[M] = s*dm_dh(t+s/2,r[i-1]+k2[R]/2,m[i-1]+k2[M]/2);
    
    k4[R] = s*dr_dh(t+s,r[i-1]+k3[R],m[i-1]+k3[M]);
    k4[M] = s*dm_dh(t+s,r[i-1]+k3[R],m[i-1]+k3[M]);
    
    /* updating the values */
    r[i] = r[i-1]+(k1[R]+2*k2[R]+2*k3[R]+k4[R])/6;
    m[i] = m[i-1]+(k1[M]+2*k2[M]+2*k3[M]+k4[M])/6;
    t = a+i*s;
    h[i] = t;
  }
  
}

/* ->return value: approximate r near center of star */
static double r_approx(const double h,const double h_c/* central enthalpy */)
{
  EoS_T *eos = initialize_EoS();
  double r = 0;
  double e,p,de_dh;
  
  eos->h = h_c;
  e = eos->energy_density(eos);
  p = eos->pressure(eos);
  de_dh = eos->de_dh(eos);
  free_EoS(eos);
  
  r = sqrt(3*(h_c-h)/(2*M_PI*(e+3*p)))*
      (1-0.25*(e-3*p-3*de_dh/5)*(h_c-h)/(e+3*p));
  
  return r;
}

/* ->return value: approximate m near center of star */
static double m_approx(const double h,const double h_c/* central enthalpy */)
{
  EoS_T *eos = initialize_EoS();
  double m = 0;
  double e,de_dh;
  
  eos->h = h_c;
  e = eos->energy_density(eos);
  de_dh = eos->de_dh(eos);
  free_EoS(eos);
  
  m = 4*M_PI/3*e*pow(r_approx(h,h_c),3)*
      (1-3./5.*de_dh*(h_c-h)/e);
  
  return m;
}

/* rbar equation:
// \frac {d\overline {r}}{dh}=\frac {\overline {r}}{r\sqrt {1-2\frac {m}{r}}}\frac {dr}{dh} */
//static double drbar_dh(const double h,const double rbar,const double r, const double m)
//{
  //return rbar/sqrt(1-2*m/r)/r*dr_dh(h,r,m);
//}

/* r equation:
// \frac {dr}{dh}=-\frac {r\left( r-2m\right) }{\left( m+4\pi r^{3}p\right) h}\\ . */
static double dr_dh(const double h,const double r, const double m)
{
  EoS_T *eos = initialize_EoS();
  double f;
  double p;/* pressure */
  const double r3 = r*r*r;
  
  eos->h = h;
  p = eos->pressure(eos); 
  f = - r*(r-2*m)/(m+4*M_PI*r3*p)/h;
  free_EoS(eos);

  return f;
}

/* m equation:
// \frac {dm}{dh}=4\pi r^{2}e\frac {dr}{dh}\\ . */
static double dm_dh(const double h,const double r, const double m)
{
  EoS_T *eos = initialize_EoS();
  double f;
  double e;/* energy density */
  const double r2 = r*r;
  
  eos->h = h;
  e = eos->energy_density(eos); 
  f = 4*M_PI*r2*e*dr_dh(h,r,m);;
  
  free_EoS(eos);

  return f;
}

/* ->return value: a pristine TOV struct */
TOV_T *TOV_init(void)
{
  TOV_T *tov = calloc(1,sizeof(*tov));
  
  return tov;
}

/* free thoroughly the given struct */
void TOV_free(TOV_T *TOV)
{
  free(TOV->m);
  free(TOV->r);
  free(TOV->p);
  free(TOV->h);
  free(TOV->phi);
  free(TOV->rbar);
  free(TOV->psi);
  free(TOV);
}
