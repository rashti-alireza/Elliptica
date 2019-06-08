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
  const double tol = 1E-5;/* tolerance in mass calculations */
  const double Fac = 10;
  double h_cent_prev = 1;
  double a = 0,b = 0;
  Flag_T increase  = YES;
  Flag_T bisection = NO;
  EoS_T *eos = 0;
  unsigned i;
  
  /* some initialization and preparation */
  TOV->calculated_baryonic_m = 0;
  TOV->h_cent = 10;
  TOV->m = alloc_double(TOV->N);
  TOV->r = alloc_double(TOV->N);
  TOV->h = alloc_double(TOV->N);
  TOV->p = alloc_double(TOV->N);
  TOV->rbar = alloc_double(TOV->N);
  
  /* find enthalpy at the center of NS such that 
  // the baryonic mass reaches the desired value.
  // the spirit of the approach is bisection method in root finders */
  while (GRT(fabs(TOV->desired_baryonic_m-TOV->calculated_baryonic_m),tol))
  {
    solve_ODE_enthalpy_approach(TOV);
    
    TOV->calculated_baryonic_m = calculate_baryonic_mass(TOV);
    
    /* as long as the mass is less than desired mass, increase enthalpy
    // to find the range of enthalpy */
    if (increase == YES)
    {
      if (LSS(TOV->calculated_baryonic_m,TOV->desired_baryonic_m))
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
      if (GRT(TOV->calculated_baryonic_m,TOV->desired_baryonic_m))
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
  
  calculate_phi(TOV);
  check_virial_relation(TOV);
  
  /* having known every thing, now populate pressure */
  eos = initialize_EoS();
  for (i = 0; i < TOV->N; ++i)
  {
    eos->h = TOV->h[i];
    TOV->p[i] = eos->pressure(eos);
  }
  
  free_EoS(eos);
  
  return TOV;
}

/* calculate phi in g_00 = - exp[2phi] for the metric of space time */
static void calculate_phi(TOV_T *const TOV)
{
  const double R = TOV->r[TOV->N],
               M = TOV->m[TOV->N],
               *const h = TOV->h;
  unsigned i;
  TOV->phi = alloc_double(TOV->N);
  
  for (i = 0; i < TOV->N; ++i)
    TOV->phi[i] = -log(h[i])+log(sqrt(1-2*M/R));
  
}

/* check if the Komar mass and ADM mass are equal */
static void check_virial_relation(const TOV_T *const TOV)
{
  Integration_T *I = init_integration();
  double *f;
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
  
  printf("Komar mass = %g, ADM mass = %g\n",Komar_mass,ADM_mass);
  if (!EQL(Komar_mass,ADM_mass))
  {
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
  
  for (i = 0; i < TOV->N; ++i)
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
  
  for (i = 0; i < TOV->N; ++i)
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
  unsigned i;
  
  for (i = 0; i < TOV->N; ++i)
  {
    eos->h = h[i];
    rho = (eos->energy_density(eos)+eos->pressure(eos))/h[i];
    f[i] = rho/sqrt(1-2*m[i]/r[i])*SQR(r[i])*dr_dh(h[i],r[i],m[i]);
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
  double *const rbar = TOV->rbar;
  const double b = 1;/* h at the NS surface */
  const double a = TOV->h_cent;/* h at the center of NS */
  const double s = (b-a)/(TOV->N-1);/* step size of h. note:
                               // it's negative since we start 
                               // from the center which has greater value
                               // of h compare to the NS's surface. */
  double t;/* independent variable h, we conventionally called it t */
  unsigned i;
  enum {R = 0,M = 1,Rbar = 2};
  
  /* initialization */
  t = a;
  r[0] = 0;/* r at center is 0 */
  m[0] = 0;/* m at center is 0 */
  h[0] = a;/* h at center */
  rbar[0] = 0;/* rbar at center is 0 */
  
  /* for all points */
  for (i = 1; i < TOV->N; ++i)
  {
    double k1[3],k2[3],k3[3],k4[3];/* variables for Runge-Kutta method of 4th order */
    
    /* dr/dh and dm/dh equations: */
    k1[R] = s*dr_dh(t,r[i-1],m[i-1]);
    k1[M] = s*dm_dh(t,r[i-1],m[i-1]);
    
    k2[R] = s*dr_dh(t+s/2,r[i-1]+k1[R]/2,m[i-1]+k1[M]/2);
    k2[M] = s*dm_dh(t+s/2,r[i-1]+k1[R]/2,m[i-1]+k1[M]/2);
    
    k3[R] = s*dr_dh(t+s/2,r[i-1]+k2[R]/2,m[i-1]+k2[M]/2);
    k3[M] = s*dm_dh(t+s/2,r[i-1]+k2[R]/2,m[i-1]+k2[M]/2);
    
    k4[R] = s*dr_dh(t+s,r[i-1]+k3[R],m[i-1]+k3[M]);
    k4[M] = s*dm_dh(t+s,r[i-1]+k3[R],m[i-1]+k3[M]);
    
    /* isotropic coordinate transformation for drbar/dh */
    k1[Rbar] = s*drbar_dh(t,rbar[i-1],r[i-1],m[i-1]);
    k2[Rbar] = s*drbar_dh(t+s/2,rbar[i-1]+k1[Rbar]/2,r[i-1]+k1[R]/2,m[i-1]+k1[M]/2);
    k3[Rbar] = s*drbar_dh(t+s/2,rbar[i-1]+k2[Rbar]/2,r[i-1]+k2[R]/2,m[i-1]+k2[M]/2);
    k4[Rbar] = s*drbar_dh(t+s,rbar[i-1]+k3[Rbar],r[i-1]+k3[R],m[i-1]+k3[M]);

    /* updating the values */
    r[i] = r[i-1]+(k1[R]+2*k2[R]+2*k3[R]+k4[R])/6;
    m[i] = m[i-1]+(k1[M]+2*k2[M]+2*k3[M]+k4[M])/6;
    h[i] = t;
    rbar[i] = rbar[i-1]+(k1[Rbar]+2*k2[Rbar]+2*k3[Rbar]+k4[Rbar])/6;
    t = a+i*s;
  }
  
}


/* rbar equation:
// \frac {d\overline {r}}{dh}=\frac {\overline {r}}{r\sqrt {1-2\frac {m}{r}}}\frac {dr}{dh} */
static double drbar_dh(const double h,const double rbar,const double r, const double m)
{
  return rbar/sqrt(1-2*m/r)/r*dr_dh(h,r,m);
}

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
