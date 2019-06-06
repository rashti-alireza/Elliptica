r/*
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
  const double Fac   = 10;/* increase factor*/
  double h_cent_prev = 1;
  double a = 0,b = 0;
  Flag_T increase  = YES;
  Flag_T bisection = NO;
  
  /* some initialization and preparation */
  TOV->ref_fac = 5;/* some odd number for refinement */
  TOV->ref_N = TOV->ref_fac*TOV->N;
  TOV->calculated_baryonic_m = 0;
  TOV->h_cent = 10;
  TOV->m = alloc_double(TOV->ref_N);
  TOV->r = alloc_double(TOV->ref_N);
  TOV->h = alloc_double(TOV->ref_N);
  
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
  
  /* now transform this solution in conformal decompostion form */
  conformal_decomposition_transformation(TOV);
  
  return TOV;
}

/* transform the TOV solution in regular format in conformal decompostion format.
// note: Ln(rbar/r) + C = integral f */
static void conformal_decomposition_transformation(TOV_T *const TOV)
{
  Integration_T *I = init_integration();
  double *f = 0,C = 0;
  const double R = TOV->r[TOV->ref_N],
               M = TOV->m[TOV->ref_N];
  const double *const h = TOV->h,
               *const r = TOV->r;
  double *rbar = 0, *psi = 0;
  unsigned i;
  
  I->type = "Composite Simpson's Rule 1D";
  I->Composite_Simpson_1D->b = 1;/* since the variable is enthalpy */
  I->Composite_Simpson_1D->a = TOV->h_cent;
  I->Composite_Simpson_1D->n = TOV->ref_N;
  f = con_dec_integrand(TOV);
  I->Composite_Simpson_1D->f = f;
  plan_integration(I);
  C  = execute_integration(I);
  C -= log((R-M+sqrt(SQR(R)-2*M*R))/(2*R));
  
  
  /* fill rbar */
  TOV->Conformal_Decomposition->r    = alloc_double(TOV->N);
  rbar = TOV->Conformal_Decomposition->r;
  rbar[0] = 0;
  for (i = 1; i < TOV->N; ++i)
  {
    I->Composite_Simpson_1D->b = h[i*TOV->ref_fac];/* since the variable is enthalpy */
    I->Composite_Simpson_1D->n = i*TOV->ref_fac;
    rbar[i] = r[i*TOV->ref_fac]*exp(execute_integration(I)-C);
  }
  
  TOV->Conformal_Decomposition->psi = alloc_double(TOV->N);
  psi = TOV->Conformal_Decomposition->psi;
  psi[0] = sqrt(exp(C));
  for (i = 1; i < TOV->N; ++i)
    psi[i] = sqrt(r[i*TOV->ref_fac]/rbar[i]);
  
  free(f);
  free_integration(I);
  
}

/* ->return value:the integrand needed in calculation r bar in conformal decomposition. */
static double *con_dec_integrand(const TOV_T *const TOV)
{
  double *f = alloc_double(TOV->ref_N);
  const double *const r = TOV->r,
               *const m = TOV->m,
               *const h = TOV->h;
  unsigned i;
  
  f[0] = 0;
  for (i = 1; i < TOV->ref_N; ++i)
  {
    f[i] = (1-sqrt(1-2*m[i]/r[i]))/(r[i]*sqrt(1-2*m[i]/r[i]))*dr_dh(h[i],r[i],m[i]);
  }
  
  return f;
}


/* calculate phi in g_00 = - exp[2phi] for the metric of space time */
static void calculate_phi(TOV_T *const TOV)
{
  const double R = TOV->r[TOV->ref_N],
               M = TOV->m[TOV->ref_N],
               *const h = TOV->h;
  unsigned i;
  TOV->phi = alloc_double(TOV->ref_N);
  
  for (i = 0; i < TOV->ref_N; ++i)
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
  I->Composite_Simpson_1D->n = TOV->ref_N;
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
  double *f = alloc_double(TOV->ref_N);
  const double *const r = TOV->r,
               *const m = TOV->m,
               *const h = TOV->h,
               *const phi = TOV->phi;
  EoS_T *eos = initialize_EoS();
  unsigned i;
  
  for (i = 0; i < TOV->ref_N; ++i)
  {
    eos->h = h[i];
    f[i] = 4*M_PI*(eos->energy_density(eos)+3*eos->pressure(eos))*
            exp(phi[i])/sqrt(1-2*m[i]/r[i])*SQR(r[i])*dr_dh(h[i],r[i],m[i]);
  }
  
  free_EoS(&eos);
  
  return f;
}

/* populate the integrand needed in calculation of ADM mass.
// ->return value: integrand for ADM mass calculation. */
static double *ADM_mass_integrand(const TOV_T *const TOV)
{
  double *f = alloc_double(TOV->ref_N);
  const double *const r = TOV->r,
               *const m = TOV->m,
               *const h = TOV->h;
  EoS_T *eos = initialize_EoS();
  unsigned i;
  
  for (i = 0; i < TOV->ref_N; ++i)
  {
    eos->h = h[i];
    f[i] = 4*M_PI*eos->energy_density(eos)*SQR(r[i])*dr_dh(h[i],r[i],m[i]);
  }
  
  free_EoS(&eos);
  
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
  I->Composite_Simpson_1D->n = TOV->ref_N;
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
  double *f = alloc_double(TOV->ref_N);
  double rho;
  const double *const r = TOV->r,
               *const m = TOV->m,
               *const h = TOV->h;
  EoS_T *eos = initialize_EoS();
  unsigned i;
  
  for (i = 0; i < TOV->ref_N; ++i)
  {
    eos->h = h[i];
    rho = (eos->energy_density(eos)+eos->pressure(eos))/h[i];
    f[i] = rho/sqrt(1-2*m[i]/r[i])*SQR(r[i])*dr_dh(h[i],r[i],m[i]);
  }
  
  free_EoS(&eos);
  
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
  const double s = (b-a)/(TOV->ref_N-1);/* step size of h. note:
                               // it's negative since we start 
                               // from the center which has greater value
                               // of h compare to the NS's surface. */
  double t;/* independent variable h, we conventionally called it t */
  unsigned i;
  enum {R = 0,M = 1};
  
  /* initialization */
  t = a;
  r[0] = 0;/* r at center is 0 */
  m[0] = 0;/* m at center is 0 */
  h[0] = a;/* h at center */
  
  /* for all points */
  for (i = 1; i < TOV->ref_N; ++i)
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
    
    r[i] = r[i-1]+(k1[R]+2*k2[R]+2*k3[R]+k4[R])/6;
    m[i] = m[i-1]+(k1[M]+2*k2[M]+2*k3[M]+k4[M])/6;
    h[i] = t;
    t = a+i*s;
  }
  
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
  
  free_EoS(&eos);

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
  
  free_EoS(&eos);

  return f;
}
