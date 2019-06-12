/*
// Alireza Rashti
// June 2019
*/

#include "TOV_solution.h"
#define N_STEPS 91 /* it should relativity big and odd number 
                   // because we are using composite Simpson's rule integral
                   // and for accuracy we choose large number and the method
                   // requires odd number */

/* solving TOV equations. it approximates the properties of a NS.
// it only deals with inside of the star, since the outside is analytic.
// ->return value: solution of TOV equations for inside of a star with give EoS */
TOV_T *TOV_solution(TOV_T *const TOV)
{
  if (!TOV->N)/* if N is not determined */
    TOV->N = N_STEPS;
  const double Fac = 0.25;
  const unsigned MAX_iter = 1000;
  double h_cent_prev = 1;
  double h_cent_new;
  double a = 0,b = 0;
  double m = 0;/* ADM mass of NS for various central enthalpy */
  Flag_T increase  = YES;
  Flag_T bisection = NO;
  EoS_T *eos = 0;
  unsigned i,iter;
 
  pr_line_custom('=');
  printf("Solving TOV equations for %s ...\n",TOV->description);
  
  /* some initialization and preparation */
  iter = 0;
  h_cent_new = 1.5;
  TOV->m = alloc_double(TOV->N);
  TOV->r = alloc_double(TOV->N);
  TOV->h = alloc_double(TOV->N);
  TOV->p = alloc_double(TOV->N);
  
  /* find enthalpy at the center of NS such that 
  // the baryonic mass reaches the desired value.
  // the spirit of the approach is bisection method in root finders */
  while (!EQL(m,TOV->bar_m) && iter < MAX_iter)
  {
    TOV->h_cent = h_cent_new;
    solve_ODE_enthalpy_approach(TOV);
    m = calculate_baryonic_mass(TOV);
    
    /* as long as the mass is less than desired mass, increase enthalpy
    // to find the range of enthalpy */
    if (increase == YES)
    {
      if (LSS(m,TOV->bar_m))
      {
        h_cent_prev = h_cent_new;
        /* increase the enthalpy: */
        h_cent_new += h_cent_new*Fac;
      }
      else
      {
        increase  = NO;
        bisection = YES;/* active bisection */
        a = h_cent_prev;
        b = h_cent_new;
      }
    }
    /* now we know the range of enthalpy and we use bisection method
    // to pin down the correct value of the central enthalpy */
    if (bisection == YES)
    {
      if (GRT(m,TOV->bar_m))
      {
        b = h_cent_new;
      }
      else
      {
        a = h_cent_new;
      }
      h_cent_new = a+0.5*(b-a);
    }
    iter++;
  }
  
  /* if it finder messed up */
  if (!isfinite(m))
  {
    abortEr("TOV solution failed!\n");
  }
  /* if it needs more step to find the root, set the last value as baryonic mass */
  if (!EQL(m,TOV->bar_m))
  {
    fprintf(stderr,"TOV root finder for enteral enthalpy needs more steps!\n"
                   "The difference between current baryonic mass and desired one is = %e\n",m-TOV->bar_m);
    TOV->bar_m = m;
  }
  TOV->ADM_m = TOV->m[TOV->N-1];
  calculate_phi(TOV);
  calculate_ADM_and_Komar_mass(TOV);/* perform some tests */
  
  /* having known every thing, now populate pressure */
  eos = initialize_EoS();
  for (i = 0; i < TOV->N; ++i)
  {
    eos->h = TOV->h[i];
    TOV->p[i] = eos->pressure(eos);
  }
  
  isotropic_coords_transformation(TOV);
  
  /* print some informations about TOV */
  printf("TOV properties:\n");
  printf("--> NS radius (Schwarzschild Coords.) = %e\n",TOV->r[TOV->N-1]);
  printf("--> NS radius (Isotropic Coords.)     = %e\n",TOV->rbar[TOV->N-1]);
  printf("--> ADM mass                          = %e\n",TOV->ADM_m);
  printf("--> psi at the center                 = %e\n",TOV->psi[0]);
  printf("--> baryonic mass                     = %e\n",TOV->bar_m);
  printf("--> central enthalpy                  = %e\n",TOV->h[0]);
  printf("--> central pressure                  = %e\n",TOV->p[0]);
  eos->h = TOV->h_cent;
  printf("--> central energy density            = %e\n",eos->energy_density(eos));
  
  free_EoS(eos);
  
  printf("Solving TOV equations for %s ==> Done.\n",TOV->description);
  pr_clock();
  pr_line_custom('=');
  
  return TOV;
}

/* finding rbar and psi in isotropic coordinates. */
static void isotropic_coords_transformation(TOV_T *const TOV)
{
  TOV->rbar = alloc_double(TOV->N);
  TOV->psi  = alloc_double(TOV->N);
  double *const rbar = TOV->rbar;
  double *const psi  = TOV->psi;
  const double M = TOV->m[TOV->N-1];
  const double R = TOV->r[TOV->N-1];
  const double a = 1;
  const double b = TOV->h_cent;
  const double s = (b-a)/(TOV->N-1);/* s > 0 */
  const double c =  c_rbar(TOV);/* constant of integration in rbar */
  double t;/* independent variable h, we conventionally called it t */
  double r,m;/* interpolated values for r and m */
  Interpolation_T *interp_r = init_interpolation();
  Interpolation_T *interp_m = init_interpolation();
  unsigned i;
  
  interp_r->method         = "Natural_Cubic_Spline_1D";
  interp_r->N_cubic_spline_1d->f   = TOV->r;
  interp_r->N_cubic_spline_1d->x   = TOV->h;
  interp_r->N_cubic_spline_1d->N   = TOV->N;
  plan_interpolation(interp_r);
  
  interp_m->method         = "Natural_Cubic_Spline_1D";
  interp_m->N_cubic_spline_1d->f   = TOV->m;
  interp_m->N_cubic_spline_1d->x   = TOV->h;
  interp_m->N_cubic_spline_1d->N   = TOV->N;
  plan_interpolation(interp_m);
  
  /* initialization */
  rbar[0] = 0;/* rabr(h=h_cent) */
  rbar[1] = r_approx(TOV->h[1],TOV->h[0])*exp(-c);/* rabr(h=h_cent-s) */
  rbar[TOV->N-1] = 0.5*(R-M+sqrt(SQR(R)-2*M*R));/* rbar(h=1) */
  t = TOV->h[TOV->N-1];/* t = 1 */
  
  /* for all points */
  for (i = TOV->N-2; i >= 2; --i)
  {
    double k1,k2,k3,k4;/* variables for Runge-Kutta method of 4th order */
    
    interp_r->N_cubic_spline_1d->h = t;
    interp_m->N_cubic_spline_1d->h = t;
    r = execute_interpolation(interp_r);
    m = execute_interpolation(interp_m);
    k1 = s*drbar_dh(t,rbar[i+1],r,m);
    
    interp_r->N_cubic_spline_1d->h = t+s/2;
    interp_m->N_cubic_spline_1d->h = t+s/2;
    r = execute_interpolation(interp_r);
    m = execute_interpolation(interp_m);
    k2 = s*drbar_dh(t+s/2,rbar[i+1]+k1/2,r,m);
    k3 = s*drbar_dh(t+s/2,rbar[i+1]+k2/2,r,m);
    
    interp_r->N_cubic_spline_1d->h = t+s;
    interp_m->N_cubic_spline_1d->h = t+s;
    r = execute_interpolation(interp_r);
    m = execute_interpolation(interp_m);
    k4 = s*drbar_dh(t+s,rbar[i+1]+k3,r,m);
    
    /* updating the values */
    rbar[i] = rbar[i+1]+(k1+2*k2+2*k3+k4)/6;
    t = a+(TOV->N-1-i)*s;
  }
  free_interpolation(interp_r);
  free_interpolation(interp_m);
  
  /* for all points */
  psi[0] = sqrt(exp(c));
  for (i = 1; i < TOV->N; ++i)
    psi[i] = sqrt(TOV->r[i]/rbar[i]);
  
}

/* ->return value - the constant of intergrartion in drbar_dr. */
static double c_rbar(TOV_T *const TOV)
{
  Integration_T *I = init_integration();
  double *f = alloc_double(TOV->N);
  const double *const r = TOV->r,
               *const m = TOV->m,
               *const h = TOV->h;
  const double M = m[TOV->N-1];
  const double R = r[TOV->N-1];

  double c,s;
  unsigned i;
  
  f[0] = 0;
  for (i = 1; i < TOV->N; ++i)
  {
    s = 1-2*m[i]/r[i];
    f[i] = (1/sqrt(s)-1)/r[i]*dr_dh(h[i],r[i],m[i]);
  }
  
  I->type = "Composite Simpson's Rule 1D";
  I->Composite_Simpson_1D->b = 1;/* since the variable is enthalpy */
  I->Composite_Simpson_1D->a = TOV->h_cent;
  I->Composite_Simpson_1D->n = TOV->N;
  I->Composite_Simpson_1D->f = f;
  plan_integration(I);
  c = execute_integration(I);
  c -= log(0.5*(R-M+sqrt(SQR(R)-2*M*R))/R);
  
  free(f);
  free_integration(I);
  
  return c;
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
  const double tol = (TOV->h_cent-1)/TOV->N;/* this is due to the integration */
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
static double drbar_dh(const double h,const double rbar,const double r, const double m)
{
  double ret;
  
  ret = rbar/sqrt(1-2*m/r)/r*dr_dh(h,r,m);
  
  if (!isfinite(ret))
    abortEr("The interpolation failed due to the high oscillation of interpolant.\n"
            " One solution could be to lower the resolution.\n");
    
  return ret;
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
