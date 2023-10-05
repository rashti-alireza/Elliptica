/*
// Alireza Rashti
// June 2019
*/

#include "eos_main.h"

/* to be used for sampling purposes */
static const double enthalpy_initial = 1.0;/* initial h */
static const double enthalpy_final   = 2.0;/* final h */

/* tutorial:
// this is how we can calculate thermodynamic quantities
// like pressure, rest mass and energy density:
//
// EoS_T *eos = init_EoS(phys);// create and populate an EoS struct based on parameter file and physics
// double pressure, rest_mass_density, energy_density;
//
// eos->h = 1.5;// the enthalpy that we are interested to find the thermodynamics quantities at
// pressure          = eos->pressure(eos);// calculate pressure 
// energy_density    = eos->energy_density(eos);// calculate energy_density
// rest_mass_density = eos->rest_mass_density(eos);// calculate rest_mass_density
// free_EoS(eos);
*/


/* allocating memory and initializing EoS structure according
// to the parameter file.
// ->return value: initialize EoS structure. */
EoS_T *init_EoS(Physics_T *const phys)
{
  EoS_T *eos = calloc(1,sizeof(*eos));
  IsNull(eos);
  
  eos->phys  = phys;
  
  assert(phys);
  AssureType(phys->ctype == NS);
  
  populate_EoS(eos);/* populating EoS based on parameter file */
  
  return eos;
}

/* clean EoS stuct */
void free_EoS(EoS_T *s)
{
  if (!s)
    return;
    
  Free(s->K);
  Free(s->rho0_th);
  Free(s->h_th);
  Free(s->n);
  Free(s->a);
  Free(s->gamma);
  
  // spline
  Free(s->spline->h_sample);
  Free(s->spline->p_sample);
  Free(s->spline->e_sample);
  Free(s->spline->rho0_sample);
  Free(s->spline->h_log);
  Free(s->spline->p_log);
  Free(s->spline->e_log);
  Free(s->spline->rho0_log);
  
  free_interpolation(s->spline->interp_p);
  free_interpolation(s->spline->interp_e);
  free_interpolation(s->spline->interp_rho0);

  free(s);
}

/* populating EoS struct based on parameter file */
static void populate_EoS(EoS_T *const eos)
{
  Physics_T *const phys = eos->phys;
  /* populate eos struct */
  sprintf(eos->description,"%s",Gets(P_"description"));
  sprintf(eos->type,"%s",       Gets(P_"type")); 
  sprintf(eos->unit,"%s",       Gets(P_"unit"));
  
  /* if it is a pwp */
  if (strcmp_i(eos->type,"piecewise_polytropic") ||
      strcmp_i(eos->type,"pwp"))
  {
    if (!strcmp_i(eos->unit,"geo"))
    {
      Error0(NO_OPTION);
    }
    
    /* quantities in polytropic EoS */
    double *K       = 0;
    double *rho0_th = 0;
    double *gamma   = 0;
    Uint N;/* number of pieces in case of pwp */
    Uint i;
    
    /* NOTE: order matters */
    gamma   = read_EoS_in_parameter_file(Gets(P_"Gamma"),&N);
    K       = read_EoS_in_parameter_file(Gets(P_"K0"),0);/* this is K0 */
    rho0_th = read_EoS_in_parameter_file(Gets(P_"rho0_th"),0);
    
    /* check if rho is in increasing order. */
    if (!rho0_th)
      Error0("rho threshold must be specified.\n");
      
    for (i = 1; i < N-1; ++i)
      if (GRT(rho0_th[i-1],rho0_th[i]))
        Error0("rho0_th for piecewise polytropic EoS "
                "must be written in increasing order.\n");
  
    eos->N      = N;
    eos->K      = K;/* this is K0 now, below all other Ks will be found. */
    eos->rho0_th = rho0_th;
    eos->gamma = gamma;
    /* NOTE: order matters */
    fill_n(eos);
    fill_K(eos);/* now complete all other Ks. */
    fill_a(eos);
    fill_h_th(eos);
    eos->pressure                 = EoS_p_h_pwp;
    eos->energy_density           = EoS_e_h_pwp;
    eos->rest_mass_density        = EoS_rho0_h_pwp;
    eos->specific_internal_energy = EoS_e0_h_pwp;
    eos->de_dh    = EoS_de_dh_h_pwp;
    eos->drho0_dh = EoS_drho0_dh_h_pwp;
  } // end if (strcmp_i(eos->type,"piecewise_polytropic") ...)
  
  /* pwp eos's are generally C^0 continuous so we use 
  // natural cubic spline method to smooth them. the idea is 
  // taking a sample of thermodynamic variables, p(h),e(h),rho0(h), 
  // and then use a cubic spline fit to these data.
  // NOTE: since we don't set the slop at the beginning and end ofn the 
  // interval the thermo. vars might get negative for h ~ 1.
  // NOTE: the required params are:
  // "sample_size"   : set how many points will be selected from the eos.
  // "enthalpy_floor": set thermodynamics vars to zero if h < enthalpy_floor. */
  else if (strcmp_i(eos->type,"pwp_natural_cubic_spline"))
  {
    if (!strcmp_i(eos->unit,"geo"))
    {
      Error0(NO_OPTION);
    }
    /* quantities in polytropic EoS */
    double *K       = 0;
    double *rho0_th = 0;
    double *gamma   = 0;
    Uint N;/* number of pieces in case of pwp */
    Uint i;
   
    /* NOTE: order matters */
    gamma   = read_EoS_in_parameter_file(Gets(P_"Gamma"),&N);
    K       = read_EoS_in_parameter_file(Gets(P_"K0"),0);/* this is K0 */
    rho0_th = read_EoS_in_parameter_file(Gets(P_"rho0_th"),0);
    
    /* check if rho is in increasing order. */
    if (!rho0_th)
      Error0("rho threshold must be specified.\n");
      
    for (i = 1; i < N-1; ++i)
      if (GRT(rho0_th[i-1],rho0_th[i]))
        Error0("rho0_th for piecewise polytropic EoS "
                "must be written in increasing order.\n");
  
    eos->N        = N;
    eos->K        = K;/* K0 */
    eos->rho0_th  = rho0_th;
    eos->gamma    = gamma;
    /* NOTE: order matters */
    fill_n(eos);
    fill_K(eos);/* now complete all other Ks. */
    fill_a(eos);
    fill_h_th(eos);
    eos->pressure          = EoS_p_h_pwp;
    eos->energy_density    = EoS_e_h_pwp;
    eos->rest_mass_density = EoS_rho0_h_pwp;
    eos->de_dh             = EoS_de_dh_h_pwp;
    eos->drho0_dh	         = EoS_drho0_dh_h_pwp;
    
    /* use pwp to find the (p, e, rho0) values and then use spline 
    // to smooth them. */
    const double h_i    = enthalpy_initial;
    const double h_f    = enthalpy_final;
    const Uint sample_s = (Uint)Geti(P_"sample_size");/* number of sample points */
    double *h_sample    = alloc_double(sample_s);
    double *p_sample    = alloc_double(sample_s);
    double *e_sample    = alloc_double(sample_s);
    double *rho0_sample = alloc_double(sample_s);
    const double dh     = (h_f-h_i)/(sample_s-1.);/* assumed equispaced */
    
    /* for sampling we use analytic pwp eqs. */
    for (i = 0; i < sample_s; ++i)
    {
      eos->h         = h_sample[i] = h_i + i*dh;
      p_sample[i]    = eos->pressure(eos);
      e_sample[i]    = eos->energy_density(eos);
      rho0_sample[i] = eos->rest_mass_density(eos);
    }
    /* save samples: */
    eos->spline->sample_size = sample_s;
    eos->spline->h_sample    = h_sample;
    eos->spline->p_sample    = p_sample;
    eos->spline->e_sample    = e_sample;
    eos->spline->rho0_sample = rho0_sample;
    eos->spline->h_floor     = Getd(P_"enthalpy_floor");
    
    /* find and save spline coeffs for (p, e, rho0).
    // NOTE: we assume each is a function of the enthalpy h. */
    // p:
    Interpolation_T *interp_p             = init_interpolation();
    interp_p->method                      = "Natural_Cubic_Spline_1D";
    interp_p->N_cubic_spline_1d->f        = p_sample;
    interp_p->N_cubic_spline_1d->x        = h_sample;
    interp_p->N_cubic_spline_1d->N        = sample_s;
    interp_p->N_cubic_spline_1d->No_Warn  = 1;/* suppress warning */
    plan_interpolation(interp_p);
    eos->spline->interp_p = interp_p;
    
    // e:
    Interpolation_T *interp_e             = init_interpolation();
    interp_e->method                      = "Natural_Cubic_Spline_1D";
    interp_e->N_cubic_spline_1d->f        = e_sample;
    interp_e->N_cubic_spline_1d->x        = h_sample;
    interp_e->N_cubic_spline_1d->N        = sample_s;
    interp_e->N_cubic_spline_1d->No_Warn  = 1;/* suppress warning */
    plan_interpolation(interp_e);
    eos->spline->interp_e = interp_e;
    
    // rho0:
    Interpolation_T *interp_rho0              = init_interpolation();
    interp_rho0->method                       = "Natural_Cubic_Spline_1D";
    interp_rho0->N_cubic_spline_1d->f         = rho0_sample;
    interp_rho0->N_cubic_spline_1d->x         = h_sample;
    interp_rho0->N_cubic_spline_1d->N         = sample_s;
    interp_rho0->N_cubic_spline_1d->No_Warn   = 1;/* suppress warning */
    plan_interpolation(interp_rho0);
    eos->spline->interp_rho0            = interp_rho0;
    
    /* assign functions for (p, e, rho0) */
    eos->pressure                 = EoS_p_h_pwp_ncs;
    eos->energy_density           = EoS_e_h_pwp_ncs;
    eos->rest_mass_density        = EoS_rho0_h_pwp_ncs;
    eos->specific_internal_energy = EoS_e0_h_pwp; /* FIXME: it uses pwp */
    
    /* FIXME: for now use analytical calculations so it isn't continuous */
    eos->de_dh    = EoS_de_dh_h_pwp;
    eos->drho0_dh = EoS_drho0_dh_h_pwp;
    
    /* set to null for precaution */
    h_sample    = 0;
    p_sample    = 0;
    e_sample    = 0;
    rho0_sample = 0;
  }// end else if (strcmp_i(eos->type,"pwp_natural_cubic_spline"))
  
  else if (strcmp_i(eos->type,"polytropic") ||
           strcmp_i(eos->type,"polytrope"))
  {
    if (!strcmp_i(eos->unit,"geo"))
    {
      Error0(NO_OPTION);
    }
    /* quantities in polytropic EoS */
    double *K       = 0;
    double *gamma   = 0;
    Uint N;/* number of pieces in case of pwp */
  
    /* NOTE: order matters */ 
    gamma   = read_EoS_in_parameter_file(Gets(P_"Gamma"),&N);
    K       = read_EoS_in_parameter_file(Gets(P_"K0"),0);/* this is K0 */
    
    if (N != 1)
    {
      Error0("This EoS is not polytropic, there is more than one piece.\n");
    }
      
    eos->N     = N;
    eos->K     = K;
    eos->gamma = gamma;
    /* NOTE: order matters */
    fill_n(eos);
    fill_a(eos);
    eos->pressure                 = EoS_p_h_p;
    eos->energy_density           = EoS_e_h_p;
    eos->rest_mass_density        = EoS_rho0_h_p;
    eos->specific_internal_energy = EoS_e0_h_p;
    eos->de_dh                    = EoS_de_dh_h_p;
    eos->drho0_dh                 = EoS_drho0_dh_h_p;
  }// end else if (strcmp_i(eos->type,"polytropic") ...)
  
  else if (strcmp_i(eos->type, "tabular") || 
           strcmp_i(eos->type, "tab")     || 
           strcmp_i(eos->type, "table"))
  {
    eos_tab_read_table(eos);
    
    if (strcmp_i(Gets(P_"Interpolation_Method"), "Hermite1D"))
    {
      if (strcmp_i(Gets(P_"interpolation_use_log"),"yes"))
      {
        eos_tab_set_hermite_log(eos);
      }
      else
      {
        eos_tab_set_hermite(eos);
      }
    }
    else
    {
      Error0(NO_OPTION);
    }
  }// end else if (strcmp_i(eos->type, "tabular") ...)
  else
  {
    Error0(NO_JOB);
  }
}

/* filling n using n = 1/(gamma-1) */
static void fill_n(EoS_T *const eos)
{
  Uint i;
  
  assert(eos->N);
  eos->n = alloc_double(eos->N);
  for (i = 0; i < eos->N; ++i)
  {
    double gm1 = eos->gamma[i] - 1.;
    assert(!EQL(gm1,0.));
    
    eos->n[i] = 1./gm1;
  }
}

/* given the initial K (K0), fills other Ks (polytropic const) using K0 (first value).
// ref: eq 11. https://arxiv.org/pdf/0812.2163.pdf */
static void fill_K(EoS_T *const eos)
{
  Uint i;
  
  assert(eos->N);
  const double *const n    = eos->n;
  const double *const rho0 = eos->rho0_th;
  double *const K = eos->K = realloc(eos->K, eos->N*sizeof(*eos->K));
  IsNull(K);
  
  for (i = 1; i < eos->N; ++i)
  {
    K[i] = K[i-1]*pow(rho0[i],(n[i]-n[i-1])/(n[i-1]*n[i]));
  }
  
  //for(i = 0; i < eos->N; ++i) { printf("K%i == %E\n", i, K[i]); }
}

/* filling a by requiring the continuity of energy density */
static void fill_a(EoS_T *const eos)
{
  assert(eos->N);
  eos->a = alloc_double(eos->N);
  
  double *const a = eos->a,
         *const K = eos->K,
         *const gamma = eos->gamma,
         *const rho0_th = eos->rho0_th;
  Uint i;
         
  a[0] = 0;
  for (i = 1; i < eos->N; ++i)
    a[i] = a[i-1] + 
           K[i-1]*pow(rho0_th[i],gamma[i-1]-1)/(gamma[i-1]-1) -
           K[i]*pow(rho0_th[i],gamma[i]-1)/(gamma[i]-1);
           
}

/* filling a thresholds of enthalpy */
static void fill_h_th(EoS_T *const eos)
{
  assert(eos->a);
  eos->h_th = alloc_double(eos->N);
  
  double *const a = eos->a,
         *const K = eos->K,
         *const gamma = eos->gamma,
         *const rho0_th = eos->rho0_th,
         *const h_th = eos->h_th;
  Uint i;
  
  h_th[0] = 1;
  
  for (i = 1; i < eos->N; ++i)
    h_th[i] = 1+a[i]+(gamma[i]/(gamma[i]-1))*K[i]*pow(rho0_th[i],gamma[i]-1);
}


/* give parameter value for EoS, it parses the string
// and returns a pointer to the value found in the string. 
// also if N != 0, it counts the number of pieses in EoS.
// ->return value: value of parameter name and number of pieses in EoS. */
static double *read_EoS_in_parameter_file(const char *const par,Uint *const N)
{
  if (par == 0)
  {
    *N = 0;
    return 0;
  }
    
  double *v = 0;
  char str[MAX_STR],*sub_tok,*save;
  Uint i = 0;
  
  if (!check_format_s(par,"[?]"))
    Error0("K, rho0_th and Gamma in EoS must be written in square brackets.\n"
    "If there are multivalues, as in piecewise polytropic EoS, the values\n"
    "must be separated by a comma ','\n");
  
  /* parsing */
  strcpy(str,par);
  sub_tok = sub_s(str,'[',']',&save);/* => v1,v2,v3,... */
  sub_tok = tok_s(sub_tok,',',&save);/* sub_tok = v1 and save = v2,v3,... */
  if (sub_tok == 0)
    Errors("There is no value in %s.\n",par);
  
  v = realloc(v,(i+1)*sizeof(*v));
  IsNull(v);
  v[i] = atof(sub_tok);
  i++;
  while (sub_tok)
  {
    sub_tok = tok_s(0,',',&save);
    if (sub_tok)
    {
      v = realloc(v,(i+1)*sizeof(*v));
      IsNull(v);
      v[i] = atof(sub_tok);
      i++;
    }
  }
  
  if (N)
    *N = i;
  
  return v;
}
