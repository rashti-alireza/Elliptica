/*
// Alireza Rashti
// June 2019
*/

#include "EoS.h"


/* tutorial:
// this is how we can calculate thermodynamic quantities
// like pressure, rest mass and energy density:
//
// EoS_T *eos = initialize_EoS();// create and populate an EoS struct based on parameter file
// double pressure, rest_mass_density, energy_density;
//
// eos->h = 10;// the enthalpy that we are interested to find the thermodynamics quantities at
// pressure          = eos->pressure(eos);// calculate pressure 
// energy_density    = eos->energy_density(eos);// calculate energy_density
// rest_mass_density = eos->rest_mass_density(eos);// calculate rest_mass_density
// free_EoS(&eos);
*/


/* allocating memory and initializing EoS structure according
// to the parameter file.
// ->return value: initialize EoS structure. */
EoS_T *initialize_EoS(void)
{
  EoS_T *eos = calloc(1,sizeof(*eos));
  pointerEr(eos);
  
  populate_EoS(eos);/* populating EoS based on parameter file */
  
  return eos;
}

/* clean EoS stuct */
void free_EoS(EoS_T *s)
{
  if (!s)
    return;
    
  if (s->K)      free(s->K);
  if (s->rho_th) free(s->rho_th);
  if (s->h_th)   free(s->h_th);
  if (s->n)      free(s->n);
  if (s->a)      free(s->a);
  if (s->gamma)  free(s->gamma);
  free(s);
}

/* populating EoS struct based on parameter file */
static void populate_EoS(EoS_T *const eos)
{
  double *K,*rho_th,*gamma;/* quantities in polytropic EoS */
  unsigned N;/* number of pieces in case of pwp */
  unsigned i;
  const char *par;
  
  /* populate eos struct */
  par = PgetsEZ("EoS_description");
  if (par)
    strcpy(eos->description,par);
  strcpy(eos->type,Pgets_E("EoS_type")); 
  strcpy(eos->unit,Pgets_E("EoS_unit"));
  
  K      = read_EoS_in_parameter_file(Pgets_E("EoS_K"),&N);
  rho_th = read_EoS_in_parameter_file(PgetsEZ("EoS_rho_th"),0);
  gamma  = read_EoS_in_parameter_file(Pgets_E("EoS_Gamma"),0);

  /* if the units are geometrised units */
  if (strcmp_i(eos->unit,"geo"))
  {
    /* if it is pwp */
    if (strcmp_i(eos->type,"piecewise_polytropic") ||
        strcmp_i(eos->type,"pwp"))
    {    
      /* check if rho is in increasing order. */
      if (!rho_th)
        abortEr("rho threshold must be specified.\n");
        
      for (i = 1; i < N-1; ++i)
        if (GRT(rho_th[i-1],rho_th[i]))
          abortEr("rho_th for piecewise polytropic EoS "
                  "must be written in increasing order.\n");
    
      eos->N      = N;
      eos->K      = K;
      eos->rho_th = rho_th;
      eos->gamma = gamma;
      fill_n(eos);
      fill_a(eos);
      fill_h_th(eos);/* this depends on a, so put it in the last */
      eos->pressure          = EoS_p_h_pwp;
      eos->energy_density    = EoS_e_h_pwp;
      eos->rest_mass_density = EoS_rho_h_pwp;
      eos->de_dh             = EoS_de_dh_h_pwp;
      eos->drho_dh	     = EoS_drho_dh_h_pwp;
    }
    else if (strcmp_i(eos->type,"polytropic") ||
             strcmp_i(eos->type,"polytrop")   ||
             strcmp_i(eos->type,"p"))
    {
      if (N != 1)
        abortEr("This EoS is not polytropic, there is more than one piece.\n");
        
      eos->N      = N;
      eos->K      = K;
      eos->gamma = gamma;
      fill_n(eos);
      fill_a(eos);
      eos->pressure          = EoS_p_h_p;
      eos->energy_density    = EoS_e_h_p;
      eos->rest_mass_density = EoS_rho_h_p;
      eos->de_dh             = EoS_de_dh_h_p;
      eos->drho_dh	     = EoS_drho_dh_h_p;
    }
    else
      abortEr(NO_JOB);
    
  }
  else
    abortEr(NO_JOB);
  
}

/* filling n using n = 1/(gamma-1) */
static void fill_n(EoS_T *const eos)
{
  unsigned i;
  
  assert(eos->N);
  eos->n = alloc_double(eos->N);
  for (i = 0; i < eos->N; ++i)
    eos->n[i] = 1/(eos->gamma[i]-1);
}

/* filling a by requiring the continuity of energy density */
static void fill_a(EoS_T *const eos)
{
  assert(eos->N);
  eos->a = alloc_double(eos->N);
  
  double *const a = eos->a,
         *const K = eos->K,
         *const gamma = eos->gamma,
         *const rho_th = eos->rho_th;
  unsigned i;
         
  a[0] = 0;
  for (i = 1; i < eos->N; ++i)
    a[i] = a[i-1] + 
           K[i-1]*pow(rho_th[i],gamma[i-1]-1)/(gamma[i-1]-1) -
           K[i]*pow(rho_th[i],gamma[i]-1)/(gamma[i]-1);
           
}

/* filling a thresholds of enthalpy */
static void fill_h_th(EoS_T *const eos)
{
  assert(eos->a);
  eos->h_th = alloc_double(eos->N);
  
  double *const a = eos->a,
         *const K = eos->K,
         *const gamma = eos->gamma,
         *const rho_th = eos->rho_th,
         *const h_th = eos->h_th;
  unsigned i;
  
  h_th[0] = 1;
  
  for (i = 1; i < eos->N; ++i)
    h_th[i] = 1+a[i]+(gamma[i]/(gamma[i]-1))*K[i]*pow(rho_th[i],gamma[i]-1);
}


/* give parameter value for EoS, it parses the string
// and returns a pointer to the value found in the string. 
// also if N != 0, it counts the number of pieses in EoS.
// ->return value: value of parameter name and number of pieses in EoS. */
static double *read_EoS_in_parameter_file(const char *const par,unsigned *const N)
{
  /* eg. if we have a single polytropic eos, then we don't have rho_th */
  if (par == 0)
    return 0;
    
  double *v = 0;
  char str[MAX_STR],*sub_tok,*save;
  unsigned i = 0;
  
  if (!check_format_s(par,"[?]"))
    abortEr("K, rho_th and Gamma in EoS must be written in square brackets.\n"
    "If there are multivalues, as in piecewise polytropic EoS, the values\n"
    "must be separated by a comma ','\n");
  
  /* parsing */
  strcpy(str,par);
  sub_tok = sub_s(str,'[',']',&save);/* => v1,v2,v3,... */
  sub_tok = tok_s(sub_tok,',',&save);/* sub_tok = v1 and save = v2,v3,... */
  if (sub_tok == 0)
    abortEr_s("There is no value in %s.\n",par);
  
  v = realloc(v,(i+1)*sizeof(*v));
  pointerEr(v);
  v[i] = atof(sub_tok);
  i++;
  while (sub_tok)
  {
    sub_tok = tok_s(0,',',&save);
    if (sub_tok)
    {
      v = realloc(v,(i+1)*sizeof(*v));
      pointerEr(v);
      v[i] = atof(sub_tok);
      i++;
    }
  }
  
  if (N)
    *N = i;
  
  return v;
}
