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
// clean_EoS(&eos);
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
void clean_EoS(EoS_T **eos)
{
  EoS_T *s = *eos;
  
  if (!s)
    return;
    
  free(s->K);
  free(s->rho_th);
  free(s->n);
  free(s->a);
  free(s);
}

/* populating EoS struct based on parameter file */
static void populate_EoS(EoS_T *const eos)
{
}

/* based on type of EoS and parameter file determines which function
// should be used for calculating themodynamics quantities like pressure and etc. */
static void plan_EoS(EoS_T *const eos)
{
  if (strstr_i(GetParameterS("EoS_type"),"piecewise_polytropic") ||
      strstr_i(GetParameterS("EoS_type"),"pwp"))
  {
    eos->pressure          = EoS_p_h_pwp;
    eos->energy_density    = EoS_e_h_pwp;
    eos->rest_mass_density = EoS_rho_h_pwp;
  }
  else
    abortEr(NO_OPTION);
}

