#include "eos_tabular.h"

//Implements thermodynamic functions for a tabular equation of state.

//FIXME: Add warning when enthalpy value out of table bounds.

Uint get_sample_size(const char* const eos_file_name)
{
    //Takes name of tabular EOS file; returns the number of lines/EOS data points in the file.
    FILE* eos_file = fopen(eos_file_name, "r");
    Uint lines = 0;
    int character;
    
    while ((character = fgetc(eos_file)) != EOF) { if (character == '\n') { lines++; } }
    fclose(eos_file);
    
    return lines;
}

//Calculates pressure from enthalpy by tabular EOS.
double EoS_p_h_tab(EoS_T* const eos)
{
    if (LSS(eos->h, eos->cubic_spline->h_floor) || GRT(eos->h, eos->cubic_spline->h_max))
    {
        printf("ERROR: tabular p(h): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->cubic_spline->h_floor, eos->cubic_spline->h_max);
        Error0("Exit");
        return 0.0;
    }
    
    double p;  
    Interpolation_T *const interp_s = eos->cubic_spline->interp_p;
    
    if (eos->cubic_spline->use_log_approach)
    {
      *interp_s->h = log(eos->h);
      p = exp(execute_interpolation(interp_s));
    }
    else
    {
      *interp_s->h = eos->h;
      p = execute_interpolation(interp_s);
    }
    return (LSSEQL(p,0.) || p == DBL_MAX ? 0. : p);
}

//Calculates rest-mass density from enthalpy by tabular EOS.
double EoS_rho0_h_tab(EoS_T* const eos)
{
    if (LSS(eos->h, eos->cubic_spline->h_floor) || GRT(eos->h, eos->cubic_spline->h_max))
    {
        printf("ERROR: tabular rho0(h): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->cubic_spline->h_floor, eos->cubic_spline->h_max);
        Error0("Exit");
        return 0.0;
    }
  
    double rho0;  
    Interpolation_T *const interp_s = eos->cubic_spline->interp_rho0;
    
    if (eos->cubic_spline->use_log_approach)
    {
      *interp_s->h = log(eos->h);
      rho0 = exp(execute_interpolation(interp_s));
    }
    else
    {
      *interp_s->h = eos->h;
      rho0 = execute_interpolation(interp_s);
    }
    return (LSSEQL(rho0,0.) || rho0 == DBL_MAX ? 0. : rho0);
}

//Calculates energy density from enthalpy by tabular EOS.
double EoS_e_h_tab(EoS_T* const eos)
{
    if (LSS(eos->h, eos->cubic_spline->h_floor) || GRT(eos->h, eos->cubic_spline->h_max))
    {
        printf("ERROR: tabular e(h): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->cubic_spline->h_floor, eos->cubic_spline->h_max);
        Error0("Exit");
        return 0.0;
    }

    double e;
    Interpolation_T *const interp_s = eos->cubic_spline->interp_e;
    
    if (eos->cubic_spline->use_log_approach)
    {
      *interp_s->h = log(eos->h);
      e = exp(execute_interpolation(interp_s));
    }
    else
    {
      *interp_s->h = eos->h;
      e = execute_interpolation(interp_s);
    }
    return (LSSEQL(e,0.) || e == DBL_MAX ? 0. : e);
}

//Calculates specific internal energy in terms of enthalpy,
//using equation e0 = e / rho0 - 1
//Dependent upon p(h) and rho0(h)
double EoS_e0_h_tab(EoS_T* const eos)
{
    if (LSS(eos->h, eos->cubic_spline->h_floor) || GRT(eos->h, eos->cubic_spline->h_max))
    {
        printf("ERROR: e0(h): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->cubic_spline->h_floor, eos->cubic_spline->h_max);
        Error0("Exit");
        return 0.0;
    }
    
    double rho0_floor = 1.0E-7; //FIXME: Make dynamic parameter
    
    return (EoS_e_h_tab(eos) / (EoS_rho0_h_tab(eos) + rho0_floor)) - 1.0; 
}

//Calculates derivative of rest-mass density wrt enthalpy, with enthalpy as the independent variable.
double EoS_drho0_dh_h_tab(EoS_T* const eos)
{
    //Check bounds for enthalpy
    if (LSS(eos->h, eos->cubic_spline->h_floor) || GRT(eos->h, eos->cubic_spline->h_max))
    {
        printf("ERROR: drho0/dh(h): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->cubic_spline->h_floor, eos->cubic_spline->h_max);
         Error0("Exit");
        return 0.0;
    }
    
    double drho0dh;
    Interpolation_T *const interp_s = eos->cubic_spline->interp_rho0;
    interp_s->FDM_derivative = 1;
    
    if (eos->cubic_spline->use_log_approach)
    {
      drho0dh = FDM_Fornberg(eos->cubic_spline->h_sample,
                             eos->cubic_spline->rho0_sample,
                             eos->h,
                             1,
                             interp_s->finite_diff_order,
                             eos->cubic_spline->sample_size);
    }
    else
    {
      *interp_s->h = eos->h;
      drho0dh = execute_derivative_interpolation(interp_s);
    }
    return (LSSEQL(drho0dh,0.) || drho0dh == DBL_MAX ? 0. : drho0dh);
}

//Calculates derivative of energy density wrt enthalpy, with enthalpy as the independent variable.
double EoS_de_dh_h_tab(EoS_T* const eos)
{
    //Check bounds for enthalpy
    if (LSS(eos->h, eos->cubic_spline->h_floor) || GRT(eos->h, eos->cubic_spline->h_max))
    {
        printf("ERROR: de/dh(h): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->cubic_spline->h_floor, eos->cubic_spline->h_max);
         Error0("Exit");
        return 0.0;
    }
    
    double dedh;
    Interpolation_T *const interp_s = eos->cubic_spline->interp_e;
    interp_s->FDM_derivative = 1;
    
    if (eos->cubic_spline->use_log_approach)
    {
      dedh = FDM_Fornberg(eos->cubic_spline->h_sample,
                             eos->cubic_spline->e_sample,
                             eos->h,
                             1,
                             interp_s->finite_diff_order,
                             eos->cubic_spline->sample_size);
    }
    else
    {
      *interp_s->h = eos->h;
      dedh = execute_derivative_interpolation(interp_s);
    }
    return (LSSEQL(dedh,0.) || dedh == DBL_MAX ? 0. : dedh);
}

///////////////Root finder approach////////////////////
// Definition of enthalpy h: (e + p)/rho0 - h = 0.
double EoS_enthalpy_def(void* eos, const double* const sol)
{
  // FIXME: See if there's a better way than eos_copy
  EoS_T* eos_copy = (EoS_T* const)eos;
  eos_copy->cubic_spline->rho0 = sol[0];
  Interpolation_T* const interp_e_rho0 = eos_copy->cubic_spline->interp_e_rho0;
  Interpolation_T* const interp_p_rho0 = eos_copy->cubic_spline->interp_p_rho0;
  *interp_e_rho0->h = eos_copy->cubic_spline->rho0;
  *interp_p_rho0->h = eos_copy->cubic_spline->rho0;
  
  double e = execute_interpolation(interp_e_rho0);
  double p = execute_interpolation(interp_p_rho0);
  
  /*
  printf("EoS_enthalpy_def:\n");
  printf("\t h == %E\n", eos_copy->h);
  printf("\t rho0 == %E\n", eos_copy->cubic_spline->rho0);
  printf("\t p == %E\n", p);
  printf("\t e == %E\n", e);
  printf("\t return == %E\n", ((e + p) / eos_copy->cubic_spline->rho0) - eos_copy->h);
  */
  
  if (eos_copy->h < Pgetd("NS_eos_enthalpy_floor")) { return 0; }
  return ((e + p) / eos_copy->cubic_spline->rho0) - eos_copy->h;
}

// Finds rest-mass density from enthalpy, using root finder method
// on enthalpy definition: (e(rho0) + p(rho0)) / rho0 - h = 0.
// Also sets value of rho0 inside EoS object for use in other
// thermodynamic functions.
double EoS_rho0_RF(EoS_T *const eos)
{
  Root_Finder_T* root_finder = (Root_Finder_T*)eos->cubic_spline->root_finder;
  double rho0 = *execute_root_finder(root_finder);
  eos->cubic_spline->rho0 = rho0;
  return rho0;
}

// Interpolates pressure from rest-mass density, using root-finder method.
double EoS_p_rho0_tab(EoS_T* const eos)
{
  // Determine rho0 using root finder, given enthalpy.
  EoS_rho0_RF(eos);
  // Having rho0, interpolate p(rho0).
  double p = execute_interpolation(eos->cubic_spline->interp_p_rho0);
  return p;
}

// Interpolates total energy density from rest-mass density,
// using root-finder method.
double EoS_e_rho0_tab(EoS_T* const eos)
{
  // Determine rho0 using root finder, given enthalpy.
  EoS_rho0_RF(eos);
  // Having rho0, interpolate e(rho0).
  double e = execute_interpolation(eos->cubic_spline->interp_e_rho0);
  return e;
}

// Interpolates specific internal energy from rest-mass density,
// using root-finder method.
// Equation: e0 = (e/rho0) - 1. 
double EoS_e0_rho0_tab(EoS_T* const eos)
{
  // Determine e and rho0.
  // (rho0 is set in eos->cubic_spline by EoS_e_rho0_tab.)
  double e = EoS_e_rho0_tab(eos);
  return (e / eos->cubic_spline->rho0) - 1;
}

/////////////////////UNFINISHED//////////////////////
// Interpolates derivative of total energy density wrt enthalpy,
// using root finder method for rho0.
// Assumes thermodynamic relation de0 = -PdV
// (i.e. T ~ 0  and composition term ~ 0).
// Thus de0/drho0 = p / rho0^2.
// May be invalid for some EoS.
double EoS_de_dh_RF(EoS_T* const eos)
{
  UNUSED(eos);
  return 0;
}

double EoS_drho0_dh_RF(EoS_T* const eos)
{
  UNUSED(eos);
  return 0;
}


