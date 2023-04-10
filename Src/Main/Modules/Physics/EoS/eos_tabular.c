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
    
    *interp_s->h = eos->h;
    p = execute_interpolation(interp_s);
    return (LSSEQL(p,0.) || p == DBL_MAX ? 0. : p);
    
    /*
    if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
    { interp_s->N_cubic_spline_1d->h  = eos->h; }
    else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline_1D"))
    { interp_s->H_cubic_spline_1d->h = eos->h; }
    else if (strstr_i(interp_s->method, "Clamped_Cubic_Spline_1D"))
    { interp_s->C_cubic_spline_1d->h = eos->h; }
    */
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
    
    
    *interp_s->h = eos->h;
    rho0 = execute_interpolation(interp_s);
    return (LSSEQL(rho0,0.) || rho0 == DBL_MAX ? 0. : rho0);
    
    /*
    if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
    { interp_s->N_cubic_spline_1d->h  = eos->h; }
    else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline_1D"))
    { interp_s->H_cubic_spline_1d->h = eos->h; }
    else if (strstr_i(interp_s->method, "Clamped_Cubic_Spline_1D"))
    { interp_s->C_cubic_spline_1d->h = eos->h; }
    */
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
    
    *interp_s->h = eos->h;
    e = execute_interpolation(interp_s);
    return (LSSEQL(e,0.) || e == DBL_MAX ? 0. : e);
    
    /*
    if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
    { interp_s->N_cubic_spline_1d->h  = eos->h; }
    else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline_1D"))
    { interp_s->H_cubic_spline_1d->h = eos->h; }
    else if (strstr_i(interp_s->method, "Clamped_Cubic_Spline_1D"))
    { interp_s->C_cubic_spline_1d->h = eos->h; }
    */
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
    //Interpolation_T *const interp_s = (Interpolation_T*)eos->cubic_spline->interp_rho0;
    
    *interp_s->h = eos->h;
    interp_s->FDM_derivative = 1;
    drho0dh = execute_derivative_interpolation(interp_s);
    return (LSSEQL(drho0dh,0.) || drho0dh == DBL_MAX ? 0. : drho0dh);
   
    /*
    if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
    { interp_s->N_cubic_spline_1d->h  = eos->h; }
    else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline_1D"))
    { interp_s->H_cubic_spline_1d->h = eos->h; }
    else if (strstr_i(interp_s->method, "Clamped_Cubic_Spline_1D"))
    { interp_s->C_cubic_spline_1d->h = eos->h; }
    */
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
    *interp_s->h = eos->h;
    interp_s->FDM_derivative = 1;
    dedh = execute_interpolation(interp_s);
    
    return (LSSEQL(dedh,0.) || dedh == DBL_MAX ? 0. : dedh);
}







