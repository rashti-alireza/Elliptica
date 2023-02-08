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
    //char character;
    //while((character = fgetc(eos_file)) != EOF) { if (character == '\n') { lines++; } }
    
    //fclose(eos_file);
    //return lines;
}

/*Removed due to removal of finite difference version of de/dh and drho0/dh.
Uint find_h_interval(EoS_T *const eos)
{
    //Finds spline interval for enthalpy value (i.e. finds j such that hj <= h < hj+1).
    //WARNING: Does not perform bounds checking for enthalpy value,
    //will break if enthalpy is out of bounds.
    Uint interval = 0;
    Interpolation_T *const interpolation = (Interpolation_T*)eos->cubic_spline->interp_e;
    
    if (GRTEQL(eos->h, eos->cubic_spline->h_max)) { return eos->cubic_spline->sample_size; }
    
    while (GRT(eos->h, interpolation->N_cubic_spline_1d->x[interval+1])) { interval++; }
    printf("Enthalpy (%E < %E < %E) interval: %i\n", 
        interpolation->N_cubic_spline_1d->x[interval],
        eos->h, interpolation->N_cubic_spline_1d->x[interval+1], interval);
    
    return interval;
}
*/

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
  
    interp_s->N_cubic_spline_1d->h  = eos->h;
    p = execute_interpolation(interp_s);
    
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
  
    interp_s->N_cubic_spline_1d->h  = eos->h;
    rho0 = execute_interpolation(interp_s);
  
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
  
    interp_s->N_cubic_spline_1d->h  = eos->h;
    e = execute_interpolation(interp_s);
  
    return (LSSEQL(e,0.) || e == DBL_MAX ? 0. : e);
}

//Calculates derivative of energy density wrt enthalpy, with enthalpy as the independent variable.
//Uses derivative of natural cubic spline:
//  e(h) ~ aj + bj(h - hj) + cj(h-hj)^2 + dj(h-hj)^3, so
//  e'(h) ~ bj*h - 2cj(h-hj) + 3dj(h-hj)^2, where j is the spline interval.
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
    
    Interpolation_T *const interpolation = (Interpolation_T*)eos->cubic_spline->interp_e;
    
    //Finds the spline interval
    Uint h_interval;
    for (h_interval=0; h_interval<eos->cubic_spline->sample_size; ++h_interval)
    {
        if (GRTEQL(eos->h,interpolation->N_cubic_spline_1d->x[h_interval]) 
            && LSSEQL(eos->h,interpolation->N_cubic_spline_1d->x[h_interval+1]))
        {
            break;
        }
    }
    
    double bj = interpolation->N_cubic_spline_1d->b[h_interval];
    double cj = interpolation->N_cubic_spline_1d->c[h_interval];
    double dj = interpolation->N_cubic_spline_1d->d[h_interval];
    double hj = interpolation->N_cubic_spline_1d->x[h_interval];

    return bj + 2*cj*(eos->h-hj) + 3*dj*(eos->h-hj) * (eos->h-hj);
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
//Uses derivative of natural cubic spline:
//  e(h) ~ aj + bj(h - hj) + cj(h-hj)^2 + dj(h-hj)^3, so
//  e'(h) ~ bj*h - 2cj(h-hj) + 3dj(h-hj)^2, where j is the spline interval.
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
    
    Interpolation_T *const interpolation = (Interpolation_T*)eos->cubic_spline->interp_rho0;
    
    //Finds the spline interval
    Uint h_interval;
    for (h_interval=0; h_interval<eos->cubic_spline->sample_size; ++h_interval)
    {
        if (GRTEQL(eos->h,interpolation->N_cubic_spline_1d->x[h_interval]) 
            && LSSEQL(eos->h,interpolation->N_cubic_spline_1d->x[h_interval+1]))
        {
            break;
        }
    }
    
    double bj = interpolation->N_cubic_spline_1d->b[h_interval];
    double cj = interpolation->N_cubic_spline_1d->c[h_interval];
    double dj = interpolation->N_cubic_spline_1d->d[h_interval];
    double hj = interpolation->N_cubic_spline_1d->x[h_interval];

    return bj + 2*cj*(eos->h-hj) + 3*dj*(eos->h-hj) * (eos->h-hj);
    /*
    if (LSSEQL(eos->h, eos->cubic_spline->h_floor)) { return 0.0; }
    
    double h_copy = eos->h;
    double h_delta = 1e-5;                  //FIXME: Make dynamic var
    
    eos->h = h_copy - h_delta;
    double rho0_0 = EoS_rho0_h_tab(eos);
    eos->h = h_copy + h_delta;
    double rho0_1 = EoS_rho0_h_tab(eos);
    
    eos->h = h_copy;
    
    return (rho0_1 - rho0_0) / (2 * h_delta);
    */
}



