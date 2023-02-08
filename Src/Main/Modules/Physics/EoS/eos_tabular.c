#include "eos_tabular.h"

//Implements thermodynamic functions for a tabular equation of state.

//Calculates pressure from enthalpy by tabular EOS.
double EoS_p_h_tab(EoS_T* const eos)
{
    if (LSSEQL(eos->h, eos->cubic_spline->h_floor)) { return 0.0; }
    
    double p;  
    Interpolation_T *const interp_s = eos->cubic_spline->interp_p;
  
    interp_s->N_cubic_spline_1d->h  = eos->h;
    p = execute_interpolation(interp_s);
    
    return (LSSEQL(p,0.) || p == DBL_MAX ? 0. : p);
}

//Calculates rest-mass density from enthalpy by tabular EOS.
double EoS_rho0_h_tab(EoS_T* const eos)
{
    if (LSSEQL(eos->h, eos->cubic_spline->h_floor)) { return 0.; }
  
    double rho0;  
    Interpolation_T *const interp_s = eos->cubic_spline->interp_rho0;
  
    interp_s->N_cubic_spline_1d->h  = eos->h;
    rho0 = execute_interpolation(interp_s);
  
    return (LSSEQL(rho0,0.) || rho0 == DBL_MAX ? 0. : rho0);
}

//Calculates energy density from enthalpy by tabular EOS.
double EoS_e_h_tab(EoS_T* const eos)
{
    if (LSSEQL(eos->h, eos->cubic_spline->h_floor))
    return 0.;

    double e;
    Interpolation_T *const interp_s = eos->cubic_spline->interp_e;
  
    interp_s->N_cubic_spline_1d->h  = eos->h;
    e = execute_interpolation(interp_s);
  
    return (LSSEQL(e,0.) || e == DBL_MAX ? 0. : e);
}

//Calculates derivative of energy density wrt enthalpy, with enthalpy as the independent variable.
//Uses central finite difference approximation from natural cubic spline interpolation of e(h).
double EoS_de_dh_h_tab(EoS_T* const eos)    //FIXME: Unfinished
{
    if (LSSEQL(eos->h, eos->cubic_spline->h_floor))
    return 0.;
    
    double h_copy = eos->h;
    double h_delta = 1e-5;                  //FIXME: Make dynamic var
    
    eos->h = h_copy - h_delta;
    double e0 = EoS_e_h_tab(eos);
    eos->h = h_copy + h_delta;
    double e1 = EoS_e_h_tab(eos);
    
    eos->h = h_copy;
    
    return (e1 - e0) / (2 * h_delta);
}

//Calculates specific internal energy in terms of enthalpy,
//using equation e0 = e / rho0 - 1
//Dependent upon p(h) and rho0(h)
double EoS_e0_h_tab(EoS_T* const eos)    //FIXME: Untested
{
    if (LSSEQL(eos->h, eos->cubic_spline->h_floor))
    return 0.;
    
    return (EoS_e_h_tab(eos) / EoS_rho0_h_tab(eos)) - 1.0;
}

//Calculates derivative of rest-mass density wrt enthalpy, with enthalpy as the independent variable.
//Uses central finite difference approximation from natural cubic spline interpolation of rho0(h).
double EoS_drho0_dh_h_tab(EoS_T* const eos)    //FIXME: Unfinished
{
    if (LSSEQL(eos->h, eos->cubic_spline->h_floor))
    return 0.;
    
    double h_copy = eos->h;
    double h_delta = 1e-5;                  //FIXME: Make dynamic var
    
    eos->h = h_copy - h_delta;
    double rho0_0 = EoS_rho0_h_tab(eos);
    eos->h = h_copy + h_delta;
    double rho0_1 = EoS_rho0_h_tab(eos);
    
    eos->h = h_copy;
    
    return (rho0_1 - rho0_0) / (2 * h_delta);
}












