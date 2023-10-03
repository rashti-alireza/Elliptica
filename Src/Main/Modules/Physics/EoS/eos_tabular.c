#include "eos_tabular.h"

//Implements thermodynamic functions for a tabular equation of state.

//FIXME: Add warning when enthalpy value out of table bounds.

//Calculates pressure from enthalpy by tabular EOS.
double EoS_p_h_tab(EoS_T* const eos)
{
    double h_copy = eos->h;
    if (LSS(eos->h, eos->spline->h_floor))
    {
        // Iff h < enthalpy floor, set h = enthalpy floor temporarily.
        eos->h = eos->spline->h_floor;
    }
    else if (GRT(eos->h, eos->spline->h_ceil))
    {
        // Iff h > enthalpy ceiling, set h = enthalpy ceiling temporarily.
        eos->h = eos->spline->h_ceil;
    }
    else if (GRT(eos->h, eos->spline->h_max) || LSS(eos->h, 0.90))
    {
        printf("ERROR: EoS_p_h_tab (eos_tabular.c): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->spline->h_floor, eos->spline->h_max);
        Error0("Exit");
        return 0.0;
    }
    
    double p;  
    Interpolation_T *const interp_s = eos->spline->interp_p;
    
    if (eos->spline->use_log_approach)
    {
      *interp_s->h = log(eos->h);
      //p = exp(execute_interpolation(interp_s)) * exp(-eos->spline->c_p);
      p = exp(execute_interpolation(interp_s)) - eos->spline->c_p;
    }
    else
    {
      *interp_s->h = eos->h;
      p = execute_interpolation(interp_s);
    }
    
    eos->h = h_copy;
    //double p_floor = 1.0E-12;
    double p_floor = 0.0;
    
    return (LSSEQL(p, p_floor) || p == DBL_MAX ? 0. : p);
}

//Calculates rest-mass density from enthalpy by tabular EOS.
double EoS_rho0_h_tab(EoS_T* const eos)
{
    double h_copy = eos->h;
    if (LSS(eos->h, eos->spline->h_floor))
    {
        // Iff h < enthalpy floor, set h = enthalpy floor temporarily.
        eos->h = eos->spline->h_floor;
    }
    else if (GRT(eos->h, eos->spline->h_ceil))
    {
        // Iff h > enthalpy ceiling, set h = enthalpy ceiling temporarily.
        eos->h = eos->spline->h_ceil;
    }
    else if (GRT(eos->h, eos->spline->h_max) || LSS(eos->h, 0.90))
    {
        printf("ERROR: EoS_p_h_tab (eos_tabular.c): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->spline->h_floor, eos->spline->h_max);
        Error0("Exit");
        return 0.0;
    }
  
    double rho0;  
    Interpolation_T *const interp_s = eos->spline->interp_rho0;
    
    if (eos->spline->use_log_approach)
    {
      *interp_s->h = log(eos->h);
      rho0 = exp(execute_interpolation(interp_s)) - eos->spline->c_e;
    }
    else
    {
      *interp_s->h = eos->h;
      rho0 = execute_interpolation(interp_s);
    }
    
    eos->h = h_copy;
    //double rho0_floor = 1.0E-12;
    double rho0_floor = 0.0;
    //if (LSSEQL(rho0, rho0_floor)) { rho0 = rho0_floor; }
    //return (LSSEQL(rho0,0.) || rho0 == DBL_MAX ? 0. : rho0);
    return (LSSEQL(rho0, rho0_floor) || rho0 == DBL_MAX ? 0. : rho0);
}

//Calculates energy density from enthalpy by tabular EOS.
double EoS_e_h_tab(EoS_T* const eos)
{
    double h_copy = eos->h;
    if (LSS(eos->h, eos->spline->h_floor))
    {
        // Iff h < enthalpy floor, set h = enthalpy floor temporarily.
        eos->h = eos->spline->h_floor;
    }
    else if (GRT(eos->h, eos->spline->h_ceil))
    {
        // Iff h > enthalpy ceiling, set h = enthalpy ceiling temporarily.
        eos->h = eos->spline->h_ceil;
    }
    else if (GRT(eos->h, eos->spline->h_max) || LSS(eos->h, 0.90))
    {
        printf("ERROR: EoS_p_h_tab (eos_tabular.c): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->spline->h_floor, eos->spline->h_max);
        Error0("Exit");
        return 0.0;
    }

    double e;
    Interpolation_T *const interp_s = eos->spline->interp_e;
    
    if (eos->spline->use_log_approach)
    {
      *interp_s->h = log(eos->h);
      e = exp(execute_interpolation(interp_s)) - eos->spline->c_e;
    }
    else
    {
      *interp_s->h = eos->h;
      e = execute_interpolation(interp_s);
    }
    
    eos->h = h_copy;
    //double e_floor = 1.0E-12;
    double e_floor = 0.0;
    //if (LSSEQL(e, e_floor)) { e = e_floor; }
    //return (LSSEQL(e,0.) || e == DBL_MAX ? 0. : e);
    return (LSSEQL(e,e_floor) || e == DBL_MAX ? 0. : e);
}

//Calculates specific internal energy in terms of enthalpy,
//using equation e0 = e / rho0 - 1
//Dependent upon p(h) and rho0(h)
double EoS_e0_h_tab(EoS_T* const eos)
{
    double h_copy = eos->h;
    if (LSS(eos->h, eos->spline->h_floor))
    {
        // Iff h < enthalpy floor, set h = enthalpy floor temporarily.
        eos->h = eos->spline->h_floor;
    }
    else if (GRT(eos->h, eos->spline->h_ceil))
    {
        // Iff h > enthalpy ceiling, set h = enthalpy ceiling temporarily.
        eos->h = eos->spline->h_ceil;
    }
    else if (GRT(eos->h, eos->spline->h_max) || LSS(eos->h, 0.90))
    {
        printf("ERROR: EoS_p_h_tab (eos_tabular.c): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->spline->h_floor, eos->spline->h_max);
        Error0("Exit");
        return 0.0;
    }
    
    double floor = 1.0E-15; //FIXME: Make dynamic parameter
    double rho0 = EoS_rho0_h_tab(eos);
    double e = EoS_e_h_tab(eos);
    
    eos->h = h_copy;
    return ((LSSEQL(e, floor) ? floor: e) / (LSSEQL(rho0, floor) ? floor : rho0)) - 1.0; 
}

//Calculates derivative of rest-mass density wrt enthalpy, with enthalpy as the independent variable.
double EoS_drho0_dh_h_tab(EoS_T* const eos)
{
    double h_copy = eos->h;
    //Check bounds for enthalpy
    if (LSS(eos->h, eos->spline->h_floor))
    {
        // Iff h < enthalpy floor, set h = enthalpy floor temporarily.
        eos->h = eos->spline->h_floor;
    }
    else if (GRT(eos->h, eos->spline->h_ceil))
    {
        // Iff h > enthalpy ceiling, set h = enthalpy ceiling temporarily.
        eos->h = eos->spline->h_ceil;
    }
    else if (GRT(eos->h, eos->spline->h_max) || LSS(eos->h, 0.90))
    {
        printf("ERROR: EoS_p_h_tab (eos_tabular.c): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->spline->h_floor, eos->spline->h_max);
        Error0("Exit");
        return 0.0;
    }
    
    double drho0dh;
    Interpolation_T *const interp_s = eos->spline->interp_rho0;
    interp_s->fd_derivative_order = 1;
    
    if (eos->spline->use_log_approach)
    {
      // Via chain rule: df/dh = (f(x)/x) * d(log(f))/d(log(x))
      // d(log(rho0))/d(log(h)):
      double dlog_log = finite_difference_Fornberg(eos->spline->h_log,
                             eos->spline->rho0_log,
                             log(eos->h),
                             1,
                             interp_s->fd_accuracy_order,
                             eos->spline->sample_size);
      *interp_s->h = eos->h;
      drho0dh = dlog_log * (EoS_rho0_h_tab(eos) + eos->spline->c_rho0) / eos->h;
    }
    else
    {
      *interp_s->h = eos->h;
      drho0dh = execute_derivative_interpolation(interp_s);
    }
    
    eos->h = h_copy;
    return (LSSEQL(drho0dh,0.) || drho0dh == DBL_MAX ? 0. : drho0dh);
}

//Calculates derivative of energy density wrt enthalpy, with enthalpy as the independent variable.
double EoS_de_dh_h_tab(EoS_T* const eos)
{
    double h_copy = eos->h;
    //Check bounds for enthalpy
    if (LSS(eos->h, eos->spline->h_floor))
    {
        // Iff h < enthalpy floor, set h = enthalpy floor temporarily.
        eos->h = eos->spline->h_floor;
    }
    else if (GRT(eos->h, eos->spline->h_ceil))
    {
        // Iff h > enthalpy ceiling, set h = enthalpy ceiling temporarily.
        eos->h = eos->spline->h_ceil;
    }
    else if (GRT(eos->h, eos->spline->h_max) || LSS(eos->h, 0.90))
    {
        printf("ERROR: EoS_p_h_tab (eos_tabular.c): enthalpy (%E) out of bounds (%E, %E).\n",
             eos->h, eos->spline->h_floor, eos->spline->h_max);
        Error0("Exit");
        return 0.0;
    }
    
    double dedh;
    Interpolation_T *const interp_s = eos->spline->interp_e;
    interp_s->fd_derivative_order = 1;
    
    if (eos->spline->use_log_approach)
    {
      // Via chain rule: df/dh = (f(x)/x) * d(log(f))/d(log(x))
      // d(log(e))/d(log(h)):
      double dlog_log = finite_difference_Fornberg(eos->spline->h_log,
                             eos->spline->e_log,
                             log(eos->h),
                             1,
                             interp_s->fd_accuracy_order,
                             eos->spline->sample_size);
      *interp_s->h = eos->h;
      dedh = dlog_log * (EoS_e_h_tab(eos) + eos->spline->c_e)/ eos->h;
    }
    else
    {
      *interp_s->h = eos->h;
      dedh = execute_derivative_interpolation(interp_s);
    }
    
    eos->h = h_copy;
    return (LSSEQL(dedh,0.) || dedh == DBL_MAX ? 0. : dedh);
}

///////////////Root finder approach////////////////////
// Definition of enthalpy h: (e + p)/rho0 - h = 0.
double EoS_enthalpy_def(void* eos, const double* const sol)
{
  // FIXME: See if there's a better way than eos_copy
  EoS_T* eos_copy = (EoS_T* const)eos;
  eos_copy->spline->rho0 = sol[0];
  Interpolation_T* const interp_e_rho0 = eos_copy->spline->interp_e_rho0;
  Interpolation_T* const interp_p_rho0 = eos_copy->spline->interp_p_rho0;
  *interp_e_rho0->h = eos_copy->spline->rho0;
  *interp_p_rho0->h = eos_copy->spline->rho0;
  
  double e = execute_interpolation(interp_e_rho0);
  double p = execute_interpolation(interp_p_rho0);
  
  /*
  printf("EoS_enthalpy_def:\n");
  printf("\t h == %E\n", eos_copy->h);
  printf("\t rho0 == %E\n", eos_copy->spline->rho0);
  printf("\t p == %E\n", p);
  printf("\t e == %E\n", e);
  printf("\t return == %E\n", ((e + p) / eos_copy->spline->rho0) - eos_copy->h);
  */
  
  if (eos_copy->h < Pgetd("NS_eos_enthalpy_floor")) { return 0; }
  return ((e + p) / eos_copy->spline->rho0) - eos_copy->h;
}

// Finds rest-mass density from enthalpy, using root finder method
// on enthalpy definition: (e(rho0) + p(rho0)) / rho0 - h = 0.
// Also sets value of rho0 inside EoS object for use in other
// thermodynamic functions.
double EoS_rho0_RF(EoS_T *const eos)
{
  Root_Finder_T* root_finder = (Root_Finder_T*)eos->spline->root_finder;
  double rho0 = *execute_root_finder(root_finder);
  eos->spline->rho0 = rho0;
  return rho0;
}

// Interpolates pressure from rest-mass density, using root-finder method.
double EoS_p_rho0_tab(EoS_T* const eos)
{
  // Determine rho0 using root finder, given enthalpy.
  EoS_rho0_RF(eos);
  // Having rho0, interpolate p(rho0).
  double p = execute_interpolation(eos->spline->interp_p_rho0);
  return p;
}

// Interpolates total energy density from rest-mass density,
// using root-finder method.
double EoS_e_rho0_tab(EoS_T* const eos)
{
  // Determine rho0 using root finder, given enthalpy.
  EoS_rho0_RF(eos);
  // Having rho0, interpolate e(rho0).
  double e = execute_interpolation(eos->spline->interp_e_rho0);
  return e;
}

// Interpolates specific internal energy from rest-mass density,
// using root-finder method.
// Equation: e0 = (e/rho0) - 1. 
double EoS_e0_rho0_tab(EoS_T* const eos)
{
  // Determine e and rho0.
  // (rho0 is set in eos->spline by EoS_e_rho0_tab.)
  double e = EoS_e_rho0_tab(eos);
  return (e / eos->spline->rho0) - 1;
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

// read thermo. vars. from a given table
void eos_tab_read_table(EoS_T* const eos)
{
  enum tab_format {
  UTab_format, // undef
  UTab_f1,     // "line,number_density,total_energy_density,pressure"
  UTab_f2,     // "rest_mass_density,specific_internal_energy,pressure"
  };
  
  Physics_T *const phys = eos->phys;
  enum tab_format format = UTab_format;
  FILE* table_file = Fopen(Gets(P_"table_path"),"r");
  double *p_tab    = alloc_double(EOS_MAX_NUM_ROW_TABLE);
  double *e_tab    = alloc_double(EOS_MAX_NUM_ROW_TABLE);
  double *rho0_tab = alloc_double(EOS_MAX_NUM_ROW_TABLE);
  double *h_tab    = 0;
  Uint num_tab_row = 0;
  double rho0_pnt = DBL_MAX;
  double e_pnt = DBL_MAX;
  double p_pnt = DBL_MAX;
  char line[EOS_MAX_NUM_COL_TABLE] = {'\0'};
  
  // open and read eos table
  if (strcmp_i(Gets(P_"table_format"), 
       "line,number_density,total_energy_density,pressure"))
  {
    format = UTab_f1;
    num_tab_row = 0;
    while (fgets(line, sizeof(line), table_file))
    {
      // skip empty lines
      if (strspn(line, " \t\r\n") == strlen(line))
      {
        continue;
      }
      // skip comment lines
      if (line[0] == '#')
      {
        continue;
      }
      // skip table num. of lines
      if (scanf(line,"%u") == 1)
      {
        continue;
      }
      // populate thermo. vars.
      Uint l;
      if (scanf(line,"%u %0.15e %0.15e %0.15e",&l, &rho0_pnt, &e_pnt, &p_pnt) == 4)
      {
        rho0_tab[num_tab_row] = rho0_pnt;
        e_tab[num_tab_row]    = e_pnt;
        p_tab[num_tab_row]    = p_pnt;
        num_tab_row++;
        if (num_tab_row >= EOS_MAX_NUM_ROW_TABLE)
        {
          Error0("number of table's row exceeds the max value. "
                 "please increase the macro EOS_MAX_NUM_ROW_TABLE.");
        }
      }
      else
      {
        Errors("Cannot read eos line: %s.",line);
      }
    }
  }
  else if (strcmp_i(Gets(P_"table_format"), 
                  "rest_mass_density,specific_internal_energy,pressure"))
  {
    format = UTab_f2;
    num_tab_row = 0;
    while (fgets(line, sizeof(line), table_file))
    {
      Uint l;
      // skip empty lines
      if (strspn(line, " \t\r\n") == strlen(line))
      {
        continue;
      }
      // skip comment lines
      if (line[0] == '#')
      {
        continue;
      }
      // skip table num. of lines
      if (scanf(line,"%u",&l) == 1)
      {
        continue;
      }
      // populate thermo. vars.
      if (scanf(line,"%0.15e %0.15e %0.15e",&rho0_pnt, &e_pnt, &p_pnt) == 3)
      {
        rho0_tab[num_tab_row] = rho0_pnt;
        e_tab[num_tab_row]    = e_pnt;
        p_tab[num_tab_row]    = p_pnt;
        num_tab_row++;
        if (num_tab_row >= EOS_MAX_NUM_ROW_TABLE)
        {
          Error0("number of table's row exceeds the max value. "
                 "please increase the macro EOS_MAX_NUM_ROW_TABLE.");
        }
      }
      else
      {
        Errors("Cannot read eos line: %s.",line);
      }
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
  Fclose(table_file);
  
  // adjust to the actual size
  p_tab = realloc(p_tab,num_tab_row*sizeof(*p_tab)); IsNull(p_tab);
  e_tab = realloc(p_tab,num_tab_row*sizeof(*e_tab)); IsNull(e_tab);
  rho0_tab = realloc(p_tab,num_tab_row*sizeof(*rho0_tab)); IsNull(rho0_tab);
  h_tab    = alloc_double(num_tab_row);
  
  // unit conversion
  if (strcmp_i(eos->unit,"compose") && format == UTab_f1)
  {
    /* compose format in columns
    *          [line number] [number density] [(total) energy density] [pressure].
    * in units               [1/fm^3]         [g/cm^3]                 [dyn/cm^2].
    * Calculates rest-mass density (total) energy density, 
    * and specific enthalpy, and converts to geometrized units.
    *
    * Physical constants from 2018 CODATA values.
    * G_const = 6.67430E-11;          // Gravitational constant in m^3/(kg s^2)
    * c_const = 299792458;            // Speed of light in m/s
    * M_const = 1.98841E30;           // Solar mass in kg
    * mn_const = 1.67492749804E-27;   // Neutron mass in kg
    * compose data is converted to SI units then geometrized.
    * Pressure conversion factor: G^3 * M_solar^2 / (10 * c^8)
    * Energy density conversion factor: G^3  * M_solar^2 / (10 * c^6)
    * Rest-mass conversion factor: 10^45 * neutron mass * G^3 * m_solar ^2 / (c^6)
    * (Rest-mass is calculated by multiplying the baryon density by the baryon mass.)
    */
    const double p_FACTOR    = 1.80162095578956E-39;
    const double rho0_FACTOR = 0.002712069678583313;
    const double e_FACTOR    = 1.619216164136643E-18;

    for (Uint i = 0; i < num_tab_row; i++)
    {
      p_tab[i]    *= p_FACTOR;
      rho0_tab[i] *= rho0_FACTOR;
      e_tab[i]    *= e_FACTOR;
    }
  }
  else if (strcmp_i(eos->unit,"geo") && format == UTab_f2)
  {
    for (Uint i = 0; i < num_tab_row; i++)
    {
      e_tab[i] = e_tab[i]*rho0_tab[i]+rho0_tab[i];
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
  
  // now everything is in correct format, set enthalpy
  for (Uint i = 0; i < num_tab_row; i++)
  {
    h_tab[i] = (p_tab[i] + e_tab[i]) / rho0_tab[i];
  }
  
  eos->spline->h_sample    = h_tab;
  eos->spline->p_sample    = p_tab;
  eos->spline->e_sample    = e_tab;
  eos->spline->rho0_sample = rho0_tab;
  eos->spline->sample_size = num_tab_row;
  
  h_tab = 0;
  p_tab = 0;
  rho0_tab = 0;
}

// set spline eos struct for Hermite interpolation method
void eos_tab_set_hermite(EoS_T* const eos)
{
  Physics_T *const phys = eos->phys;
  const Uint sample_size = eos->spline->sample_size;
  
  eos->pressure                 = EoS_p_h_tab;
  eos->energy_density           = EoS_e_h_tab;
  eos->rest_mass_density        = EoS_rho0_h_tab;
  eos->specific_internal_energy = EoS_e0_h_tab;
  eos->de_dh                    = EoS_de_dh_h_tab;
  eos->drho0_dh                 = EoS_drho0_dh_h_tab;
  eos->spline->h_floor = Getd(P_"enthalpy_floor");
  eos->spline->h_ceil  = Getd(P_"enthalpy_ceiling");

  // Fill arrays with log(values) if we're using log-log interpolation.
  if (Pcmps(Gets(P_"log_approach"),"yes"))
  {
    eos->spline->use_log_approach = 1;
    double *h_log = alloc_double(sample_size);
    double *p_log = alloc_double(sample_size);
    double *e_log = alloc_double(sample_size);
    double *rho0_log = alloc_double(sample_size);

    // shifting to avoid log(0)
    // TODO: are they the best number?
    eos->spline->c_p = 1E-3;
    eos->spline->c_e = 1E-3;
    eos->spline->c_rho0 = 1E-3;
    
    for (Uint i = 0; i < sample_size; i++)
    {
      p_log[i]    = log(eos->spline->p_sample[i]    + eos->spline->c_p);
      rho0_log[i] = log(eos->spline->rho0_sample[i] + eos->spline->c_rho0);
      e_log[i]    = log(eos->spline->e_sample[i]    + eos->spline->c_e);
      h_log[i]    = log(eos->spline->h_sample[i]);
    }
    
    eos->spline->h_log    = h_log;
    eos->spline->p_log    = p_log;
    eos->spline->e_log    = e_log;
    eos->spline->rho0_log = rho0_log;
    p_log = 0;
    e_log = 0;
    rho0_log = 0;
    h_log = 0;
  }

  // NOTE: we assume each is a function of the enthalpy h. */
  // p:
  Interpolation_T *interp_p = init_interpolation();
  interp_p->method = Gets(P_"Interpolation_Method");
  interp_p->Hermite_1d->fd_accuracy_order = (Uint)Geti(P_"fd_accuracy_order");
  interp_p->Hermite_1d->spline_order      = (Uint)Geti(P_"spline_order");
  interp_p->Hermite_1d->f = eos->spline->p_sample;
  interp_p->Hermite_1d->x = eos->spline->h_sample;
  if (eos->spline->use_log_approach)
  {
    interp_p->Hermite_1d->f = eos->spline->p_log;
    interp_p->Hermite_1d->x = eos->spline->h_log;
  }
  interp_p->Hermite_1d->N        = sample_size;
  interp_p->Hermite_1d->No_Warn  = 1;/* suppress warning */
  eos->spline->interp_p = interp_p;
  plan_interpolation(interp_p);
  
  // e:
  Interpolation_T *interp_e = init_interpolation();
  interp_e->method = Gets(P_"Interpolation_Method");
  interp_e->Hermite_1d->fd_accuracy_order = (Uint)Geti(P_"fd_accuracy_order");
  interp_e->Hermite_1d->spline_order      = (Uint)Geti(P_"spline_order");
  interp_e->Hermite_1d->f = eos->spline->e_sample;
  interp_e->Hermite_1d->x = eos->spline->h_sample;
  if (eos->spline->use_log_approach)
  {
    interp_e->Hermite_1d->f = eos->spline->e_log;
    interp_e->Hermite_1d->x = eos->spline->h_log;
  }
  interp_e->Hermite_1d->N = sample_size;
  interp_e->Hermite_1d->No_Warn = 1;/* suppress warning */
  eos->spline->interp_e = interp_e;
  plan_interpolation(interp_e);
  
  // rho0:
  Interpolation_T *interp_rho0 = init_interpolation();
  interp_rho0->method  = Gets(P_"Interpolation_Method");
  interp_rho0->Hermite_1d->fd_accuracy_order = (Uint)Geti(P_"fd_accuracy_order");
  interp_rho0->Hermite_1d->spline_order      = (Uint)Geti(P_"spline_order");
  interp_rho0->Hermite_1d->f = eos->spline->rho0_sample;
  interp_rho0->Hermite_1d->x = eos->spline->h_sample;
  if (eos->spline->use_log_approach)
  {
    interp_rho0->Hermite_1d->f = eos->spline->rho0_log;
    interp_rho0->Hermite_1d->x = eos->spline->h_log;
  }
  interp_rho0->Hermite_1d->N       = sample_size;
  interp_rho0->Hermite_1d->No_Warn = 1;/* suppress warning */
  eos->spline->interp_rho0 = interp_rho0;
  plan_interpolation(interp_rho0);
}
