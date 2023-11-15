/*
// Andrew Noe & Alireza Rashti
// October 2023
*/

#include "eos_tabular.h"

// Implements thermodynamic functions for a tabular equation of state.

// log for table data
#define Table_Log(x) log10(x)
// log inverse for table data(note exp10 is only part of GNU_SOURCE)
#define Table_Log_Inv(x) pow(10.0,x)

// assumes a generic Hermite interpolant in Table_Log, Table_Log(y) = Table_Log(y(Table_Log(h)))
// ->return: y(in non Table_Log format)
static double logy_of_logh_hermite(EoS_T* const eos, 
                                 Interpolation_T *const interp_s,
                                 const double y_shift/* shifting constant */,
                                 const double y_floor)
{
  const double h_floor = eos->spline->h_floor;
  const double h_ceil  = eos->spline->h_ceil;
  double h = eos->h;
  double y = DBL_MAX;
  
  // TODO: why this happens(if any)?
  h = (h <= h_floor) ? h_floor : h;
  h = (h >= h_ceil)  ? h_ceil  : h;
  // NOTE: we modify eos->h here:
  eos->h = h;
  interp_s->Hermite_1d->h = Table_Log(h);
  
  y = Table_Log_Inv(execute_interpolation(interp_s)) - y_shift;
  
  // TODO: DEBUG, why this happens(if any)?
  assert(y >= y_floor);
  
  return (y < y_floor) ? y_floor : y;
}

// assumes Table_Log(p) = hermite(Table_Log(h))
// ->return: p (pressure)
static double p_of_h_hermite_log(EoS_T* const eos)
{
  return
    logy_of_logh_hermite(eos, eos->spline->interp_p,
                         eos->spline->p_shift, eos->spline->p_floor);
}

// assumes Table_Log(e) = hermite(Table_Log(h))
// ->return: e (energy density)
static double e_of_h_hermite_log(EoS_T* const eos)
{
  return
    logy_of_logh_hermite(eos, eos->spline->interp_e,
                         eos->spline->e_shift, eos->spline->e_floor);
}

// assumes Table_Log(rho0) = hermite(Table_Log(h))
// ->return: rho0 (rest mass density)
static double rho0_of_h_hermite_log(EoS_T* const eos)
{
  return 
    logy_of_logh_hermite(eos, eos->spline->interp_rho0,
                         eos->spline->rho0_shift, eos->spline->rho0_floor);
}

// assumes rho0 = (e+p)/h
// ->return: rho0 (rest mass density)
static double rho0_e_p_h(EoS_T* const eos)
{
  return ( eos->pressure(eos) + eos->energy_density(eos) ) / eos->h;
}

// Calculates specific internal energy(e0) in terms of enthalpy,
// using equation e0 = e / rho0 - 1
// Dependent upon e(h) and rho0(h)
static double e0_of_e_and_rho0(EoS_T* const eos)
{
  double rho0 = eos->rest_mass_density(eos);
  double e    = eos->energy_density(eos);
  double e0   = e/rho0 - 1.;

  // the following should also capture when rho0 = 0.
  e0 = EQL(rho0,e) ? 0.0 : e0;

  return e0;
}

// dp/dh when we have interpolant Table_Log(p) = Hermite(Table_Log(h))
// chain rule: df/dh = ( f(h) / h ) * d(Table_Log(f)) / d(Table_Log(h)
static double dp_dh_hermite_log(EoS_T* const eos)
{
  const double f = eos->pressure(eos);
  const double h = eos->h;
  const double dlogf_dlogh = execute_1st_deriv_interpolation(eos->spline->interp_p);
  
  return (f / h ) * dlogf_dlogh;
}

// de/dh when we have interpolant Table_Log(e) = Hermite(Table_Log(h))
// chain rule: df/dh = ( f(h) / h ) * d(Table_Log(f)) / d(Table_Log(h)
static double de_dh_hermite_log(EoS_T* const eos)
{
  const double f = eos->energy_density(eos);
  const double h = eos->h;
  const double dlogf_dlogh = execute_1st_deriv_interpolation(eos->spline->interp_e);
  
  return (f / h ) * dlogf_dlogh;
}

// Calculates derivative of rest-mass density wrt enthalpy, 
// when drho0/dh = d{(e+p)/h}/dh.
static double drho0_dh_e_p_h(EoS_T* const eos)
{
  const double e     = eos->energy_density(eos);
  const double de_dh = eos->de_dh(eos);
  const double p     = eos->pressure(eos);
  const double dp_dh = dp_dh_hermite_log(eos);
  const double h     = eos->h;
  double drho0_dh    = ( de_dh + dp_dh - ( e + p ) / h ) / h;
  
  return drho0_dh;
}

// 1) read thermo. vars. from a given table
// 2) convert to geometrised units
// 3) populate pressure, rest mass density, total energy density, and enthalpy
void eos_tab_read_table(EoS_T* const eos)
{
  enum tab_format {
  UTab_format, // undef
  Tab_Format1, // "line,number_density,total_energy_density,pressure"
  Tab_Format2, // "rest_mass_density,specific_internal_energy,pressure"
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
    format = Tab_Format1;
    num_tab_row = 0;
    while (fgets(line, sizeof(line), table_file))
    {
      Uint l;
      // skip comment lines
      if (line[0] == '#')
      {
        continue;
      }
      // populate thermo. vars.
      // NOTE: it only reads the matching lines.
      if (sscanf(line,"%u %lf %lf %lf",&l, &rho0_pnt, &e_pnt, &p_pnt) == 4)
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
    }
  }
  else if (strcmp_i(Gets(P_"table_format"), 
                  "rest_mass_density,specific_internal_energy,pressure"))
  {
    format = Tab_Format2;
    num_tab_row = 0;
    while (fgets(line, sizeof(line), table_file))
    {
      // skip comment lines
      if (line[0] == '#')
      {
        continue;
      }
      // populate thermo. vars.
      if (sscanf(line,"%lf %lf %lf",&rho0_pnt, &e_pnt, &p_pnt) == 3)
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
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
  Fclose(table_file);
  
  if (num_tab_row == 0)
  {
    Error0("nothing read from the table.");
  }
  // adjust to the actual size
  p_tab = realloc(p_tab,num_tab_row*sizeof(*p_tab)); IsNull(p_tab);
  e_tab = realloc(e_tab,num_tab_row*sizeof(*e_tab)); IsNull(e_tab);
  rho0_tab = realloc(rho0_tab,num_tab_row*sizeof(*rho0_tab)); IsNull(rho0_tab);
  h_tab    = alloc_double(num_tab_row);
  
  // unit conversion
  if (strcmp_i(eos->unit,"compose") && format == Tab_Format1)
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
  else if (strcmp_i(eos->unit,"geo") && format == Tab_Format2)
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
void eos_tab_set_hermite_log(EoS_T* const eos)
{
  Physics_T *const phys = eos->phys;
  const Uint sample_size = eos->spline->sample_size;
  const Uint shift_row   = 3; // pick this row to shift the numbers
  double *h_log = alloc_double(sample_size);
  double *p_log = alloc_double(sample_size);
  double *e_log = alloc_double(sample_size);
  double *rho0_log = alloc_double(sample_size);

  eos->pressure                 = p_of_h_hermite_log;
  eos->energy_density           = e_of_h_hermite_log;
  eos->rest_mass_density        = rho0_e_p_h;// NOTE: using consistency eq.
  eos->specific_internal_energy = e0_of_e_and_rho0;// NOTE: using consistency eq.
  eos->de_dh                    = de_dh_hermite_log;
  eos->drho0_dh                 = drho0_dh_e_p_h;// NOTE: using consistency eq.

  eos->spline->use_log = 1;
  // TODO: one can set these values from the table?
  eos->spline->h_floor = Getd(P_"enthalpy_floor");
  eos->spline->h_ceil  = Getd(P_"enthalpy_ceiling");
  // TODO: are they fine?
  eos->spline->p_floor = 0.;
  eos->spline->e_floor = 0.;
  eos->spline->rho0_floor = 0.;
  
  // shifting to avoid Table_Log(0)
  eos->spline->p_shift    = eos->spline->p_sample[shift_row];
  eos->spline->e_shift    = eos->spline->e_sample[shift_row];
  eos->spline->rho0_shift = eos->spline->rho0_sample[shift_row];
  assert(!EQL(eos->spline->p_shift,0.0));
  assert(!EQL(eos->spline->e_shift,0.0));
  assert(!EQL(eos->spline->rho0_shift,0.0));
  
  for (Uint i = 0; i < sample_size; i++)
  {
    p_log[i]    = Table_Log(eos->spline->p_sample[i]    + eos->spline->p_shift);
    rho0_log[i] = Table_Log(eos->spline->rho0_sample[i] + eos->spline->rho0_shift);
    e_log[i]    = Table_Log(eos->spline->e_sample[i]    + eos->spline->e_shift);
    h_log[i]    = Table_Log(eos->spline->h_sample[i]);
  }
  eos->spline->h_log    = h_log;
  eos->spline->p_log    = p_log;
  eos->spline->e_log    = e_log;
  eos->spline->rho0_log = rho0_log;
  p_log = 0;
  e_log = 0;
  rho0_log = 0;
  h_log = 0;

  // NOTE: we assume each is a function of the enthalpy h. */
  // p:
  Interpolation_T *interp_p = init_interpolation();
  interp_p->method = "Hermite1D";
  interp_p->Hermite_1d->fd_accuracy_order = (Uint)Geti(P_"Hermite1D_FD_accuracy");
  interp_p->Hermite_1d->num_points = (Uint)Geti(P_"Hermite1D_num_points");
  interp_p->Hermite_1d->f = eos->spline->p_log;
  interp_p->Hermite_1d->x = eos->spline->h_log;
  interp_p->Hermite_1d->N = sample_size;
  //interp_p->Hermite_1d->No_Warn  = 1;/* suppress warning */
  eos->spline->interp_p = interp_p;
  plan_interpolation(interp_p);
  
  // e:
  Interpolation_T *interp_e = init_interpolation();
  interp_e->method = "Hermite1D";
  interp_e->Hermite_1d->fd_accuracy_order = (Uint)Geti(P_"Hermite1D_FD_accuracy");
  interp_e->Hermite_1d->num_points = (Uint)Geti(P_"Hermite1D_num_points");
  interp_e->Hermite_1d->f = eos->spline->e_log;
  interp_e->Hermite_1d->x = eos->spline->h_log;
  interp_e->Hermite_1d->N = sample_size;
  //interp_e->Hermite_1d->No_Warn = 1;/* suppress warning */
  eos->spline->interp_e = interp_e;
  plan_interpolation(interp_e);
  
  // rho0:
  Interpolation_T *interp_rho0 = init_interpolation();
  interp_rho0->method = "Hermite1D";
  interp_rho0->Hermite_1d->fd_accuracy_order = (Uint)Geti(P_"Hermite1D_FD_accuracy");
  interp_rho0->Hermite_1d->num_points = (Uint)Geti(P_"Hermite1D_num_points");
  interp_rho0->Hermite_1d->f = eos->spline->rho0_log;
  interp_rho0->Hermite_1d->x = eos->spline->h_log;
  interp_rho0->Hermite_1d->N = sample_size;
  ///interp_rho0->Hermite_1d->No_Warn = 1;/* suppress warning */
  eos->spline->interp_rho0 = interp_rho0;
  plan_interpolation(interp_rho0);
}

// set spline eos struct for Hermite interpolation method
void eos_tab_set_hermite(EoS_T* const eos)
{
  Physics_T *const phys = eos->phys;
  const Uint sample_size = eos->spline->sample_size;
  
  // TODO: if needed in future set these:
  Error0(NO_OPTION);
  eos->pressure                 = 0;
  eos->energy_density           = 0;
  eos->rest_mass_density        = 0;
  eos->specific_internal_energy = 0;
  eos->de_dh                    = 0;
  eos->drho0_dh                 = 0;
  eos->spline->h_floor = Getd(P_"enthalpy_floor");
  eos->spline->h_ceil  = Getd(P_"enthalpy_ceiling");
  // TODO: are they fine?
  eos->spline->p_floor = 0.;
  eos->spline->e_floor = 0.;
  eos->spline->rho0_floor = 0.;
  
  // NOTE: we assume each is a function of the enthalpy h. */
  // p:
  Interpolation_T *interp_p = init_interpolation();
  interp_p->method = "Hermite1D";
  interp_p->Hermite_1d->fd_accuracy_order = (Uint)Geti(P_"Hermite1D_FD_accuracy");
  interp_p->Hermite_1d->num_points = (Uint)Geti(P_"Hermite1D_num_points");
  interp_p->Hermite_1d->f = eos->spline->p_sample;
  interp_p->Hermite_1d->x = eos->spline->h_sample;
  interp_p->Hermite_1d->N = sample_size;
  //interp_p->Hermite_1d->No_Warn  = 1;/* suppress warning */
  eos->spline->interp_p = interp_p;
  plan_interpolation(interp_p);
  
  // e:
  Interpolation_T *interp_e = init_interpolation();
  interp_e->method = "Hermite1D";
  interp_e->Hermite_1d->fd_accuracy_order = (Uint)Geti(P_"Hermite1D_FD_accuracy");
  interp_e->Hermite_1d->num_points = (Uint)Geti(P_"Hermite1D_num_points");
  interp_e->Hermite_1d->f = eos->spline->e_sample;
  interp_e->Hermite_1d->x = eos->spline->h_sample;
  interp_e->Hermite_1d->N = sample_size;
  //interp_e->Hermite_1d->No_Warn = 1;/* suppress warning */
  eos->spline->interp_e = interp_e;
  plan_interpolation(interp_e);
  
  // rho0:
  Interpolation_T *interp_rho0 = init_interpolation();
  interp_rho0->method = "Hermite1D";
  interp_rho0->Hermite_1d->fd_accuracy_order = (Uint)Geti(P_"Hermite1D_FD_accuracy");
  interp_rho0->Hermite_1d->num_points = (Uint)Geti(P_"Hermite1D_num_points");
  interp_rho0->Hermite_1d->f = eos->spline->rho0_sample;
  interp_rho0->Hermite_1d->x = eos->spline->h_sample;
  interp_rho0->Hermite_1d->N = sample_size;
  //interp_rho0->Hermite_1d->No_Warn = 1;/* suppress warning */
  eos->spline->interp_rho0 = interp_rho0;
  plan_interpolation(interp_rho0);
}
