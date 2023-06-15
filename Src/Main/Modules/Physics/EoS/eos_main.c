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
    
  if (s->K)      free(s->K);
  if (s->rho0_th) free(s->rho0_th);
  if (s->h_th)   free(s->h_th);
  if (s->n)      free(s->n);
  if (s->a)      free(s->a);
  if (s->gamma)  free(s->gamma);
  
  if (strstr_i(PgetsEZ(P_"Approach"), "Root_Finder"))
  {
    free_interpolation(s->cubic_spline->interp_e_rho0);
    free_interpolation(s->cubic_spline->interp_p_rho0);
    Root_Finder_T* rf_copy = (Root_Finder_T*)(s->cubic_spline->root_finder);
    free_root_finder(rf_copy);
  }
  
  /* free cubic spline struct */
  Free(s->cubic_spline->h_sample);
  Free(s->cubic_spline->p_sample);
  Free(s->cubic_spline->e_sample);
  Free(s->cubic_spline->rho0_sample);
  set_interp_alloc_mem_flag(s->cubic_spline->interp_p, 0);
  set_interp_alloc_mem_flag(s->cubic_spline->interp_e, 0);
  set_interp_alloc_mem_flag(s->cubic_spline->interp_rho0, 0);
  free_interpolation(s->cubic_spline->interp_p);
  free_interpolation(s->cubic_spline->interp_e);
  free_interpolation(s->cubic_spline->interp_rho0);
  
  free(s);
}

/* populating EoS struct based on parameter file */
static void populate_EoS(EoS_T *const eos)
{
  Physics_T *const phys = eos->phys;
  double *K,*rho0_th,*gamma;/* quantities in polytropic EoS */
  Uint N;/* number of pieces in case of pwp */
  Uint i;
  
  /* populate eos struct */
  sprintf(eos->description,"%s",Gets(P_"description"));
  sprintf(eos->type,"%s",       Gets(P_"type")); 
  sprintf(eos->unit,"%s",       Gets(P_"unit"));
  
  /* NOTE: order matters */
  if (strcmp_i(eos->type,"polytropic") || strcmp_i(eos->type,"piecewise_polytropic") || 
        strcmp_i(eos->type,"polytrop") || strcmp_i(eos->type, "pwp") || 
        strcmp_i(eos->type,"pwp_natural_cubic_spline"))
  {
    gamma  = read_EoS_in_parameter_file(Gets(P_"Gamma"),&N);
    K      = read_EoS_in_parameter_file(Gets(P_"K0"),0);/* this is K0 */
  }
  
  /* if we have a single polytropic eos, then we don't have rho0_th  */
  if (!strcmp_i(eos->type,"polytropic") &&
      !strcmp_i(eos->type,"polytrop") &&
      !strcmp_i(eos->type,"tabular"))
    rho0_th = read_EoS_in_parameter_file(Gets(P_"rho0_th"),0);   
  else
    rho0_th = 0;
    
  /* if the units are geometrised units */
  if (strcmp_i(eos->unit,"geo"))
  {
    /* if it is pwp */
    if (strcmp_i(eos->type,"piecewise_polytropic") ||
        strcmp_i(eos->type,"pwp"))
    {    
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
    }
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
      eos->cubic_spline->sample_size = sample_s;
      eos->cubic_spline->h_sample    = h_sample;
      eos->cubic_spline->p_sample    = p_sample;
      eos->cubic_spline->e_sample    = e_sample;
      eos->cubic_spline->rho0_sample = rho0_sample;
      eos->cubic_spline->h_floor     = Getd(P_"enthalpy_floor");
      
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
      eos->cubic_spline->interp_p           = interp_p;
      
      // e:
      Interpolation_T *interp_e             = init_interpolation();
      interp_e->method                      = "Natural_Cubic_Spline_1D";
      interp_e->N_cubic_spline_1d->f        = e_sample;
      interp_e->N_cubic_spline_1d->x        = h_sample;
      interp_e->N_cubic_spline_1d->N        = sample_s;
      interp_e->N_cubic_spline_1d->No_Warn  = 1;/* suppress warning */
      plan_interpolation(interp_e);
      eos->cubic_spline->interp_e           = interp_e;
      
      // rho0:
      Interpolation_T *interp_rho0              = init_interpolation();
      interp_rho0->method                       = "Natural_Cubic_Spline_1D";
      interp_rho0->N_cubic_spline_1d->f         = rho0_sample;
      interp_rho0->N_cubic_spline_1d->x         = h_sample;
      interp_rho0->N_cubic_spline_1d->N         = sample_s;
      interp_rho0->N_cubic_spline_1d->No_Warn   = 1;/* suppress warning */
      plan_interpolation(interp_rho0);
      eos->cubic_spline->interp_rho0            = interp_rho0;
      
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
    }
    else if (strcmp_i(eos->type,"polytropic") ||
             strcmp_i(eos->type,"polytrop"))
    {
      if (N != 1)
        Error0("This EoS is not polytropic, there is more than one piece.\n");
        
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
    }
    
    ////////////////////////////////////////////////////////////////////////////////Tabular EOS
    //Generates tabular EOS
    else if (strcmp_i(eos->type, "tabular") || strcmp_i(eos->type, "tab"))
    {
        FILE* eos_table = fopen(Pgets("eos_table_name"),"r");        //Name of EOS table
        if (!eos_table) { Error0("ERROR: Could not open EOS table."); }
        
        //Reads number of data points in EOS file.
        Uint sample_s = get_sample_size(Pgets("eos_table_name"));    //Defined in eos_tabular.c
        if (strstr_i(Pgets("EOS_table_format"), "Lorene")) { sample_s -= 9; }
        
        eos->cubic_spline->sample_size = sample_s;
        double *h_sample    = alloc_double(sample_s); //Arrays for tabular data points
        double *p_sample    = alloc_double(sample_s);
        double *e_sample    = alloc_double(sample_s);
        double *rho0_sample = alloc_double(sample_s);
        
        // Sets flag for log-log interpolation based on parameter file.
        if (strstr_i(PgetsEZ("EOS_log_approach"), "yes"))
        { eos->cubic_spline->use_log_approach = 1; }
        else { eos->cubic_spline->use_log_approach = 0; }
        
        //if (eos->cubic_spline->use_log_approach)
        //if (1) ///////////////FIXME: Won't compile unless arrays allocated
        //{
          double* h_log    = alloc_double(sample_s);
          double* p_log    = alloc_double(sample_s);
          double* e_log    = alloc_double(sample_s);
          double* rho0_log = alloc_double(sample_s);
        //}
        
        //Reads EOS data from text file.
        double h_point;
        double p_point;
        double rho0_point;
        double e_point;
        
        if (!PgetsEZ("EOS_table_format"))
        { Error0("ERROR: EOS table format not specified.\n"); }
        else if (strstr_i(PgetsEZ("EOS_table_format"), "Elliptica"))
        {
            // Default format in columns [pressure] [rest-mass density] [energy density] [enthalpy]
            // in geometrized (G = c = solar mass = 1) units.
            for (unsigned int line=0; line<sample_s; line++)
            {
                if (fscanf(eos_table, "%lf %lf %lf %lf\n", &p_point, &rho0_point, &e_point, &h_point) != 4)
                {
                    Error0("ERROR reading EOS data table.");
                    return;
                }

                p_sample[line]      = p_point;
                rho0_sample[line]   = rho0_point;
                e_sample[line]      = e_point;
                h_sample[line]      = h_point;
                
                // If we're using log-log interpolation, fill arrays with log(values).
                if (eos->cubic_spline->use_log_approach)
                {
                  p_log[line]      = log(p_point);
                  rho0_log[line]   = log(rho0_point);
                  e_log[line]      = log(e_point);
                  h_log[line]      = log(h_point);
                }
            }
        }
        else if (strstr_i(PgetsEZ("EOS_table_format"), "Lorene"))
        {
            // Lorene format in columns
            //          [line number] [number density] [(total) energy density] [pressure].
            // in units               [1/fm^3]         [g/cm^3]                 [dyn/cm^2].
            // Calculates rest-mass density (total) energy density, 
            // and specific enthalpy, and converts to geometrized units.
            
            // Physical constants from 2018 CODATA values.
            //double G_const = 6.67430E-11;          // Gravitational constant in m^3/(kg s^2)
            //double c_const = 299792458;            // Speed of light in m/s
            //double M_const = 1.98841E30;           // Solar mass in kg
            //double mn_const = 1.67492749804E-27;   // Neutron mass in kg
            // Lorene data is converted to SI units then geometrized.
            
            double n_point;
            char dummy_var_1;
            Uint dummy_var;
            // Pressure conversion factor: G^3 * M_solar^2 / (10 * c^8)
            // Energy density conversion factor: G^3  * M_solar^2 / (10 * c^6)
            // Rest-mass conversion factor: 10^45 * neutron mass * G^3 * m_solar ^2 / (c^6)
            // (Rest-mass is calculated by multiplying the baryon density by the baryon mass.)
            const double p_FACTOR = 1.80162095578956E-39;
            const double rho0_FACTOR = 0.002712069678583313;
            const double e_FACTOR = 1.619216164136643E-18;
               
            Uint dummy_var_2;
            for (Uint ctr=0; ctr<5; ctr++)
            {
                dummy_var_2 = (Uint)fscanf(eos_table, "%c\n", &dummy_var_1);
                dummy_var_2 += 1;
            }
            dummy_var_2 = (Uint)fscanf(eos_table, "%u\n", &dummy_var);
            for (Uint ctr=0; ctr<3; ctr++)
            {
                dummy_var_2 = (Uint)fscanf(eos_table, "%c\n", &dummy_var_1);
                dummy_var_2 += 1;
            }
            for (Uint line=0; line<sample_s; line++)
            {
                if (fscanf(eos_table, "%u %lf %lf %lf\n", &dummy_var, &n_point, &e_point, &p_point) != 4)
                {
                    Error0("ERROR reading EOS data table.");
                    return;
                }
                
                p_sample[line]      = p_point * p_FACTOR;
                rho0_sample[line]   = n_point * rho0_FACTOR;
                e_sample[line]      = e_point * e_FACTOR;
                h_sample[line]      = (p_sample[line] + e_sample[line]) / rho0_sample[line];
                
                if (eos->cubic_spline->use_log_approach)
                {
                  p_log[line]      = log(p_point * p_FACTOR);
                  rho0_log[line]   = log(n_point * rho0_FACTOR);
                  e_log[line]      = log(e_point * e_FACTOR);
                  h_log[line]      = log((p_point * p_FACTOR + e_point * e_FACTOR) / (n_point * rho0_FACTOR));
                }
            }
        }
        else
        { Error0("ERROR: Unrecognized EOS table format.\n"); }
        
        fclose(eos_table);
        //Sets interpolation bounds.
        
        if (eos->cubic_spline->use_log_approach)
        {
          eos->cubic_spline->h_floor = exp(h_sample[0]);
          eos->cubic_spline->h_max = exp(h_sample[sample_s-1]);
        }
        else
        {
          eos->cubic_spline->h_floor = h_sample[0];
          eos->cubic_spline->h_max = h_sample[sample_s-1];
        }
        
        //Saves data points in EOS.
        eos->cubic_spline->sample_size = sample_s;
        eos->cubic_spline->h_sample    = h_sample;
        eos->cubic_spline->p_sample    = p_sample;
        eos->cubic_spline->e_sample    = e_sample;
        eos->cubic_spline->rho0_sample = rho0_sample;
        eos->cubic_spline->h_floor     = Getd(P_"enthalpy_floor");
        if (eos->cubic_spline->use_log_approach)
        {
          eos->cubic_spline->h_log    = h_log;
          eos->cubic_spline->p_log    = p_log;
          eos->cubic_spline->e_log    = e_log;
          eos->cubic_spline->rho0_log = rho0_log;
        }
        
        Interpolation_T *interp_p    = init_interpolation(); //pressure
        Interpolation_T *interp_e    = init_interpolation(); //energy density
        Interpolation_T *interp_rho0 = init_interpolation(); //rest-mass density


        interp_p->method    = Gets(P_"Interpolation_Method");
        interp_e->method    = Gets(P_"Interpolation_Method");
        interp_rho0->method = Gets(P_"Interpolation_Method");
        assign_interpolation_ptrs(interp_p);
        assign_interpolation_ptrs(interp_e);
        assign_interpolation_ptrs(interp_rho0);
        if (strstr_i(Gets(P_"Interpolation_Method"), "Hermite") || 
            strstr_i(Gets(P_"Interpolation_Method"), "Clamped_Cubic_Spline"))
        {
            interp_p->finite_diff_order    = (Uint)Pgeti("interpolation_finite_diff_order");
            interp_e->finite_diff_order    = (Uint)Pgeti("interpolation_finite_diff_order");
            interp_rho0->finite_diff_order = (Uint)Pgeti("interpolation_finite_diff_order");
        }
        
        if (eos->cubic_spline->use_log_approach)
        {
          *interp_p->f    = p_log;
          *interp_p->x    = h_log;
          *interp_e->f    = e_log;
          *interp_e->x    = h_log;
          *interp_rho0->f = rho0_log;
          *interp_rho0->x = h_log;
        }
        else
        {
          *interp_p->f    = p_sample;
          *interp_p->x    = h_sample;
          *interp_e->f    = e_sample;
          *interp_e->x    = h_sample;
          *interp_rho0->f = rho0_sample;
          *interp_rho0->x = h_sample;
        }
        
        *interp_p->N    = sample_s;
        *interp_e->N    = sample_s;
        *interp_rho0->N = sample_s;
        set_interp_warn_flag(interp_p, 1); // suppress warning
        set_interp_warn_flag(interp_e, 1);
        set_interp_warn_flag(interp_rho0, 1);
        
        plan_interpolation(interp_p);
        plan_interpolation(interp_e);
        plan_interpolation(interp_rho0);
        eos->cubic_spline->interp_p = interp_p;
        eos->cubic_spline->interp_e = interp_e;
        eos->cubic_spline->interp_rho0 = interp_rho0;
        
        // assign functions for (p, e, rho0)
        eos->pressure                 = EoS_p_h_tab;
        eos->energy_density           = EoS_e_h_tab;
        eos->rest_mass_density        = EoS_rho0_h_tab;
        eos->specific_internal_energy = EoS_e0_h_tab;
        eos->de_dh                    = EoS_de_dh_h_tab;
        eos->drho0_dh                 = EoS_drho0_dh_h_tab;
        
        /////////////////////////Root finder approach to EOS//////////////////
        // The general approach is to invert rho0(h) using the equation
        // (e(rho0) + p(rho0)) / rho0 - h = 0,
        // and use rho0 as the independent variable by interpolating
        // e(rho0) and p(rho0).
        if (strstr_i(PgetsEZ(P_"Approach"), "Root_Finder"))
        {
          // Set up interpolation for p(rho0) and e(rho0)
          Interpolation_T* interp_p_rho0 = init_interpolation(); // Pressure from rest-mass
          Interpolation_T* interp_e_rho0 = init_interpolation(); // Total energy from rest-mass
          interp_p_rho0->method = Gets(P_"Interpolation_Method");
          interp_e_rho0->method = Gets(P_"Interpolation_Method");
          assign_interpolation_ptrs(interp_p_rho0);
          assign_interpolation_ptrs(interp_e_rho0);
          *interp_p_rho0->x = rho0_sample;
          *interp_e_rho0->x = rho0_sample;
          *interp_p_rho0->f = p_sample;
          *interp_e_rho0->f = e_sample;
          *interp_p_rho0->N = sample_s;
          *interp_e_rho0->N = sample_s;
          set_interp_warn_flag(interp_p_rho0, 1); // suppress warnings
          set_interp_warn_flag(interp_e_rho0, 1);
          plan_interpolation(interp_p_rho0);
          plan_interpolation(interp_e_rho0);
          eos->cubic_spline->interp_p_rho0    = interp_p_rho0;
          eos->cubic_spline->interp_e_rho0    = interp_e_rho0;
          eos->cubic_spline->enthalpy_eqn[0]  = EoS_enthalpy_def;
          
          // Set up root finder for (e(rho0) + p(rho0)) / rho0 - h = 0
          Root_Finder_T* root_finder  = init_root_finder(1);
          root_finder->type           = "Bisect_Single";
          root_finder->f[0]           = *eos->cubic_spline->enthalpy_eqn;
          root_finder->params         = eos;
          root_finder->tolerance      = 1E-8;  //FIXME: Make dynamic pars
          root_finder->MaxIter        = 100;
          // FIXME: Errors near h ~ 1.0, temporary solution /////////////
          root_finder->a_bisect       = rho0_sample[1];///////////////////////////
          /////////////////////////////////////////////////////
          root_finder->b_bisect       = rho0_sample[sample_s-1];
          root_finder->verbose = 0;
          root_finder->x_sol[0]       = eos->cubic_spline->rho0;
          plan_root_finder(root_finder);
          
          eos->cubic_spline->root_finder    = root_finder;
          eos->rest_mass_density            = EoS_rho0_RF;
          eos->pressure                     = EoS_p_rho0_tab;
          eos->energy_density               = EoS_e_rho0_tab;
          eos->specific_internal_energy     = EoS_e0_rho0_tab;
          eos->de_dh                        = EoS_de_dh_RF;
          eos->drho0_dh                     = EoS_drho0_dh_RF;
        }
        
        h_sample = 0;
        p_sample = 0;
        e_sample = 0;
        rho0_sample = 0;
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    
    else
    {
      Error0(NO_JOB);
    }
    
  }
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
  
  for(i = 0; i < eos->N; ++i) { printf("K%i == %E\n", i, K[i]); }
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
