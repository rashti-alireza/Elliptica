/*
// Alireza Rashti
// August 2019
*/

#include "bbn_study_bbn.h"

/* study and analyze initial data. whenever it is called, 
// it increments the solving_iteration_number by 1. */
void bbn_study_initial_data(Grid_T *const grid)
{
  if (Pgeti("STOP"))
    return;
    
  pr_clock();
  pr_line_custom('=');
  printf("{ Studying initial data for binary BH and NS ...\n");
  
  if (!grid)
  {
    printf("~> grid is empty.\n");
    printf("} Studying initial data for binary BH and NS ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
    return;
  }
  
  const char *const folder = Pgets("Diagnostics");
  int solving_iter         = PgetiEZ("solving_iteration_number");
  
  /* calculating the constraints */
  bbn_calculate_constraints_1st(grid);
  bbn_calculate_constraints_2nd(grid);
  
  /* make sure solving_iter is set */
  if (solving_iter != INT_MAX)
  {
    /* calculating ADM , Kommar masses ratios, spins, errors etc. */
    if (Pcmps("Elliptic_Convergence_Test","no"))
      bbn_measures(grid);
    
    /* prints */
    bbn_print_fields(grid,(unsigned)solving_iter,folder);
    bbn_print_residual_norms(grid,(unsigned)solving_iter,folder);
    bbn_print_properties(grid,(unsigned)solving_iter,folder,"a",1);
    Pseti("solving_iteration_number",++solving_iter);
  }
  
  printf("} Studying initial data for binary BH and NS ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* calculating ADM and Kommar masses, error, erc. */
void bbn_measures(Grid_T *const grid)
{
  Observable_T *obs = 0;
  double adm_mass,kommar_mass;
  const int solving_iter   = PgetiEZ("solving_iteration_number");
  const double ns_adm_mass = Pgetd("NS_ADM_mass");
  const double bh_irr_mass = Pgetd("BH_irreducible_mass");
  double virial_error, binding_energy,mass_ratio;
  const double NS_R_avg = NS_r_average(grid);
  double Rc_NS[3] = {0}, Rc_BH[3] = {0};/* centers */
  double S_NS[3]  = {0}, S_BH[3]  = {0};/* spins */
  double S_NS_JRP_method[3]  = {0}, S_BH_JRP_method[3]  = {0};/* spins */
  
  double chris_bh_mass = 1;/* Christodoulou mass, 1 to avoid division by 0 */
  
  /* adm mass */
  obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "ADM(M)|BBN";
  plan_observable(obs);
  adm_mass = obs->M(obs);
  free_observable(obs);
 
  /* kommar mass */
  obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "Kommar(M)|BBN";
  plan_observable(obs);
  kommar_mass = obs->M(obs);
  free_observable(obs);
  
  /* mass ratio, E_b , error */
  mass_ratio     = bh_irr_mass/ns_adm_mass;
  binding_energy = adm_mass-bh_irr_mass-ns_adm_mass;
  virial_error   = fabs(1-kommar_mass/adm_mass);
  
  /* some pars are defined after 0 iter */
  if (solving_iter > 0 && solving_iter != INT_MAX)
  {
    /* NS spin */
    bbn_define_spin_JRP(S_NS_JRP_method,grid,"NS");
    bbn_define_spin_integral(S_NS,grid,"NS");
    /* NS center */
    bbn_Rc_NS(Rc_NS,grid);
    
    /* BH spin */
    bbn_define_spin_JRP(S_BH_JRP_method,grid,"BH");
    bbn_define_spin_integral(S_BH,grid,"BH");
    /* BH center */
    bbn_Rc_BH(Rc_BH,grid);
    
    double s_BH2  = Pow2(S_BH[0])+Pow2(S_BH[1])+Pow2(S_BH[2]);
    double irrm2  = Pow2(bh_irr_mass);
    chris_bh_mass = sqrt(irrm2+s_BH2/(4*irrm2));
    Psetd("Christodoulou_mass",chris_bh_mass);
  }
  
  /* update parameters */ 
  Psetd("BBN_ADM_mass"   ,adm_mass);
  Psetd("BBN_Kommar_mass",kommar_mass);
  Psetd("Binding_energy" ,binding_energy);
  Psetd("Virial_error"   ,virial_error);
  Psetd("mass_ratio"     ,mass_ratio);
  Psetd("NS_compactness" ,Pgetd("NS_ADM_mass")/NS_R_avg);
  Psetd("NS_average_proper_radius" ,NS_R_avg);
  
  Psetd("NS_Sx_Campanelli",S_NS[0]);
  Psetd("NS_Sy_Campanelli",S_NS[1]);
  Psetd("NS_Sz_Campanelli",S_NS[2]);
  
  Psetd("NS_Sx_JRP",S_NS_JRP_method[0]);
  Psetd("NS_Sy_JRP",S_NS_JRP_method[1]);
  Psetd("NS_Sz_JRP",S_NS_JRP_method[2]);
  
  Psetd("NS_chi_x",S_NS[0]/Pow2(ns_adm_mass));
  Psetd("NS_chi_y",S_NS[1]/Pow2(ns_adm_mass));
  Psetd("NS_chi_z",S_NS[2]/Pow2(ns_adm_mass));
  
  Psetd("NS_Rcenter_x",Rc_NS[0]);
  Psetd("NS_Rcenter_y",Rc_NS[1]);
  Psetd("NS_Rcenter_z",Rc_NS[2]);
  
  Psetd("BH_Sx_Campanelli",S_BH[0]);
  Psetd("BH_Sy_Campanelli",S_BH[1]);
  Psetd("BH_Sz_Campanelli",S_BH[2]);
  
  Psetd("BH_Sx_JRP",S_BH_JRP_method[0]);
  Psetd("BH_Sy_JRP",S_BH_JRP_method[1]);
  Psetd("BH_Sz_JRP",S_BH_JRP_method[2]);
  
  Psetd("BH_chi_x",S_BH[0]/Pow2(chris_bh_mass));
  Psetd("BH_chi_y",S_BH[1]/Pow2(chris_bh_mass));
  Psetd("BH_chi_z",S_BH[2]/Pow2(chris_bh_mass));
  
  Psetd("BH_Rcenter_x",Rc_BH[0]);
  Psetd("BH_Rcenter_y",Rc_BH[1]);
  Psetd("BH_Rcenter_z",Rc_BH[2]);
  
}

/* print the properites of the system for instance:
// mass, spin, momentum, distance and etc.
// folder        : output directory 
// open_file_mode: "a" for append and "w" for newfile
// pr_flg        : if pr_flg = 1 it also prints them at stout otherwise none. */
void bbn_print_properties(Grid_T *const grid,const unsigned iteration, const char *const folder,const char *const open_file_mode,const int pr_flg)
{
  const char *const file_name = "bbn_geometry_and_physics.txt";
  FILE *file = 0;
  char str[MAX_STR_LEN];
  const char *const geometry_logo = "\n############\n"
                                    "# Geometry #\n"
                                    "############\n";
  const char *const physics_logo = "\n###########\n"
                                   "# Physics #\n"
                                   "###########\n";
  
  /* open file */
  sprintf(str,"%s/%s",folder,file_name);
  file = Fopen(str,open_file_mode);
  fprintf(file,"%s\n",LINE_STR);
  fprintf(file,"# iteration = %u\n",iteration);
  

  if (pr_flg)
  {
    pr_line_custom('=');
    printf("{ Geometry and physics of binary BH and NS ...\n");
  }
  
  /* { geometry */
  if (pr_flg)
    printf(geometry_logo);
    
  fprintf(file,geometry_logo);
  
  /* these centers needed for ID reader */
  PR_PARAMETR_IN_FILE(NS_center_x)
  PR_PARAMETR_IN_FILE(NS_center_y)
  PR_PARAMETR_IN_FILE(NS_center_z)
  
  PR_PARAMETR_IN_FILE(BH_center_x)
  PR_PARAMETR_IN_FILE(BH_center_y)
  PR_PARAMETR_IN_FILE(BH_center_z)
  
  PR_PARAMETR_IN_FILE(r_excision)
  PR_PARAMETR_IN_FILE(NS_max_radius)
  PR_PARAMETR_IN_FILE(NS_min_radius)
  PR_PARAMETR_IN_FILE(BH_NS_separation)
  /* } geometry */
  
  /* { physics */
  if (pr_flg)
    printf(physics_logo);
  
  fprintf(file,physics_logo);
  
  PR_PARAMETR_IN_FILE_s(EoS_description)
  PR_PARAMETR_IN_FILE_s(EoS_type)
  PR_PARAMETR_IN_FILE_s(EoS_unit)
  PR_PARAMETR_IN_FILE_s(EoS_K)
  PR_PARAMETR_IN_FILE_s(EoS_rho_th)
  PR_PARAMETR_IN_FILE_s(EoS_Gamma)
  
  PR_PARAMETR_IN_FILE(BBN_ADM_mass)
  PR_PARAMETR_IN_FILE(BBN_Kommar_mass)
  PR_PARAMETR_IN_FILE(Binding_energy)
  
  PR_PARAMETR_IN_FILE(NS_baryonic_mass)
  PR_PARAMETR_IN_FILE(NS_baryonic_mass_current)
  PR_PARAMETR_IN_FILE(NS_ADM_mass)
  PR_PARAMETR_IN_FILE(NS_Kommar_mass)
  PR_PARAMETR_IN_FILE(NS_average_proper_radius)
  PR_PARAMETR_IN_FILE(NS_Sx_Campanelli)
  PR_PARAMETR_IN_FILE(NS_Sy_Campanelli)
  PR_PARAMETR_IN_FILE(NS_Sz_Campanelli)
  PR_PARAMETR_IN_FILE(NS_Sx_JRP)
  PR_PARAMETR_IN_FILE(NS_Sy_JRP)
  PR_PARAMETR_IN_FILE(NS_Sz_JRP)
  
  PR_PARAMETR_IN_FILE(NS_chi_x)
  PR_PARAMETR_IN_FILE(NS_chi_y)
  PR_PARAMETR_IN_FILE(NS_chi_z)
  PR_PARAMETR_IN_FILE(NS_Rcenter_x)
  PR_PARAMETR_IN_FILE(NS_Rcenter_y)
  PR_PARAMETR_IN_FILE(NS_Rcenter_z)
  
  PR_PARAMETR_IN_FILE(BH_irreducible_mass)
  PR_PARAMETR_IN_FILE(BH_irreducible_mass_current)
  PR_PARAMETR_IN_FILE(BH_ADM_mass)
  PR_PARAMETR_IN_FILE(BH_Kommar_mass)
  PR_PARAMETR_IN_FILE(Christodoulou_mass)
  
  PR_PARAMETR_IN_FILE(BH_AH_area)
  PR_PARAMETR_IN_FILE(BH_Sx_Campanelli)
  PR_PARAMETR_IN_FILE(BH_Sy_Campanelli)
  PR_PARAMETR_IN_FILE(BH_Sz_Campanelli)
  PR_PARAMETR_IN_FILE(BH_Sx_JRP)
  PR_PARAMETR_IN_FILE(BH_Sy_JRP)
  PR_PARAMETR_IN_FILE(BH_Sz_JRP)
  PR_PARAMETR_IN_FILE(BH_chi_x)
  PR_PARAMETR_IN_FILE(BH_chi_y)
  PR_PARAMETR_IN_FILE(BH_chi_z)
  PR_PARAMETR_IN_FILE(BH_Omega_U0)
  PR_PARAMETR_IN_FILE(BH_Omega_U1)
  PR_PARAMETR_IN_FILE(BH_Omega_U2)
  PR_PARAMETR_IN_FILE(BH_Rcenter_x)
  PR_PARAMETR_IN_FILE(BH_Rcenter_y)
  PR_PARAMETR_IN_FILE(BH_Rcenter_z)
    
  PR_PARAMETR_IN_FILE(mass_ratio)
  PR_PARAMETR_IN_FILE(NS_compactness)
  
  PR_PARAMETR_IN_FILE(Px_ADM)
  PR_PARAMETR_IN_FILE(Py_ADM)
  PR_PARAMETR_IN_FILE(Pz_ADM)
  PR_PARAMETR_IN_FILE(Jx_ADM)
  PR_PARAMETR_IN_FILE(Jy_ADM)
  PR_PARAMETR_IN_FILE(Jz_ADM)
  
  PR_PARAMETR_IN_FILE(NS_Px_ADM)
  PR_PARAMETR_IN_FILE(NS_Py_ADM)
  PR_PARAMETR_IN_FILE(NS_Pz_ADM)
  PR_PARAMETR_IN_FILE(NS_Jx_ADM)
  PR_PARAMETR_IN_FILE(NS_Jy_ADM)
  PR_PARAMETR_IN_FILE(NS_Jz_ADM)
  
  PR_PARAMETR_IN_FILE(BH_Px_ADM)
  PR_PARAMETR_IN_FILE(BH_Py_ADM)
  PR_PARAMETR_IN_FILE(BH_Pz_ADM)
  PR_PARAMETR_IN_FILE(BH_Jx_ADM)
  PR_PARAMETR_IN_FILE(BH_Jy_ADM)
  PR_PARAMETR_IN_FILE(BH_Jz_ADM)
  
  PR_PARAMETR_IN_FILE(BH_NS_angular_velocity)
  PR_PARAMETR_IN_FILE(BH_NS_infall_velocity)
  PR_PARAMETR_IN_FILE(x_CM)
  PR_PARAMETR_IN_FILE(y_CM)
  PR_PARAMETR_IN_FILE(z_CM)
  
  PR_PARAMETR_IN_FILE(rho_center)
  PR_PARAMETR_IN_FILE(pressure_center)
  PR_PARAMETR_IN_FILE(energy_density_center)
  PR_PARAMETR_IN_FILE(Euler_equation_constant)
  
  PR_PARAMETR_IN_FILE(Virial_error)
  PR_PARAMETR_IN_FILE(largest_L2norm_error)
   
  /* } physics */
  
  fprintf(file,"%s\n",LINE_STR);
  fclose(file);
  
  if (pr_flg)
  {
    printf("} Geometry and physics of binary BH and NS ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
  }
  
  UNUSED(grid);
}

/* print residual norms L2, L1 and L_inf of the specified fields. */
void bbn_print_residual_norms(Grid_T *const grid,const unsigned iteration, const char *const folder)
{
  /* list of the fields to be printed out */
  char **f = 
     read_separated_items_in_string(PgetsEZ("output_2d_txt"),',');
  double largest_L2_error = 0;                   
  unsigned i,p;
  
  if(f)         
  for (i = 0; f[i]; ++i)
  {
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      unsigned nn    = patch->nn;
      Field_T *field;
      int field_ind = _Ind(f[i]);
      double Linf,L2,L1;
      FILE *file_Linf,*file_L1,*file_L2;
      char file_name_Linf[1000];
      char file_name_L1[1000];
      char file_name_L2[1000];
      char *stem = strstr(patch->name,"_");
      stem++;
      
      if (field_ind < 0)
        continue;
      
      field = patch->pool[field_ind];
      if (!field->v)/* if field empty */
        continue;
      
      sprintf(file_name_Linf,"%s/%s_Linf_%s.table",folder,f[i],stem);
      sprintf(file_name_L1,  "%s/%s_L1_%s.table",folder,f[i],stem);
      sprintf(file_name_L2,  "%s/%s_L2_%s.table",folder,f[i],stem);
      
      if (access(file_name_Linf,F_OK) != -1)/* if file exists */
      {
        file_Linf = Fopen(file_name_Linf,"a");
      }
      else
      {
        file_Linf = Fopen(file_name_Linf,"w");
        fprintf(file_Linf,"#iteration  %s\n",f[i]);
      }
      
      if (access(file_name_L1,F_OK) != -1)
      {
        file_L1 = Fopen(file_name_L1,"a");
      }
      else
      {
        file_L1 = Fopen(file_name_L1,"w");
        fprintf(file_L1,"#iteration  %s\n",f[i]);
      }
      
      if (access(file_name_L2,F_OK) != -1)
      {
        file_L2 = Fopen(file_name_L2,"a");
      }
      else
      {
        file_L2 = Fopen(file_name_L2,"w");
        fprintf(file_L2,"#iteration  %s\n",f[i]);
      }
        
      Linf  = L_inf(nn,field->v);
      L2    = L2_norm(nn,field->v,0);
      L1    = L1_norm(nn,field->v,0);
      
      fprintf(file_Linf,"%-11u %0.15f\n",iteration,Linf);
      fprintf(file_L1,  "%-11u %0.15f\n",iteration,L1);
      fprintf(file_L2,  "%-11u %0.15f\n",iteration,L2);
      
      /* since we do all sort of approx. inside the BH don't 
      // consider this as the numerical error. */
      if (!IsItInsideBHPatch(patch))
        largest_L2_error = L2 > largest_L2_error ? L2 : largest_L2_error;
      
      fclose(file_Linf);
      fclose(file_L1);
      fclose(file_L2);
    }
  }
  free_2d(f);
  /* update numeric error */
  Psetd("largest_L2norm_error",largest_L2_error);
}

/* printing fields determined in parameter file */
void bbn_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder)
{
  pr_line_custom('=');
  printf("{ Printing Specified Fields for Binary BH and NS ...\n");

  Pr_Field_T *pr  = init_PrField(grid);
  pr->folder = folder;
  pr->cycle  = (int)iteration;
  pr_fields(pr);
  free_PrField(pr);
  
  printf("} Printing Specified Fields for Binary BH and NS ==> Done.\n");
  pr_line_custom('=');
}

/* calculate the proper area of the NS then using area = 4 pi R^2 to find R
// as the avarage of NS radius.
// ->return value: avarage NS radius */
static double NS_r_average(Grid_T *const grid)
{
  double R = 0, area = 0;
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItNSSurroundingPatch(patch))
      continue;
      
    unsigned ijk;
    unsigned nn = patch->nn;
    ADD_FIELD(NS_R_average_integrand)
    READ_v(_gamma_D2D2)
    READ_v(_gamma_D0D2)
    READ_v(_gamma_D0D0)
    READ_v(_gamma_D0D1)
    READ_v(_gamma_D1D2)
    READ_v(_gamma_D1D1)
    READ_v(psi);
    double *g00 = alloc_double(nn);
    double *g01 = alloc_double(nn);
    double *g02 = alloc_double(nn);
    double *g11 = alloc_double(nn);
    double *g12 = alloc_double(nn);
    double *g22 = alloc_double(nn);
    
    {/* local variables */
      REALLOC_v_WRITE_v(NS_R_average_integrand)

      FOR_ALL_POINTS(ijk,patch)
      {
        NS_R_average_integrand[ijk] = 1;
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*_gamma_D0D0[ijk];
        g01[ijk] = psi4*_gamma_D0D1[ijk];
        g02[ijk] = psi4*_gamma_D0D2[ijk];
        g11[ijk] = psi4*_gamma_D1D1[ijk];
        g12[ijk] = psi4*_gamma_D1D2[ijk];
        g22[ijk] = psi4*_gamma_D2D2[ijk];
      }
    }
    DECLARE_FIELD(NS_R_average_integrand)
    Integration_T *I = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->Spectral->f = NS_R_average_integrand;
    I->g00 = g00;
    I->g01 = g01;
    I->g02 = g02;
    I->g11 = g11;
    I->g12 = g12;
    I->g22 = g22;
    I->Spectral->Z_surface = 1;
    I->Spectral->K         = 0;
    plan_integration(I);
    area += execute_integration(I);
    free_integration(I);
    REMOVE_FIELD(NS_R_average_integrand)
    free(g00);
    free(g01);
    free(g02);
    free(g11);
    free(g12);
    free(g22);
  }
  R = sqrt(area/(4*M_PI));
  return R;
}

/* calculating physical center of BH to be used in spin calculations */
void bbn_Rc_BH(double Rc[3],Grid_T *const grid)
{
  const double AH_area = Pgetd("BH_AH_area");
  const double x_CM = Pgetd("x_CM");
  const double y_CM = Pgetd("y_CM");
  unsigned p;

  Rc[0] = 0;
  Rc[1] = 0;
  Rc[2] = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if (!IsItHorizonPatch(patch))
      continue;
      
    unsigned ijk;
    unsigned nn = patch->nn;
    
    READ_v(_gamma_D2D2)
    READ_v(_gamma_D0D2)
    READ_v(_gamma_D0D0)
    READ_v(_gamma_D0D1)
    READ_v(_gamma_D1D2)
    READ_v(_gamma_D1D1)
    READ_v(psi);
    ADD_FIELD(Rc_integrandx)
    ADD_FIELD(Rc_integrandy)
    ADD_FIELD(Rc_integrandz)
    
    double *g00 = alloc_double(nn);
    double *g01 = alloc_double(nn);
    double *g02 = alloc_double(nn);
    double *g11 = alloc_double(nn);
    double *g12 = alloc_double(nn);
    double *g22 = alloc_double(nn);
    
    {/* local variables */
      REALLOC_v_WRITE_v(Rc_integrandx)
      REALLOC_v_WRITE_v(Rc_integrandy)
      REALLOC_v_WRITE_v(Rc_integrandz)

      FOR_ALL_POINTS(ijk,patch)
      {
        double x = patch->node[ijk]->x[0];
        double y = patch->node[ijk]->x[1];
        double z = patch->node[ijk]->x[2];
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*_gamma_D0D0[ijk];
        g01[ijk] = psi4*_gamma_D0D1[ijk];
        g02[ijk] = psi4*_gamma_D0D2[ijk];
        g11[ijk] = psi4*_gamma_D1D1[ijk];
        g12[ijk] = psi4*_gamma_D1D2[ijk];
        g22[ijk] = psi4*_gamma_D2D2[ijk];
        Rc_integrandx[ijk] = x-x_CM;
        Rc_integrandy[ijk] = y-y_CM;
        Rc_integrandz[ijk] = z;
      }
    }
    DECLARE_FIELD(Rc_integrandx)
    DECLARE_FIELD(Rc_integrandy)
    DECLARE_FIELD(Rc_integrandz)
    Integration_T *I = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->g00 = g00;
    I->g01 = g01;
    I->g02 = g02;
    I->g11 = g11;
    I->g12 = g12;
    I->g22 = g22;
    I->Spectral->Z_surface = 1;
    I->Spectral->K         = 0;
    
    I->Spectral->f = Rc_integrandx;
    plan_integration(I);
    Rc[0] += execute_integration(I);
    
    I->Spectral->f = Rc_integrandy;
    plan_integration(I);
    Rc[1] += execute_integration(I);
    
    I->Spectral->f = Rc_integrandz;
    plan_integration(I);
    Rc[2] += execute_integration(I);
    
    free_integration(I);
    REMOVE_FIELD(Rc_integrandx)
    REMOVE_FIELD(Rc_integrandy)
    REMOVE_FIELD(Rc_integrandz)
    free(g00);
    free(g01);
    free(g02);
    free(g11);
    free(g12);
    free(g22);
  }
  Rc[0] /= (AH_area);
  Rc[1] /= (AH_area);
  Rc[2] /= (AH_area);
}


