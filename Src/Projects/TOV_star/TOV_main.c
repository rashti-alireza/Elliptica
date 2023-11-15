/*
// Alireza Rashti
// June 2019
*/

#include "TOV_main.h"
#define STR_LEN_MAX (900)

/* solving the problem of TOV stars
// ->return value: if succeeds EXIT_SUCCESS */
int TOV_star(void *vp)
{
  FUNC_TIC

  Physics_T *const phys = init_physics(0,NS);

  /* making output directory for this project */
  char folder[STR_LEN_MAX] = {'\0'};
  char *outdir = 0;
  sprintf(folder,"%s",Pgets("parameter_file_name_stem"));
  outdir = make_directory(Pgets("relative_root_path"),folder);
  add_parameter("top_directory",outdir);
  free(outdir);

  /* set some default parameters: */
  
  /* number of points for composite Simpson's rule integral.
  // for accuracy we choose large number and the method
  // requires the number be an odd number. */
  Pset_default("TOV_Star_n","1001");
  
  if (strstr_i(PgetsEZ("TOV_star_mass_radius_curve"),"yes"))
  {
    multiple_tov(phys);
  }
  else // solve only one tov
  {
    single_tov(phys);
  }
  
  free_physics(phys);
  UNUSED(vp);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

// solving a single tov star for a given baryonic mass
static void single_tov(Physics_T *const phys)
{
  TOV_T *tov = TOV_init();
  tov->phys  = phys;
  const char *const path_par = Pgets("top_directory");
  char *path =  make_directory(path_par,"TOV_Star");
  char file_name[STR_LEN_MAX];
  FILE *file;
  const double *p,*m,*r,*rbar;
  const Uint N = (Uint)Pgeti("TOV_Star_n");
  Uint i;
  
  tov->bar_m = Pgetd("TOV_star_baryonic_mass");
  tov->description = "A TOV star";
  tov = TOV_solution(tov);
  r = tov->r;
  
  /* print pressure */
  p = tov->p;
  sprintf(file_name,"%s/pressure.2d",path);
  file = Fopen(file_name,"w+");
  fprintf(file,"#Schwarzchild_r  pressure\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],p[i]);
  
  Fclose(file);
  
  /* print mass */
  m = tov->m;
  sprintf(file_name,"%s/mass.2d",path);
  file = Fopen(file_name,"w+");
  fprintf(file,"#Schwarzchild_r  mass\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],m[i]);
    
  Fclose(file);
  
  /* print isotropic r */
  rbar = tov->rbar;
  sprintf(file_name,"%s/isotropic_radius.2d",path);
  file = Fopen(file_name,"w+");
  fprintf(file,"#Schwarzchild_r  isotropic_radius\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],rbar[i]);
    
  Fclose(file);
  free(path);
  TOV_free(tov);
}

// solving tov stars for a range of masses to find mass radius curve
static void multiple_tov(Physics_T *const phys)
{
  const double m_i = Pgetd("TOV_star_baryonic_mass_initial");
  const double m_f = Pgetd("TOV_star_baryonic_mass_final");
  const double dm  = Pgetd("TOV_star_baryonic_mass_step");
  const Uint n_tot  = (Uint)ceil((m_f - m_i)/dm) + 1;
  double *radii_sch = 0;// radius of each tov (Schwarzschild Coords.)
  double *radii_iso = 0;// radius of each tov (Isotropic Coords.)
  double *ms_bar    = 0;// baryonic mass of each tov
  double *ms_adm    = 0;// adm mass of each tov
  double *ps_r0    = 0;// pressure at r = 0
  double *rho0s_r0 = 0;// rho0 at r = 0
  double *es_r0    = 0;// e    at r = 0
  double *hs_r0    = 0;// enthalpy at r = 0
  const char *const path_par = Pgets("top_directory");
  char *path = make_directory(path_par,"TOV_Star");
  char file_name[STR_LEN_MAX];
  EoS_T *eos = init_EoS(phys);
  FILE *file;
  Uint i;
  
  radii_sch  = alloc_double(n_tot);
  radii_iso  = alloc_double(n_tot);
  ms_bar     = alloc_double(n_tot);
  ms_adm     = alloc_double(n_tot);
  ps_r0      = alloc_double(n_tot);
  rho0s_r0   = alloc_double(n_tot);
  es_r0      = alloc_double(n_tot);
  hs_r0      = alloc_double(n_tot);
  
  for (i = 0; i < n_tot; ++i)
  {
    double bar_m = m_i + i*dm;
    char des[STR_LEN_MAX];
    snprintf(des,100,"A TOV star with baryonic mass = %0.4e.",bar_m);
    ms_bar[i] = bar_m;
    
    TOV_T *tov = TOV_init();
    tov->phys  = phys;
    tov->bar_m = bar_m;
    tov->description = des;
    tov = TOV_solution(tov);
    radii_sch[i] = tov->r[tov->N-1];
    radii_iso[i] = tov->rbar[tov->N-1];
    ms_adm[i]    = tov->ADM_m;
    hs_r0[i]     = tov->h[0];
    eos->h       = tov->h[0];
    rho0s_r0[i]  = eos->rest_mass_density(eos);
    es_r0[i]     = eos->energy_density(eos);
    ps_r0[i]     = tov->p[0];
    TOV_free(tov);
  }
  
  // write results into a file
  sprintf(file_name,"%s/tov_properties.txt",path);
  file = Fopen(file_name,"w+");
  fprintf(file,"# 1:baryonic_mass 2:adm_mass 3:radius(Schwarzschild Coords.) "
          "4:radius(Isotropic Coords.) 5:pressure(r=0) 6:rest_mass_density(r=0) "
          "7:energy_density(r=0) 8:enthalpy(r=0) 9:compactness\n");
  for (i = 0; i < n_tot; ++i)
  {
    fprintf(file,"%0.15e",ms_bar[i]);
    fprintf(file," %0.15e",ms_adm[i]);
    fprintf(file," %0.15e",radii_sch[i]);
    fprintf(file," %0.15e",radii_iso[i]);
    fprintf(file," %0.15e",ps_r0[i]);
    fprintf(file," %0.15e",rho0s_r0[i]);
    fprintf(file," %0.15e",es_r0[i]);
    fprintf(file," %0.15e",hs_r0[i]);
    fprintf(file," %0.15e",ms_adm[i]/radii_sch[i]);
    fprintf(file,"\n");
  }
  Fclose(file);
  
  // cleanup
  Free(path);
  Free(radii_sch);
  Free(radii_iso);
  Free(ms_bar);
  Free(ms_adm);
  Free(ps_r0)
  Free(rho0s_r0)
  Free(es_r0)
  Free(hs_r0)
  free_EoS(eos);
}
