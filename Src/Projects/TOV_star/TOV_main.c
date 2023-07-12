/*
// Alireza Rashti
// June 2019
*/

#include "TOV_main.h"
#define STR_LEN_MAX 900

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
  
  TOV_T *tov = TOV_init();
  tov->phys  = phys;
  const char *const path_par = Pgets("top_directory");
  char *path =  make_directory(path_par,"TOV_Star");
  char file_name[STR_LEN_MAX];
  FILE *file;
  const double *p,*e,*rho0, *eps,*m,*m0,*r,*rbar;
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
  
  /* print energy density */
  e = tov->e;
  sprintf(file_name,"%s/energy_density.2d",path);
  file = Fopen(file_name,"w+");
  fprintf(file,"#Schwarzchild_r  energy_density\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],e[i]);
  
  Fclose(file);
  
  /* print rest mass density */
  rho0 = tov->rho0;
  sprintf(file_name,"%s/rest_mass_density.2d",path);
  file = Fopen(file_name,"w+");
  fprintf(file,"#Schwarzchild_r  rest_mass_density\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],rho0[i]);
  
  Fclose(file);
  
  /* print specific internal energy */
  eps = tov->eps;
  sprintf(file_name,"%s/internal_energy.2d",path);
  file = Fopen(file_name,"w+");
  fprintf(file,"#Schwarzchild_r  specific internal energy\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],eps[i]);
  
  Fclose(file);
  
  /* print mass */
  m = tov->m;
  sprintf(file_name,"%s/mass.2d",path);
  file = Fopen(file_name,"w+");
  fprintf(file,"#Schwarzchild_r  mass\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],m[i]);
    
  Fclose(file);
  
  /* print rest mass */
  m0 = tov->m0;
  sprintf(file_name,"%s/rest_mass.2d",path);
  file = Fopen(file_name,"w+");
  fprintf(file,"#Schwarzchild_r  rest_mass\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],m0[i]);
    
  Fclose(file);
  
  /* print isotropic r */
  rbar = tov->rbar;
  sprintf(file_name,"%s/isotropic_radius.2d",path);
  file = Fopen(file_name,"w+");
  fprintf(file,"#Schwarzchild_r  isotropic_radius\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],rbar[i]);
    
  Fclose(file);
  
  //////////////////////////////Mass-Radius curve//////////////////
  // Option to generate mass radius curve. (ADM mass vs Schwarzschild radius
  // in solar masses and kilometers).
  if (strstr_i(PgetsEZ("mass_radius_curve"),"yes"))
  {
    double mass_initial = Pgetd("initial_mass");
    double mass_final   = Pgetd("final_mass");
    Uint stars          = (Uint)Pgeti("stars");
    if (mass_final <= mass_initial || stars == 0)
    { Error0("Error in mass-radius curve parameters.\n"); }
    
    double test_mass           = mass_initial;
    double delta_m             = (mass_final - mass_initial) / stars;
    double* radii              = alloc_double(stars);
    double* masses             = alloc_double(stars);
    double* central_enthalpies = alloc_double(stars);
    
    // Geometric units to km conversion factor: (G * Msolar / c^2) / (10^3)
    //double r_FACTOR = 1.47667;
    
    TOV_T* tov_star;
    for (Uint star = 0; star < stars; star++)
    {
      tov_star              = TOV_init();
      tov_star->description = 0; // Must be 'off' to avoid excessive print statements.
      tov_star->phys        = phys;
      tov_star->bar_m       = test_mass;
      tov_star              = TOV_solution(tov_star);
      radii[star]           = tov_star->r[tov->N-1];
      masses[star]          = tov_star->ADM_m;
      central_enthalpies[star] = tov_star->h[0];
      
      TOV_free(tov_star);
      test_mass += delta_m;
    }
    
    // Print results to file
    // mass vs radius
    sprintf(file_name,"%s/mass_radius.txt",path);
    file = Fopen(file_name,"w+");
    fprintf(file,"#Radius (km) \t ADM Mass\n");
    for (i = 0; i < stars; ++i)
      { fprintf(file,"%E \t %E\n",radii[i],masses[i]); }
    Fclose(file);
    
    // mass vs central enthalpy
    sprintf(file_name,"%s/mass_enthalpy.txt",path);
    file = Fopen(file_name,"w+");
    fprintf(file,"#Central Enthalpy \t ADM Mass\n");
    for (i = 0; i < stars; ++i)
      { fprintf(file,"%E \t %E\n",central_enthalpies[i],masses[i]); }
    Fclose(file);
    
    free(radii);
    free(masses);
  }
    
  
  free_physics(phys);
  TOV_free(tov);
  free(path);
  UNUSED(vp);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}
