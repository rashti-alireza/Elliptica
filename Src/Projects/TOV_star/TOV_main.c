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
    
    double test_mass = mass_initial;
    double delta_m = (mass_final - mass_initial) / stars;
    double* radii = alloc_double(stars);
    double* masses = alloc_double(stars);
    tov->description = 0; // Must be 'off' to avoid print statements.
    
    // Geometric units to km conversion factor: (G * Msolar / c^2) / (10^3)
    double r_FACTOR = 1.47667;
    
    for (Uint star = 0; star < stars; star++)
    {
      tov->bar_m = test_mass;
      tov = TOV_solution(tov);
      radii[star] = tov->r[tov->N-1] * r_FACTOR;
      masses[star] = tov->ADM_m;
      test_mass += delta_m;
    }
    
    // Print results to file
    sprintf(file_name,"%s/mass_radius.txt",path);
    file = Fopen(file_name,"w+");
    fprintf(file,"#Radius (km) \t ADM Mass (solar masses)\n");
    for (i = 0; i < stars; ++i)
      { fprintf(file,"%E \t %E\n",radii[i],masses[i]); }
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
