/*
// Alireza Rashti
// June 2019
*/

#include "TOV_main.h"
#define MAXSTR 400

/* solving the problem of TOV stars
// ->return value: if succeeds EXIT_SUCCESS */
int TOV_star(void)
{
  TOV_T *tov = TOV_init();
  const char *const path_par = Pgets_E("output_directory_path");
  char *path =  make_directory(path_par,"TOV_Star");
  char file_name[MAXSTR];
  FILE *file;
  const double *p,*m,*r,*rbar;
  unsigned N = (unsigned)GetParameterI_E("TOV_star_n");
  unsigned i;
  
  if (N %2 == 0)/* the method need odd number */
    N += 1;
    
  tov->N = N;
  tov->bar_m = GetParameterD_E("TOV_star_baryonic_mass");
  tov->description = "A TOV star";
  tov = TOV_solution(tov);
  r = tov->r;
  
  /* print pressure */
  p = tov->p;
  sprintf(file_name,"%s/pressure.2d",path);
  file = fopen(file_name,"w+");
  pointerEr(file);
  fprintf(file,"#Schwarzchild_r  pressure\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],p[i]);
  
  fclose(file);
  
  /* print mass */
  m = tov->m;
  sprintf(file_name,"%s/mass.2d",path);
  file = fopen(file_name,"w+");
  pointerEr(file);
  fprintf(file,"#Schwarzchild_r  mass\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],m[i]);
    
  fclose(file);
  
  /* print isotropic r */
  rbar = tov->rbar;
  sprintf(file_name,"%s/isotropic_radius.2d",path);
  file = fopen(file_name,"w+");
  pointerEr(file);
  fprintf(file,"#Schwarzchild_r  isotropic_radius\n");
  
  for (i = 0; i < N; ++i)
    fprintf(file,"%-15e  %e\n",r[i],rbar[i]);
    
  fclose(file);
  
  TOV_free(tov);
  free(path);
  
  return EXIT_SUCCESS;
}
