/*
// Alireza Rashti
// June 2019
*/

#include "eos_test.h"

void test_EoS(Physics_T *const phys)
{
  Grid_T *const grid = phys->grid;
  EoS_T *eos = init_EoS(phys);
  const char *const path_par = Pgets("top_directory");
  char *path,file_name[1000];
  FILE *file = 0;
  Uint N = 1000;
  const double h_max = 2.;// eos->h_th != 0 ? eos->h_th[eos->N-1]+1: 2;
  const double h_min = 1.;
  double s = (h_max-h_min)/(N-1);
  Uint i;
  
  path = make_directory(path_par,"EoS_Tests");
  
  /* values */
  /* if it is pwp */
  if (strstr_i(eos->type,"piecewise_polytropic") ||
      strstr_i(eos->type,"pwp"))
  {
    sprintf(file_name,"%s/%s.pwp",path,eos->description);
    file = Fopen(file_name,"w+");
    
    fprintf(file,"[1]:piece  [2]:Kappa  [3]:rho  [4]:gamma  [5]:a  [6]:h-1\n");
    for (i = 0; i < eos->N; ++i)  
      fprintf(file,"%u %0.15f %0.15f %0.15f %0.15f %0.15f\n",i,eos->K[i],eos->rho0_th[i],eos->gamma[i],eos->a[i],eos->h_th[i]-1);
    Fclose(file);
  }
    
  /* continuity */
  sprintf(file_name,"%s/%s",path,"pressure");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  pressure\n");  
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"%0.15f %0.15f\n",eos->h,eos->pressure(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s",path,"rest_mass_density");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  rest_mass_density\n");  
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"%0.15f %0.15f\n",eos->h,eos->rest_mass_density(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s",path,"energy_density");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  energy_density\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"%0.15f %0.15f\n",eos->h,eos->energy_density(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s",path,"specific_internal_energy");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  specific_internal_energy\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"%0.15f %0.15f\n",eos->h,eos->specific_internal_energy(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s",path,"drho0_dh");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  drho0_dh\n");  
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"%0.15f %0.15f\n",eos->h,eos->drho0_dh(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s",path,"de_dh");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  de_dh\n");  
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"%0.15f %0.15f\n",eos->h,eos->de_dh(eos));
  }
  Fclose(file);
  
  
  free_EoS(eos);
  free(path);
  UNUSED(grid);
  
}
