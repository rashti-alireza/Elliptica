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
  const double h_max = PgetdEZ(Ftype("NS_EoS_enthalpy_ceiling")) == DBL_MAX ?
                       3.00 : Getd(P_"enthalpy_ceiling");
  const double h_min = PgetdEZ(Ftype("NS_EoS_enthalpy_floor")) == DBL_MAX ?
                       1.00 : Getd(P_"enthalpy_floor");
  double s = (h_max-h_min)/(N-1);
  Uint i;
  path = make_directory(path_par,"EoS_Tests");
  
  /* values */
  /* if it is pwp */
  if (strstr_i(eos->type,"piecewise_polytropic") ||
      strstr_i(eos->type,"pwp"))
  {
    sprintf(file_name,"%s/%s_pwp.txt",path,eos->description);
    file = Fopen(file_name,"w+");
    
    fprintf(file,"[1]:piece  [2]:Kappa  [3]:rho  [4]:gamma  [5]:a  [6]:h-1\n");
    for (i = 0; i < eos->N; ++i)  
    {
      fprintf(file,"%u %0.15f %0.15f %0.15f %0.15f %0.15f\n",
        i,eos->K[i],eos->rho0_th[i],eos->gamma[i],eos->a[i],eos->h_th[i]-1);
    }
    Fclose(file);
  }
    
  /* continuity */
  sprintf(file_name,"%s/%s.txt",path,"pressure");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  pressure\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->pressure(eos));
  }
  Fclose(file);

  sprintf(file_name,"%s/%s.txt",path,"rest_mass_density");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  rest_mass_density\n");  
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->rest_mass_density(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s.txt",path,"energy_density");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  energy_density\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->energy_density(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s.txt",path,"specific_internal_energy");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  specific_internal_energy\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->specific_internal_energy(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s.txt",path,"drho0_dh");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  drho0_dh\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->drho0_dh(eos));
  }
  Fclose(file);
      
  sprintf(file_name,"%s/%s.txt",path,"de_dh");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  de_dh\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->de_dh(eos));
  }
  Fclose(file);
 
  if (strcmp_i(eos->type, "tabular") || 
      strcmp_i(eos->type, "tab")     || 
      strcmp_i(eos->type, "table"))
  {
    sprintf(file_name,"%s/%s_tabular.txt",path,eos->description);
    file = Fopen(file_name,"w");
    // header
    fprintf(file,"# [1]:enthalpy(h) [2]:pressure(p) [3]:rest_mass_density(rho0) "
                 "[4]:energy_density(e) [5]:specific_internal_energy(e0) "
                 "[6]:drho0/dh [7]:de/dh\n");
                 
    for (i = 0; i < N; ++i)
    {
      eos->h = h_min+s*i;
      fprintf(file,"%+0.15e ",eos->h);
      fprintf(file,"%+0.15e ",eos->pressure(eos));
      fprintf(file,"%+0.15e ",eos->rest_mass_density(eos));
      fprintf(file,"%+0.15e ",eos->energy_density(eos));
      fprintf(file,"%+0.15e ",eos->specific_internal_energy(eos));
      fprintf(file,"%+0.15e ",eos->drho0_dh(eos));
      fprintf(file,"%+0.15e",eos->de_dh(eos));
      fprintf(file,"\n");
    }
    Fclose(file);
  }
 
  free_EoS(eos);
  free(path);
  UNUSED(grid);
}
