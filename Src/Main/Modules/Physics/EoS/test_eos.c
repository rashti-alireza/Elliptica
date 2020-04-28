/*
// Alireza Rashti
// June 2019
*/

#include "test_eos.h"

void test_EoS(Grid_T *const grid)
{
  EoS_T *eos = initialize_EoS();
  const char *const path_par = Pgets("output_directory_path");
  char *path,file_name[1000];
  FILE *file = 0;
  unsigned N = 100;
  const double h_max = eos->h_th != 0 ? eos->h_th[eos->N-1]+10: 10;
  const double h_min = 1;
  double s = (h_max-h_min)/(N-1);
  unsigned i;
  
  path = make_directory(path_par,"EoS_Tests");
  
  /* values */
  /* if it is pwp */
  if (strstr_i(eos->type,"piecewise_polytropic") ||
      strstr_i(eos->type,"pwp"))
  {
    sprintf(file_name,"%s/%s.pwp",path,eos->description);
    file = Fopen(file_name,"w+");
    IsNull(file);
    
    fprintf(file,"piece  Kappa         rho           gamma         a             h-1\n");
    for (i = 0; i < eos->N; ++i)  
      fprintf(file,"%u      %e  %e  %e  %e  %e\n",i,eos->K[i],eos->rho_th[i],eos->gamma[i],eos->a[i],eos->h_th[i]-1);
    fclose(file);
  }
    
  /* continuity */
  sprintf(file_name,"%s/%s",path,"pressure");
  file = Fopen(file_name,"w+");
  IsNull(file);
  fprintf(file,"# enthalpy   pressure\n");  
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"  %-7g    %-7g\n",eos->h,eos->pressure(eos));
  }
  fclose(file);
  
  sprintf(file_name,"%s/%s",path,"rest_mass_density");
  file = Fopen(file_name,"w+");
  IsNull(file);
  fprintf(file,"# enthalpy   rest_mass_density\n");  
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"  %-7g    %-7g\n",eos->h,eos->rest_mass_density(eos));
  }
  fclose(file);
  
  sprintf(file_name,"%s/%s",path,"energy_density");
  file = Fopen(file_name,"w+");
  IsNull(file);
  fprintf(file,"# enthalpy   energy_density\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"  %-7g    %-7g\n",eos->h,eos->energy_density(eos));
  }
  fclose(file);
  
  sprintf(file_name,"%s/%s",path,"drho_dh");
  file = Fopen(file_name,"w+");
  IsNull(file);
  fprintf(file,"# enthalpy   drho_dh\n");  
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"  %-7g    %-7g\n",eos->h,eos->drho_dh(eos));
  }
  fclose(file);
  
  sprintf(file_name,"%s/%s",path,"de_dh");
  file = Fopen(file_name,"w+");
  IsNull(file);
  fprintf(file,"# enthalpy   de_dh\n");  
  for (i = 0; i < N; ++i)
  {
    eos->h = 1+s*i;
    fprintf(file,"  %-7g    %-7g\n",eos->h,eos->de_dh(eos));
  }
  fclose(file);
  
  
  free_EoS(eos);
  free(path);
  UNUSED(grid);
  
}
