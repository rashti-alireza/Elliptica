/*
// Alireza Rashti
// June 2019
*/

#include "bbn_initialize.h"

/* initialize this system according to the parameter file. 
// ->return value: the grid as a result of this initialization. */
Grid_T *bbn_initialize_fields_and_grid(void)
{
  Grid_T *grid = 0;
  
  /* if we use TOV and Kerr-Schil black hole approximation */
  if (strcmp_i(GetParameterS_E("BBHNS_initialization"),"TOV_KerrShild"))
    grid = TOV_KerrShild_approximation();
  else
    abortEr(NO_OPTION);
  
  return grid;   
}

/* use TOV and Kerr-Schil black hole approximation.
// ->return value: resultant grid from this approximation */
static Grid_T *TOV_KerrShild_approximation(void)
{
  Grid_T *grid = 0;
  
  /* solve fields for a TOV star located at left side of y axis */
  TOV_T *tov = TOV_init();
  tov->N = (unsigned)GetParameterI("left_NS_n_a");
  if (tov->N == INT_MAX)/* if not specifed override it */
    tov->N = (unsigned)GetParameterI("n_a");
  tov->N *= 10;/* using Composite Simpson's rule integral, for accuracy we increas N */
  tov->N += 1;/* using Composite Simpson's rule integral, N must be odd */
  
  /* for the TOV guess we approximatly take adm mass as the barionic mass of the NS */
  tov->ADM_m = GetParameterD_E("NS_initial_baryonic_mass");
  tov->description = "Estimating NS";
  tov = TOV_solution(tov);

  /* solve fields for Kerr Shild black hole */
  //KerrShild_solution();
  
  /* combining these two solution to initialize the grid */
  //grid = creat_grid_TOV_KerrShild(...,...);
  
  TOV_free(tov);
  //KerrSchil_free();
  return grid;
}

