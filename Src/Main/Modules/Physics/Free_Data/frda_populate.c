/*
// Alireza Rashti
// November 2020
*/

/* various functions and ways to populate free data */

#include "frda_populate.h"


/* compute confromal Ricci_{ij} and its trace */
void frda_conformal_Ricci(Physics_T *const phys,
                          const char *const ig,
                          const char *const Chris,
                          const char *const dChris,
                          const char *const RicciConf,
                          const char *const trRicciConf)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,".*");
  unsigned p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Ricci_3d(patch,ig,Chris,dChris,RicciConf,trRicciConf);
  }
  
  FUNC_TOC
}

/* compute Christoffel symbol compatible with given metric */
void frda_compatible_Christoffel_symbol(Physics_T *const phys,
                                        const char *const ig,
                                        const char *const dg,
                                        const char *const Chris
                                       )
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,".*");
  unsigned p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Christoffel_symbol_3d(patch,ig,dg,Chris);
  }
  
  FUNC_TOC
}

/* compute 1st derivative Christoffel symbol */
void frda_1st_derivative_Christoffel_symbol(Physics_T *const phys,
                                            const char *const dChris)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,".*");
  unsigned p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    dChristoffel_symbol_3d(patch,dChris);
  }
  
  FUNC_TOC
}

/* populate confromal metric, inverse of confromal metric 
// and first order derivative of confromal metric. */
void frda_populate_gConf_dgConf_igConf_KerrSchild(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,".*");
  const double BHx = Getd("center_x");
  const double BHy = Getd("center_y");
  const double BHz = Getd("center_z");
  unsigned p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    frda_kerr_schild_g_analytic(patch,BHx,BHy,BHz,"gConf");
    frda_kerr_schild_dg_analytic(patch,BHx,BHy,BHz,"dgConf");
    
    READ_v(gConf_D2D2)
    READ_v(gConf_D0D2)
    READ_v(gConf_D0D0)
    READ_v(gConf_D0D1)
    READ_v(gConf_D1D2)
    READ_v(gConf_D1D1)
    
    REALLOC_v_WRITE_v(igConf_U2U2)
    REALLOC_v_WRITE_v(igConf_U0U2)
    REALLOC_v_WRITE_v(igConf_U0U0)
    REALLOC_v_WRITE_v(igConf_U0U1)
    REALLOC_v_WRITE_v(igConf_U1U2)
    REALLOC_v_WRITE_v(igConf_U1U1)
    
    FOR_ALL_ijk
    {
      Matrix_Inverse_3x3_Symmetric_Field(gConf,D,igConf,U,ijk);
    }
  }
  
  FUNC_TOC
}
