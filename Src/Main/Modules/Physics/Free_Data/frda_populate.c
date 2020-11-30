/*
// Alireza Rashti
// November 2020
*/

/* various functions and ways to populate free data */

#include "frda_populate.h"


void frda_populate_gConf_dgConf_igConf_KerrSchild(Physics_T *const phys)
{
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
}
