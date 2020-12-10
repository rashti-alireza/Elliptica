/*
// Alireza Rashti
// December 2020
*/


/* various doc test to make sure things works correctly */

#include "adm_tests.h"

/* test Aconf^{ij}.
// NOTE: it overwrites many fields thus fields are not supposed to be 
// used again and one must exit from project or reset fields. */
void adm_doctest_AConfIJ(Physics_T *const phys)
{
  FUNC_TIC
  
  if (phys->sys                       == SBH         &&
      Pcmps(P_"doctest_AConfIJ_compare","KerrSchild"))
  {
    /* important to have dedicated BH physics to read correct parameters */
    Physics_T *const bh = init_physics(phys,BH);
    Grid_T *const grid  = mygrid(phys,".*");
    
    /* calculate analytic adm_K_{ij} name test_adm_Kij */
    add_3x3_symmetric_field(grid,"test_adm_Kij","down");
    
    fd_populate_gConf_dgConf_igConf_KerrSchild(bh,".*","gConf",
                                                "igConf","dgConf");
                                                
    fd_1st_derivative_Christoffel_symbol(bh,".*","dChrisConf");
    
    fd_conformal_Ricci(bh,".*","igConf","ChrisConf","dChrisConf",
                         "RicciConf","trRicciConf");
    
    fd_extrinsic_curvature_KerrSchild(bh,".*","igConf","ChrisConf",
                                        "test_adm_Kij","trK","dtrK");
    
    /* compute adm_K_{ij} using AConf^{ij} */
    physics(phys,ADM_UPDATE_Kij);
    
    /* compare */
    diff_3x3_symmetric_fields
      (grid,"test_adm_Kij","adm_Kij","down",1);
    
    /* remove test_adm_Kij */
    FOR_ALL_p(grid->np)
    {
      Patch_T *patch = grid->patch[p];
      remove_field_with_regex(patch,"^test_adm_Kij_D.*");
    }
    
    free_physics(bh);
  }
  else
    Error0(NO_OPTION);
    
  FUNC_TOC
}
