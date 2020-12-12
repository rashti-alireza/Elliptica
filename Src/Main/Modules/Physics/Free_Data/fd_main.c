/*
// Alireza Rashti
// November 2020
*/

/* free data general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. */

#include "fd_main.h"

/* main function to issue command */
int fd_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case FREE_DATA_SET_PARAMS:
      ret = set_free_data_params(phys);
    break;
    
    case FREE_DATA_ADD_FIELDS:
      ret = add_free_data_fields(phys);
    break;
    
    case FREE_DATA_POPULATE:
      ret = populate_free_data(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}


/* set default paramters */
static int set_free_data_params(Physics_T *const phys)
{
  FUNC_TIC
  
  /* how to set confromal metic:
  // options:
  // flat:       gConf = delta_{ij}
  // KerrSchild: gConf = Kerr-Schild black hole
  // IsoSchild:  gConf = delta_{ij} for Schwarzchild in isotropic coords.
  // PGSchild:   gConf = delta_{ij} for Schwarzchild in Painleve-Gullstrand coords */
  Pset_default(P_"conformal_metric","KerrSchild");
  
  /* how to set Christoffel symbol:
  // options:
  // flat:       ChrisConf = 0
  // KerrSchild: ChrisConf made of gConf of Kerr-Schild black hole
  // IsoSchild:  ChrisConf = 0 for Schwarzchild in isotropic coords.
  // PGSchild:   ChrisConf = 0 for Schwarzchild in Painleve-Gullstrand coords */
  Pset_default(P_"conformal_Christoffel_symbol","KerrSchild");
  
  /* how to set trK = Tr(K_{ij})
  // options:
  // maximal:    trK = 0.
  // KerrSchild: trK = trK of Kerr-Schild black hole K_{ij} 
  // IsoSchild:  trK = 0 for Schwarzchild in isotropic coordinates.
  // PGSchild:   trK for Schwarzchild in Painleve-Gullstrand coords */
  Pset_default(P_"trK","KerrSchild");
  
  /* how to set conformal Ricci tensor
  // options:
  // flat:       RicciConf_{ij} = 0
  // KerrSchild: use Kerr-Schild black hole metric 
  // IsoSchild:  use Schwarzchild in isotropic coordinates.
  // PGSchild:   RicciConf_{ij} = 0 for Schwarzchild in Painleve-Gullstrand coords */
  Pset_default(P_"conformal_Ricci","KerrSchild");
  
 
  /* how to set MConf^{ij} in AConf^{ij} = 1/sigma (LConf W)^{ij} + MConf^{ij}
  // options:
  // zero:       MConf^{ij} = 0 */
  Pset_default(P_"MConfIJ","zero");
  
  
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* add fields parameters */
static int add_free_data_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  assert(phys->grid);
  
  fd_add_fields_gConf_dgConf_igConf(phys->grid);
  fd_add_fields_ChrisConf_dChrisConf(phys->grid);
  fd_add_fields_trK_dtrK(phys->grid);
  fd_add_fields_RicciConf(phys->grid);
  fd_add_fields_MConfIJ(phys->grid);
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}


/* populate free data */
static int populate_free_data(Physics_T *const phys)
{
  FUNC_TIC
  
  if (phys->sys                             == SBH         && 
      Pcmps(P_"conformal_metric"            ,"KerrSchild") &&
      Pcmps(P_"conformal_Christoffel_symbol","KerrSchild") &&
      Pcmps(P_"conformal_Ricci"             ,"KerrSchild") &&
      Pcmps(P_"trK"                         ,"KerrSchild") &&
      Pcmps(P_"MConfIJ"                     ,"zero"      )
     )
  {
    /* important to have dedicated BH physics to read correct parameters */
    Physics_T *const bh = init_physics(phys,BH);

    fd_populate_gConf_dgConf_igConf_KerrSchild(bh,".*","gConf",
                                              "igConf","dgConf");
    fd_compatible_Christoffel_symbol(bh,".*","igConf",
                                     "dgConf","ChrisConf");
    fd_1st_derivative_Christoffel_symbol(bh,".*","dChrisConf");

    fd_conformal_Ricci(bh,".*","igConf","ChrisConf","dChrisConf",
                       "RicciConf","trRicciConf");
    fd_extrinsic_curvature_KerrSchild(bh,".*","igConf","ChrisConf",
                                      "adm_Kij","trK","dtrK");
    free_physics(bh);
  }
  else if 
    (phys->sys                             == SBH        && 
     Pcmps(P_"conformal_metric"            ,"IsoSchild") &&
     Pcmps(P_"conformal_Christoffel_symbol","IsoSchild") &&
     Pcmps(P_"conformal_Ricci"             ,"IsoSchild") &&
     Pcmps(P_"trK"                         ,"IsoSchild") &&
     Pcmps(P_"MConfIJ"                     ,"zero"     )
    )
  {
    /* important to have dedicated BH physics to read correct parameters */
    Physics_T *const bh = init_physics(phys,BH);

    fd_populate_gConf_dgConf_igConf_IsoSchild(bh,".*","gConf",
                                              "igConf","dgConf");
    fd_compatible_Christoffel_symbol(bh,".*","igConf",
                                    "dgConf","ChrisConf");
    fd_1st_derivative_Christoffel_symbol(bh,".*","dChrisConf");

    fd_conformal_Ricci(bh,".*","igConf","ChrisConf","dChrisConf",
                       "RicciConf","trRicciConf");
    fd_extrinsic_curvature_IsoSchild(bh,".*","igConf","ChrisConf",
                                     "adm_Kij","trK","dtrK");
    free_physics(bh);
  }
  else if 
    (phys->sys                             == SBH       && 
     Pcmps(P_"conformal_metric"            ,"PGSchild") &&
     Pcmps(P_"conformal_Christoffel_symbol","PGSchild") &&
     Pcmps(P_"conformal_Ricci"             ,"PGSchild") &&
     Pcmps(P_"trK"                         ,"PGSchild") &&
     Pcmps(P_"MConfIJ"                     ,"zero"    )
    )
  {
    /* important to have dedicated BH physics to read correct parameters */
    Physics_T *const bh = init_physics(phys,BH);

    fd_populate_gConf_dgConf_igConf_PGSchild(bh,".*","gConf",
                                             "igConf","dgConf");
    fd_compatible_Christoffel_symbol(bh,".*","igConf",
                                    "dgConf","ChrisConf");
    fd_1st_derivative_Christoffel_symbol(bh,".*","dChrisConf");

    fd_conformal_Ricci(bh,".*","igConf","ChrisConf","dChrisConf",
                       "RicciConf","trRicciConf");
    fd_extrinsic_curvature_PGSchild(bh,".*","igConf","ChrisConf",
                                    "adm_Kij","trK","dtrK");
    free_physics(bh);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}
