/*
// Alireza Rashti
// November 2020
*/

/* free data general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. */

#include "frda_main.h"

/* main function to issue command */
int frda_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case FREE_DATA_ADD_PARAMS:
      ret = add_free_data_params(phys);
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


/* add default parameters */
static int add_free_data_params(Physics_T *const phys)
{
  FUNC_TIC
  
  /* how to set confromal metic:
  // options:
  // flat:       gConf = delta_{ij}
  // KerrSchild: gConf = Kerr-Schild black hole */
  Pset_default(P_"conformal_metric","KerrSchild");
  
  /* how to set Christoffel symbol:
  // options:
  // flat:       GhrisConf = 0
  // KerrSchild: GhrisConf made of gConf of Kerr-Schild black hole */
  Pset_default(P_"conformal_Christoffel_symbol","KerrSchild");
  
  /* how to set trK = Tr(K_{ij})
  // options:
  // maximal:    trK = 0.
  // KerrSchild: trK = trK of Kerr-Schild black hole K_{ij} */
  Pset_default(P_"trK","KerrSchild");
  
  /* how to set conformal Ricci tensor
  // options:
  // flat:       RicciConf_{ij} = 0
  // KerrSchild: RicciConf = made out of Kerr-Schild black hole metric */
  Pset_default(P_"conformal_Ricci","KerrSchild");
  
  
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* add fields parameters */
static int add_free_data_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  assert(phys->grid);
  
  frda_add_fields_gConf_dgConf_igConf(phys->grid);
  frda_add_fields_ChrisConf_dChrisConf(phys->grid);
  frda_add_fields_trK_dtrK(phys->grid);
  frda_add_fields_RicciConf(phys->grid);
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}


/* populate free data */
static int populate_free_data(Physics_T *const phys)
{
  FUNC_TIC
  
  if (Pcmps(P_"conformal_metric"            ,"KerrSchild") &&
      Pcmps(P_"conformal_Christoffel_symbol","KerrSchild") &&
      Pcmps(P_"trK"                         ,"KerrSchild") &&
      Pcmps(P_"conformal_Ricci"             ,"KerrSchild") 
     )
  {
    frda_populate_gConf_dgConf_igConf_KerrSchild(phys);
    
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}
