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
  // flat:       gConf = delta_{ij}.
  // KerrSchild: gConf = Kerr-Schild black hole.
  // ConfKerrSchild: gConf = Kerr-Schild black hole which decomposed conformally so det(gConf) = 1.
  // IsoSchild:  gConf = delta_{ij} for Schwarzchild in isotropic coords.
  // PGSchild:   gConf = delta_{ij} for Schwarzchild in Painleve-Gullstrand coords.
  // w1*flat+w2*KerrSchild: gConf_{ij} = w1*delta_{ij} + w2*gKS_{ij}. */
  Pset_default(P_"conformal_metric","KerrSchild");
  
  /* how to set Christoffel symbol:
  // options:
  // flat:       ChrisConf = 0.
  // KerrSchild: ChrisConf made of gConf of Kerr-Schild black hole.
  // ConfKerrSchild: ChrisConf made of gConf of ConfKerrSchild black hole.
  // IsoSchild:  ChrisConf = 0 for Schwarzchild in isotropic coords.
  // PGSchild:   ChrisConf = 0 for Schwarzchild in Painleve-Gullstrand coords.
  // w1*flat+w2*KerrSchild: for gConf_{ij} = w1*delta_{ij} + w2*gKS_{ij}. */
  Pset_default(P_"conformal_Christoffel_symbol","KerrSchild");
  
  /* how to set trK = Tr(K_{ij})
  // options:
  // zero:       trK = 0.
  // KerrSchild: trK = trK of Kerr-Schild black hole K_{ij}.
  // IsoSchild:  trK = 0 for Schwarzchild in isotropic coordinates.
  // PGSchild:   trK for Schwarzchild in Painleve-Gullstrand coords.
  // w*KerrSchild: trK = w*(KerrSchild trK). */
  Pset_default(P_"trK","KerrSchild");
  
  /* how to set conformal Ricci tensor
  // options:
  // flat:       RicciConf_{ij} = 0.
  // KerrSchild: use Kerr-Schild black hole metric .
  // ConfKerrSchild: use ConfKerrSchild black hole metric .
  // IsoSchild:  use Schwarzchild in isotropic coordinates.
  // PGSchild:   RicciConf_{ij} = 0 for Schwarzchild in Painleve-Gullstrand coords.
  // w1*flat+w2*KerrSchild: for gConf_{ij} = w1*delta_{ij} + w2*gKS_{ij}. */
  Pset_default(P_"conformal_Ricci","KerrSchild");
  
  /* how to set MConf^{ij} in AConf^{ij} = 1/sigma (LConf W)^{ij} + MConf^{ij}
  // options:
  // zero:       MConf^{ij} = 0. */
  Pset_default(P_"MConfIJ","zero");
  
  /* soft parameters:
  // ================
  //
  // 
  // name: "RollOff_function"
  // the transition function to stich free data (same as w1 and w2).
  // options:
  // o.  "exp(-lambda*(r/rmax)^p):r<rmax". # 0 for r >= rmax.
  // o.  "exp(-lambda*(r/rmax)^p)". # no condition on r.
  //
  //
  // name "RollOff_lambda"
  // lambda in "RollOff_function".
  // options: 
  // o. "|(r-rmin)/(rmax-r)|". # rmin,rmax are apparent horizon and 
  //                           # roll-off radii, respectively.
  // o. "constant_1". # a constant function
  //
  // name: "RollOff_rmax"
  // rmax in "RollOff_function".
  //
  //
  // name: "RollOff_power"
  // p in "RollOff_function".
  //
  //
  //
  */ 
  
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* add fields parameters */
static int add_free_data_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  assert(phys->grid);
  
  fd_add_fields_gConf_igConf_dgConf(phys->grid);
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

    fd_populate_gConf_igConf_dgConf_KerrSchild(bh,".*","gConf",
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

    fd_populate_gConf_igConf_dgConf_IsoSchild(bh,".*","gConf",
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

    fd_populate_gConf_igConf_dgConf_PGSchild(bh,".*","gConf",
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
  else if (phys->sys                             == SBH             && 
           Pcmps(P_"conformal_metric"            ,"ConfKerrSchild") &&
           Pcmps(P_"conformal_Christoffel_symbol","ConfKerrSchild") &&
           Pcmps(P_"conformal_Ricci"             ,"ConfKerrSchild") &&
           Pcmps(P_"trK"                         ,"KerrSchild")     &&
           Pcmps(P_"MConfIJ"                     ,"zero"      ))
  {
    /* important to have dedicated BH physics to read correct parameters */
    Physics_T *const bh = init_physics(phys,BH);
    
    /* first make trK which is KerrSchild:  */
    fd_populate_gConf_igConf_dgConf_KerrSchild(bh,".*","gConf",
                                               "igConf","dgConf");
    fd_compatible_Christoffel_symbol(bh,".*","igConf",
                                     "dgConf","ChrisConf");
    fd_extrinsic_curvature_KerrSchild(bh,".*","igConf","ChrisConf",
                                      "adm_Kij","trK","dtrK");
    /* then make the others, NOTE gConf is ready */
    fd_populate_gConf_igConf_dgConf_ConfKerrSchild
                 (bh,".*","gConf","gConf","igConf","dgConf");
    fd_compatible_Christoffel_symbol(bh,".*","igConf",
                                     "dgConf","ChrisConf");
    fd_1st_derivative_Christoffel_symbol(bh,".*","dChrisConf");

    fd_conformal_Ricci(bh,".*","igConf","ChrisConf","dChrisConf",
                       "RicciConf","trRicciConf");
    free_physics(bh);
  }
  else if 
   (phys->sys                             == BHNS                  && 
   Pcmps(P_"conformal_metric"            ,"w1*flat+w2*KerrSchild") &&
   Pcmps(P_"conformal_Christoffel_symbol","w1*flat+w2*KerrSchild") &&
   Pcmps(P_"conformal_Ricci"             ,"w1*flat+w2*KerrSchild") &&
   Pcmps(P_"trK"                         ,"w*KerrSchild")          &&
   Pcmps(P_"MConfIJ"                     ,"zero")                   )
  {
    /* important to have dedicated BH physics to read correct parameters */
    Physics_T *const bh = init_physics(phys,BH);

    /* first must make KerrSchild */
    fd_populate_gConf_igConf_dgConf_KerrSchild(bh,".*","gConf",
                                              "igConf","dgConf");
    fd_compatible_Christoffel_symbol(bh,".*","igConf",
                                     "dgConf","ChrisConf");
    fd_extrinsic_curvature_KerrSchild(bh,".*","igConf","ChrisConf",
                                      "adm_Kij","trK",0);
    
    /* modify metric to be "w1*flat+w2*KerrSchild" */
    fd_modify_gConf_igConf_dgConf_to_w1flat_w2KS(bh,".*","gConf",
                                                 "igConf","dgConf");
    
    fd_compatible_Christoffel_symbol(phys,".*","igConf",
                                     "dgConf","ChrisConf");
    fd_1st_derivative_Christoffel_symbol(phys,".*","dChrisConf");

    fd_conformal_Ricci(phys,".*","igConf","ChrisConf","dChrisConf",
                       "RicciConf","trRicciConf");
    
    /* modify trK to w*trK and computer its derivatives */
    fd_modify_trK_to_wtrK_compute_dtrK(bh,".*","trK","dtrK");

    free_physics(bh);
  }
  else if 
   (
   Pcmps(P_"conformal_metric"            ,"flat") &&
   Pcmps(P_"conformal_Christoffel_symbol","flat") &&
   Pcmps(P_"conformal_Ricci"             ,"flat") &&
   Pcmps(P_"trK"                         ,"zero") &&
   Pcmps(P_"MConfIJ"                     ,"zero")
   )
  {
    fd_populate_gConf_igConf_dgConf_flat(phys,".*","gConf",
                                              "igConf","dgConf");
    fd_compatible_Christoffel_symbol(phys,".*","igConf",
                                     "dgConf","ChrisConf");
    fd_1st_derivative_Christoffel_symbol(phys,".*","dChrisConf");

    fd_conformal_Ricci(phys,".*","igConf","ChrisConf","dChrisConf",
                       "RicciConf","trRicciConf");

    fd_trace_extrinsic_curvature_zero(phys,".*","trK","dtrK");
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}
