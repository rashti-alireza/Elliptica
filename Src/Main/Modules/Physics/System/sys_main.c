/*
// Alireza Rashti
// November 2020
*/

/* physical system general affairs. one can add new different 
// function readily by adding new parameter and the name of 
// the function as shown. */

#include "sys_main.h"

/* update stress energy tensor */
int sys_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case SYS_SET_PARAMS:
      ret = set_system_params(phys);
    break;
    
    case SYS_ADD_FIELDS:
      ret = add_system_fields(phys);
    break;
    
    case SYS_TUNE_P_ADM:
      ret = tune_system_ADM_momenta(phys);
    break;
    
    case SYS_INITIALIZE_FIELDS:
      ret = initialize_fields(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* tune adm momenta */
static int tune_system_ADM_momenta(Physics_T *const phys)
{
  FUNC_TIC
  
  int ret = EXIT_SUCCESS;
  
  ret = sys_tune_ADM_momenta(phys);
  
  FUNC_TOC
  return ret;
}

/* initialize fields for the very first time to start off the 
// initial data procedure, using known cases, like TOV, KerrSchild etc. */
static int initialize_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  if(phys->sys == SBH                            &&
     Pcmps(P_"initialize","one_exact_KerrSchild"))
  {
    if(Pcmps(P_"initialize_fields","XCTS"))
    {
      /* important to have dedicated BH physics to read correct parameters */
      Physics_T *const bh = init_physics(phys,BH);
      fd_populate_psi_alphaPsi_beta_KerrSchild
        (bh,".*","psi","alphaPsi","beta",0);
      free_physics(bh);
    }
    else
      Error0(NO_OPTION);
  }
  else if(phys->sys == SBH                            &&
          Pcmps(P_"initialize","one_exact_ConfKerrSchild"))
  {
    if(Pcmps(P_"initialize_fields","XCTS"))
    {
      /* important to have dedicated BH physics to read correct parameters */
      Physics_T *const bh = init_physics(phys,BH);
      fd_populate_psi_alphaPsi_beta_ConfKerrSchild
        (bh,".*","psi","alphaPsi","beta");
      free_physics(bh);
    }
    else
      Error0(NO_OPTION);
  }
  else if(phys->sys == SBH                          &&
          Pcmps(P_"initialize","one_exact_IsoSchild"))
  {
    if(Pcmps(P_"initialize_fields","XCTS"))
    {
      /* important to have dedicated BH physics to read correct parameters */
      Physics_T *const bh = init_physics(phys,BH);
      fd_populate_psi_alphaPsi_beta_IsoSchild
        (bh,".*","psi","alphaPsi","beta",0);
      free_physics(bh);
    }
    else
      Error0(NO_OPTION);
  }
  else if(phys->sys == SBH                          &&
          Pcmps(P_"initialize","one_exact_PGSchild"))
  {
    if(Pcmps(P_"initialize_fields","XCTS"))
    {
      /* important to have dedicated BH physics to read correct parameters */
      Physics_T *const bh = init_physics(phys,BH);
      fd_populate_psi_alphaPsi_beta_PGSchild
        (bh,".*","psi","alphaPsi","beta",0);
      free_physics(bh);
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* set default parameters. */
static int set_system_params(Physics_T *const phys)
{
  FUNC_TIC
  
  /* which fields to initialize:
  // options:
  // XCTS: psi, alphaPsi and beta.
  */
  Pset_default(P_"initialize_fields","XCTS");
  
  /* how to initialize fields
  // options:
  // one_exact_KerrSchild: use analytic values of KerrSchild BH
  // one_exact_IsoSchild : use analytic values of Schawrzchild in isotropic coords.
  // one_exact_PGSchild  : use analytic values of Schawrzchild in Painleve-Gullstrand coords.
  */
  Pset_default(P_"initialize","one_exact_KerrSchild");

  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* add fields. */
static int add_system_fields(Physics_T *const phys)
{
  FUNC_TIC
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

