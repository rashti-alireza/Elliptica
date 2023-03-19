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
  else if(phys->sys == BHNS                    &&
          Pcmps(P_"initialize","TOV+KerrSchild"))
  {
    if(Pcmps(P_"initialize_fields","XCTS"))
    {
      /* add some auxiliary fields */
      add_aux_fields(mygrid(phys,".*"),
        "psi_ks,alphaPsi_ks,psi_tov,alphaPsi_tov");
      
      /* important to have dedicated BH physics to read correct parameters */
      Physics_T *const bh = init_physics(phys,BH);
      fd_populate_psi_alphaPsi_beta_KerrSchild
        (bh,".*","psi_ks","alphaPsi_ks","beta",0);
      free_physics(bh);
      
      /* important to have dedicated NS physics to read correct parameters */
      Physics_T *const ns = init_physics(phys,NS);
      star_populate_psi_alphaPsi_matter_fields_TOV
        (ns,".*","psi_tov","alphaPsi_tov","enthalpy","rho0","phi","W");
      /* alse we need NS spin vector */
      star_W_spin_vector_idealfluid_update(ns,"NS");
      free_physics(ns);
      
      /* beta, phi,W and rho0 remain intact */
      /* superimpose add f = f1 + f2 -1. */
      superimpose_simple(mygrid(phys,".*"),
                         "psi","psi_tov","psi_ks",-1.);
      superimpose_simple(mygrid(phys,".*"),
                        "alphaPsi","alphaPsi_tov","alphaPsi_ks",-1.);
      
      /* remove auxiliary fields */
      remove_aux_fields(mygrid(phys,".*"),
        "psi_ks,alphaPsi_ks,psi_tov,alphaPsi_tov");
      
    }
    else
        Error0(NO_OPTION);
      
  }
  else if(phys->sys == BHNS                    &&
          Pcmps(P_"initialize","TOV+IsoSchild"))
  {
    if(Pcmps(P_"initialize_fields","XCTS"))
    {
      /* add some auxiliary fields */
      add_aux_fields(mygrid(phys,".*"),
        "psi_is,alphaPsi_is,psi_tov,alphaPsi_tov");
      
      /* important to have dedicated BH physics to read correct parameters */
      Physics_T *const bh = init_physics(phys,BH);
      fd_populate_psi_alphaPsi_beta_IsoSchild
        (bh,".*","psi_is","alphaPsi_is","beta",0);
      free_physics(bh);
      
      /* important to have dedicated NS physics to read correct parameters */
      Physics_T *const ns = init_physics(phys,NS);
      star_populate_psi_alphaPsi_matter_fields_TOV
        (ns,".*","psi_tov","alphaPsi_tov","enthalpy","rho0","phi","W");
      /* alse we need NS spin vector */
      star_W_spin_vector_idealfluid_update(ns,"NS");
      free_physics(ns);
      
      /* beta, phi,W and rho0 remain intact */
      /* superimpose add f = f1 + f2 -1. */
      superimpose_simple(mygrid(phys,".*"),
                         "psi","psi_tov","psi_is",-1.);
      superimpose_simple(mygrid(phys,".*"),
                        "alphaPsi","alphaPsi_tov","alphaPsi_is",-1.);
      
      /* remove auxiliary fields */
      remove_aux_fields(mygrid(phys,".*"),
        "psi_is,alphaPsi_is,psi_tov,alphaPsi_tov");
    }
    else
      Error0(NO_OPTION);
  }
  else if(phys->sys == NSNS                &&
          Pcmps(P_"initialize","TOV+TOV"))
  {
    if(Pcmps(P_"initialize_fields","XCTS"))
    {
      /* add some auxiliary fields */
      add_aux_fields(mygrid(phys,".*"),
        "psi_tov1,alphaPsi_tov1,psi_tov2,alphaPsi_tov2");
      
      /* important to have dedicated NS physics to read correct parameters */
      Physics_T *const ns1 = init_physics(phys,NS1);
      star_populate_psi_alphaPsi_matter_fields_TOV
        (ns1,".*","psi_tov1","alphaPsi_tov1","enthalpy","rho0","phi","W");
      /* alse we need NS spin vector */
      star_W_spin_vector_idealfluid_update(ns1,"NS1");
      free_physics(ns1);

      Physics_T *const ns2 = init_physics(phys,NS2);
      star_populate_psi_alphaPsi_matter_fields_TOV
        (ns2,".*","psi_tov2","alphaPsi_tov2","enthalpy","rho0","phi","W");
      /* alse we need NS spin vector */
      star_W_spin_vector_idealfluid_update(ns2,"NS2");
      free_physics(ns2);
      
      /* phi,W and rho0 remain intact */
      /* superimpose add f = f1 + f2 -1. */
      superimpose_simple(mygrid(phys,".*"),
                         "psi","psi_tov1","psi_tov2",-1.);
      superimpose_simple(mygrid(phys,".*"),
                        "alphaPsi","alphaPsi_tov1","alphaPsi_tov2",-1.);
      /* set beta^i = 0. */
      superimpose_simple(mygrid(phys,".*"),
                        "beta_U0",0,0,0);
      superimpose_simple(mygrid(phys,".*"),
                        "beta_U1",0,0,0);
      superimpose_simple(mygrid(phys,".*"),
                        "beta_U2",0,0,0);
      
      /* remove auxiliary fields */
      remove_aux_fields(mygrid(phys,".*"),
        "psi_tov1,alphaPsi_tov1,psi_tov2,alphaPsi_tov2");
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
  // TOV+KerrSchild      : superpose TOV solution + KerrSchild solution.
  // TOV+IsoSchild       : superpose TOV solution + IsoSchild solution.
  */
  Pset_default(P_"initialize","one_exact_KerrSchild");

  /* SOFT params: */
  
  /* how to adjust P_ADM:
  // name: "P_ADM_control_method"
  // options:
  // none: do nothing
  // adjust(x_CM): adjust x_CM to drive P_y = 0.
  // adjust(y_CM): adjust y_CM to drive P_x = 0.
  // also one can combine things for instance: adjust(x_CM,y_CM).
  */
  
  /* update weight for P_ADM
  // name:
  // "P_ADM_control_update_weight"
  */
  
  /* how sensitive gets to adjust P_ADM
  // name:
  // "P_ADM_control_tolerance"
  */
  
  
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

