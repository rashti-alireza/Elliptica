/*
// Alireza Rashti
// November 2020
*/

/* ADM general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. */

#include "adm_main.h"

/* main function to issue command */
int adm_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case ADM_SET_PARAMS:
      ret = set_adm_params(phys);
    break;
    
    case ADM_ADD_FIELDS:
      ret = add_adm_fields(phys);
    break;
    
    case ADM_COMPUTE_CONSTRAINTS:
      ret = compute_ham_and_mom_constrains(phys);
    break;
    
    case ADM_UPDATE_AConfIJ:
      ret = compute_AConfIJ(phys);
    break;
   
    case ADM_UPDATE_KIJ:
      ret = compute_adm_KIJ(phys);
    break;
   
    case ADM_UPDATE_Kij:
      ret = compute_adm_Kij(phys);
    break;
    
    case ADM_UPDATE_gij:
      ret = compute_adm_gij(phys);
    break;
    
    case ADM_UPDATE_B1I:
      ret = compute_B1I(phys);
    break;
    
    case ADM_UPDATE_beta:
      ret = compute_beta(phys);
    break;
    
    case ADM_DOCTEST:
      ret = preform_adm_doctest(phys);
    break;

    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}


/* set default paramters */
static int set_adm_params(Physics_T *const phys)
{
  FUNC_TIC
  
  /* how to compute constraints:
  // options:
  // from_scratch: using beta,psi,alphaPsi and trK to make Kij, 
  //               then using psi, gConf to make g
  //               and then finally make R and using sources 
  //               to calucalte constraints.
  // from_identities: using AConfIJ and various other identities
  //                  to calculate constraints.
  // from_residuals : using residual of elliptic eqs to compute. */
  Pset_default(P_"constraints_method","from_scratch");
  
  /* computing AConf^{ij} = 1./(sigma)*(LConf(W)^{ij}) + MConf^{ij}:
  // options:
  // XCTS_MConfIJ0: => \bar{sigma} = alpha*psi^-6, W = beta and MConf^{ij} = 0.
  // SCTT_MConfIJ0: => \bar{sigma} = 1, W = beta and MConf^{ij} = 0. */
  Pset_default(P_"compute_AConfIJ","XCTS_MConfIJ0");
  
  /* computing K_{ij}:
  // options:
  // use_AIJ: => use formula K_{ij} = g_{il} g_{jk} A^{lk} +1/3 trK g_{ij}. */
  Pset_default(P_"compute_adm_Kdd_method","use_AIJ");
  
  /* computing K^{ij}:
  // options:
  // use_AIJ: => use formula K^{ij} = A^{ij} +1/3 trK g^{ij}.
  // use_Kij: => use formula K^{ij} = g^{il} g^{jk} *K_{lk}*/
  Pset_default(P_"compute_adm_Kuu_method","use_AIJ");
  
  /* computing B1^i for beta^i = B0^i + B1^i:
  // options:
  // zero:     => B1^i = 0
  // inspiral: => B1^i = omega*(r-r_CM) + v/D*(r-r_CM) */
  Pset_default(P_"B1I_form","inspiral");
  
  /* doctest for AConf^{ij} using adm_K{ij}:
  // options:
  // KerrSchild: compare computed adm_K_{ij} with analytic Kerr-Schild */
  Pset_default(P_"doctest_AConfIJ_compare","KerrSchild");
  
  if(Pcmps(P_"B1I_form","inspiral"))
  {
    adm_update_B1I_patch = adm_update_B1I_inspiral;
  }
  else if (Pcmps(P_"B1I_form","zero"))
  {
    adm_update_B1I_patch = adm_update_B1I_zero;
  }
  else
    Error0(NO_OPTION);
  
  
  if(Pcmps(P_"compute_AConfIJ","XCTS_MConfIJ0"))
  {
    adm_update_AConfIJ_patch = adm_update_AConfIJ_XCTS_MConfIJ0;
  }
  else if(Pcmps(P_"compute_AConfIJ","SCTT_MConfIJ0"))
  {
    adm_update_AConfIJ_patch = adm_update_AConfIJ_SCTT_MConfIJ0;
  }
  else
    Error0(NO_OPTION);
  
  if(Pcmps(P_"compute_adm_Kdd_method","use_AIJ"))
  {
    adm_update_adm_Kij_patch = adm_update_adm_Kij_useAIJ;
  }
  else
    Error0(NO_OPTION);
  
  if(Pcmps(P_"compute_adm_Kuu_method","use_AIJ"))
  {
    adm_update_adm_KIJ_patch = adm_update_adm_KIJ_useAIJ;
  }
  else if(Pcmps(P_"compute_adm_Kuu_method","use_Kij"))
  {
    adm_update_adm_KIJ_patch = adm_update_adm_KIJ_useKij;
  }
  else
    Error0(NO_OPTION);
  
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* add fields parameters */
static int add_adm_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  assert(phys->grid);
  
  adm_add_3plus1_fields(phys->grid);
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}

/* compute constraints */
static int compute_ham_and_mom_constrains(Physics_T *const phys)
{
  FUNC_TIC
  
  if(strstr_i(Pgets(P_"constraints_method"),"from_scratch"))
    adm_compute_constraints(phys,".*","from_scratch"  ,"ham1","mom1");
    
  if(strstr_i(Pgets(P_"constraints_method"),"from_identities"))
    adm_compute_constraints(phys,".*","from_identities","ham2","mom2");
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}

/* compute AConf^{ij} */
static int compute_AConfIJ(Physics_T *const phys)
{
  FUNC_TIC
  
  adm_update_AConfIJ(phys,".*");
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}

/* compute adm_K_{ij} */
static int compute_adm_Kij(Physics_T *const phys)
{
  FUNC_TIC
  
  adm_update_adm_Kij(phys,".*");
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}

/* compute adm_K^{ij} */
static int compute_adm_KIJ(Physics_T *const phys)
{
  FUNC_TIC
  
  adm_update_adm_KIJ(phys,".*");
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}

/* compute adm_g_{ij} */
static int compute_adm_gij(Physics_T *const phys)
{
  FUNC_TIC
  
  adm_update_adm_gij(phys,".*");
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}

/* compute B1 in beta = B0+B1 */
static int compute_B1I(Physics_T *const phys)
{
  FUNC_TIC
  
  adm_update_adm_B1I(phys,".*");
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}

/* compute shifts, beta = B0+B1 */
static int compute_beta(Physics_T *const phys)
{
  FUNC_TIC
  
  adm_update_beta(phys,".*");
  
  FUNC_TOC
  return EXIT_SUCCESS; 
}


/* performing some internal tests to make sure 
// things implemented correctly. */
static int preform_adm_doctest(Physics_T *const phys)
{
  FUNC_TIC
  
  /* test AConf^{ij} */
  if (DO) adm_doctest_AConfIJ(phys);
  
  
  FUNC_TOC
  return EXIT_SUCCESS;
}
