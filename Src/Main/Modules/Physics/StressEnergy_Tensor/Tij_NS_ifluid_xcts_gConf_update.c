/*
// Alireza Rashti
// August 2019
*/

#include "Tij_NS_ifluid_xcts_gConf_update.h"


/* updating all matter related and their derivatives. */
void Tij_NS_idealfluid_XCTS_gConf_update(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid  = mygrid(phys,"NS");
  const double W  = Getd("enthalpy_update_weight");
  const int  neat = strstr_i(Gets("enthalpy_neat"),"yes");
  const double Euler_const = Getd("Euler_equation_constant");
  
  printf(Pretty0"weight update  = %e\n"
         Pretty0"neat it?       = %d\n"
         Pretty0"Euler constant = %e\n",
         W,neat,Euler_const);
  fflush(stdout);
         
  FOR_ALL_p(grid->np)
  {
    Patch_T *patch = grid->patch[p];
    EoS_T *eos     = init_EoS(phys);
    
    RELAX_UPDATE_FUNC(Tij_NS_IF_XCTS_gConf_enthalpy(patch,Euler_const),
                      patch,enthalpy,W);
    
    if (neat)
      Tij_NS_neat_enthalpy(patch);
    
    Tij_NS_eos_update_rho0(patch,eos);
    Tij_NS_IF_XCTS_gConf_u0(patch);
    Tij_NS_IF_XCTS_gConf_derives(patch);
    /* sources */
    Tij_NS_IF_XCTS_gConf_psi6J_Ui(patch);
    Tij_NS_IF_XCTS_gConf_psi6E(patch,eos);
    Tij_NS_IF_XCTS_gConf_psi6S(patch,eos);
    
    free_EoS(eos);
  }
  
  FUNC_TOC
}

/* updating rho0 using eos */
void Tij_NS_eos_update_rho0(Patch_T *const patch,EoS_T *const eos)
{
  READ_v(enthalpy)
  REALLOC_v_WRITE_v(rho0)
  const Uint nn = patch->nn;
  Uint ijk;

  for (ijk = 0; ijk < nn; ++ijk)
  {
    eos->h    = enthalpy[ijk];
    rho0[ijk] = eos->rest_mass_density(eos);
    
    /* make sure rho0 won't get nan due to Gibbs phenomena or 
    // when finding NS surface. */
    if (!isfinite(rho0[ijk]))
    {
      //printf("Put rho0(h = %g) = 0.\n",enthalpy[ijk]);
      rho0[ijk] = 0;
    }
  }
}

/* after finding new NS surface, root finder might find h ~ 1
// so put it h = 1, to prevent nan in matter fields  */
void Tij_NS_neat_enthalpy(Patch_T *const patch)
{
  WRITE_v(enthalpy)
  
  const Uint nn = patch->nn;
  Uint ijk;
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    if (LSSEQL(enthalpy[ijk],1))
      enthalpy[ijk] = 1;
  }

}

