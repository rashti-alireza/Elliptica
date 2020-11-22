/*
// Alireza Rashti
// August 2019
*/

#include "Tij_ifluid_cts_gConf_update.h"


/* updating all matter related and their derivatives. */
void Tij_idealfluid_CTS_gConf_update(Physics_T *const obj)
{
  FUNC_TIC
  
  Grid_T *const grid  = obj->grid;
  char opar[OPAR_LEN] = {'\0'};
  const double W  = Getd("enthalpy_update_weight");
  const int  neat = Geti("enthalpy_neat");
  const double Euler_const = Getd("Euler_equation_constant");
  unsigned p;
  
  printf(Pretty0"weight update  = %e\n"
         Pretty0"neat it?       = %d\n"
         Pretty0"Euler constant = %e\n",
         W,neat,Euler_const);
  fflush(stdout);
         
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if_not_cover(patch,obj)  continue;
    
    RELAX_UPDATE_FUNC(Tij_IF_CTS_gConf_enthalpy(patch,Euler_const),
                      patch,enthalpy,W);
    
    if (neat)
      Tij_neat_enthalpy(patch);
    
    Tij_eos_update_rho0(patch);
    Tij_IF_CTS_gConf_u0(patch);
    Tij_IF_CTS_gConf_derives(patch);
    /* sources */
    Tij_IF_CTS_gConf_psi6J_Ui(patch);
    Tij_IF_CTS_gConf_psi6E(patch);
    Tij_IF_CTS_gConf_psi6S(patch);
  }
  
  FUNC_TOC
}

/* updating rho0 using eos */
void Tij_eos_update_rho0(Patch_T *const patch)
{
  EoS_T *eos = initialize_EoS();
  READ_v(enthalpy)
  REALLOC_v_WRITE_v(rho0)
  const unsigned nn = patch->nn;
  unsigned ijk;

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
  free_EoS(eos);
}

/* after finding new NS surface, root finder might find h ~ 1
// so put it h = 1, to prevent nan in matter fields  */
void Tij_neat_enthalpy(Patch_T *const patch)
{
  WRITE_v(enthalpy)
  
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    if (LSSEQL(enthalpy[ijk],1))
      enthalpy[ijk] = 1;
  }

}

