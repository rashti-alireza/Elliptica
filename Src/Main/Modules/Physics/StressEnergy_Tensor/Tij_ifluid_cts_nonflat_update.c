/*
// Alireza Rashti
// August 2019
*/

#include "Tij_header.h"


/* updating all matter related and their derivatives. */
void Tij_idealfluid_CTS_nonflat_update(Obj_Man_T *const obj)
{
  FUNC_TIC
  
  Grid_T *const grid  = obj->grid;
  char par[OPAR_SIZE] = {'\0'};
  const double W  = OPgetd(par,obj,"enthalpy_update_weight");
  const int  neat = OPgeti(par,obj,"enthalpy_neat");
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    ONLY_IF_COVER(patch,obj);
    
    RELAX_UPDATE_FUNC(Tij_IF_CTS_nonflat_enthalpy,patch,enthalpy,W);
    
    if (neat)
      Tij_neat_enthalpy(patch);
    
    Tij_eos_update_rho0(patch);
    Tij_IF_CTS_nonflat_u0(patch);
    Tij_IF_CTS_nonflat_derives(patch);
    /* sources */
    Tij_IF_CTS_nonflat_psi6J_Ui(patch);
    Tij_IF_CTS_nonflat_psi6E(patch);
    Tij_IF_CTS_nonflat_psi6S(patch);
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