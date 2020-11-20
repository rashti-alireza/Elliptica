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
  //const int  neat = OPgeti(par,obj,"enthalpy_neat");
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    ONLY_IF_COVER(patch,obj);
    
    RELAX_UPDATE_FUNC(Tij_IF_CTS_nonflat_enthalpy,patch,enthalpy,W);
    
    //if (neat)
      //neat_enthalpy(patch);
    
    //Tij_IF_CTS_nonflat_rho0(patch);
    Tij_IF_CTS_nonflat_u0(patch);
    Tij_IF_CTS_nonflat_psi6J_Ui(patch);
    Tij_IF_CTS_nonflat_psi6E(patch);
    Tij_IF_CTS_nonflat_psi6S(patch);
    
    //update_partial_derivatives(enthalpy);
    //update_partial_derivatives(rho0);
  }
  
  FUNC_TOC
}
