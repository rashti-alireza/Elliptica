/*
// Alireza Rashti
// August 2019
*/

#include "Tij_header.h"


/* updating all matter related and their derivatives. */
void Tij_ideal_fluid_CTS_update(Obj_Man_T *const obj)
{
  FUNC_TIC
  
  Grid_T *const grid = obj->grid;  
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    ONLY_IF_COVER(patch,obj);
    
    Tij_IF_CTS_enthalpy(patch);
    Tij_IF_CTS_u0(patch);
    Tij_IF_CTS_psi6J_Ui(patch);
    Tij_IF_CTS_psi6E(patch);
    Tij_IF_CTS_psi6S(patch);
  }
  
  FUNC_TOC
}
