/*
// Alireza Rashti
// August 2019
*/

#include "Tij_IdealFluid_3plus1_decomposition.h"

/* building u0, _J^i, _E and _S */
void Tij_IF_CTS_psi6Sources(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Building u^0, _J^i, _E and _S ...\n");
  
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    Tij_IF_CTS_u0(patch);
    Tij_IF_CTS_psi6J_Ui(patch);
    Tij_IF_CTS_psi6E(patch);
    Tij_IF_CTS_psi6S(patch);
  }
  
  printf("} Building u^0, _J^i, _E and _S ==> Done.\n");
  pr_clock();
  pr_line_custom('=');

}
