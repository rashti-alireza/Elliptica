/*
// Alireza Rashti
// November 2020
*/

/* Tij general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. 
// also in a project first mount must be called and then update. */

#include "Tij_main.h"

/* update stress energy tensor */
int Tij_update(Obj_Man_T *const obj)
{
  if (Pcmps("Tij_fluid","ideal_fluid") && 
      Pcmps("Tij_decomposition","CTS") &&
      Pcmps("Tij_gConf","non_flat"))
  {
    Tij_idealfluid_CTS_nonflat_update(obj);
  }
  else
    Error0(NO_OPTION);
  
  return EXIT_SUCCESS;
}

/* adding default parameters and fields. */
int Tij_mount(Obj_Man_T *const obj)
{
  /* decomposition type:
  // CTS (conformal thin sandwich): like: Phys. Rev. D 100, 124046  */
  Pset_default("Tij_decomposition","CTS");
  
  /* fluid type:
  // ideal_fluid: like: Phys. Rev. D 100, 124046  */
  Pset_default("Tij_fluid","ideal_fluid");
  
  /* conformal metric type: 
  // flat     => gConf  = delta_{ij},
  // non_flat => gConf != delta_{ij}. */
  Pset_default("Tij_gConf","non_flat");
  
  if (Pcmps("Tij_fluid","ideal_fluid") && 
      Pcmps("Tij_decomposition","CTS") &&
      Pcmps("Tij_gConf","non_flat"))
  {
    Tij_idealfluid_CTS_nonflat_add_fields(obj->grid);
  }
  else
    Error0(NO_OPTION);
  
  return EXIT_SUCCESS;
}
