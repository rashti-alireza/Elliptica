/*
// Alireza Rashti
// November 2020
*/


/* computing Hamiltonian and momentum constraints */


#include "adm_constraints.h"

/* measuring Hamiltonian constraint.
// method:
// =======
// from_scratch: using beta,psi,alphaPsi and trK to make Kij, 
//               then using psi, gConf to make g
//               and then finally make R and using sources 
//               to calucalte constraints.
// from_identities: using AConfIJ and various other identities
//                  to calculate constraints.
// */
void adm_compute_constraints(Physics_T *const phys,
                                  const char *const region,
                                  const char *const method,
                                  const char *const ham,
                                  const char *const mom)
{
  if (strcmp_i(method,"from_identities"))
  {
    printf(Pretty0"method: from_identities.\n");
    fflush(stdout);
    
    Grid_T *const grid = mygrid(phys,region);
    unsigned p;
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *const patch = grid->patch[p];
      
      adm_ham_and_mom_from_identities(patch,ham,mom);
      
    }
  }
  else if (strcmp_i(method,"from_scratch"))
  {
    Grid_T *const grid = mygrid(phys,region);
    unsigned p;
    
    printf(Pretty0"method: from_scratch.\n");
    fflush(stdout);
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *const patch = grid->patch[p];
      
      adm_ham_and_mom_from_scratch(patch,ham,mom);
      
    }
  }
  else
    Error0(NO_OPTION);
  
}
