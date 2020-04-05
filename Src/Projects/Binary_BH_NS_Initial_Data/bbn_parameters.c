/*
// Alireza Rashti
// February 2020
*/
#include "bbn_parameters.h"

/* set default values of some of the parameters, unless they have been
// set by the input parameter file. */
void bbn_set_default_parameters(void)
{
  /* how the grid is covered */
  Pset_default("grid_kind","BBN_CubedSpherical_grid");
  
  /* how to take derivative */
  Pset_default("Derivative_Method","Spectral");
  
  /* how to perform interpolation */
  Pset_default("Interpolation_Method","Spectral");
  
  /* how to calculate Fourier coeffs */
  Pset_default("Fourier_Transformation_Method","RFT");
  
  /* how to calculate functional variations */
  Pset_default("dF/du_for_Newton_Method","Spectral");
  
  /* how to calcualte J in Jx = -F in Newton Method */
  Pset_default("Making_Jacobian_For_Newton_Method","Spectral");
  
  /* how to solve Jx = -F equations */
  Pset_default("Solving_Method","DDM_Schur_Complement");
  
  /* root finder method */
  Pset_default("RootFinder_Method","Steepest_Descent");
  
  /* maximum number of iteration for the root finder */
  Pset_default("RootFinder_Max_Number_of_Iteration","2000");
  
  /* stopping criteria for the root finder */
  Pset_default("RootFinder_Tolerance","10E-10");
 
  /* if stop = 1, the main iteration loop is halted */
  Pseti("STOP",0);
  
  /* total number of iterations that have been executed */
  Pseti("iteration_number",0);
  
  /* if _Aij is active for jacobian of psi equation it is 1
  // otherwise it is 0. the default is 0. it seems when it is active
  // the elliptic solve won't converge. */
  Pset_default("is_Aij_on_at_AH?","0");
  
  /* the tolerance of BH mass while adjusting the AH excision radius to
  // meet the target BH mass */
  Pset_default("BH_mass_tolerance","1e-2");
  
  /* the weight we use to update enthalpy */
  Pset_default("enthalpy_update_weight","1e-3");
  
  /* the weight we use to update NS surface */
  Pset_default("NS_surface_update_weight","1e-3");
  
  /* the weight we use to update BH excision radius */
  Pset_default("BH_r_excision_update_weight","1e-2");
  
  /* how to extrapolate fluid_fields outside the NS 
  // options:
  // "phi:exp_exp_continuity,enthalpy:Ylm" -->
  // 	for phi requires to be exponentially decreasing and be C^2
  //    and for enthalpy it requires to be sum{Clm * r^-(l+1) * Ylm}. */
  Pset_default("extrapolate_fluid_fields_method","phi:exp_continuity,enthalpy:Ylm");
}
