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
  
  /* umfpack settings for size int (~2GB matrix) or long int 
  // options: "0" for int, otherwise long int */
  Pset_default("Solving_UMFPACK_size","0");
  
  /* umfpack max iter. refinement step: the bigger the slower and precise
  // options: integer number 0,1,2,... */
  Pset_default("Solving_UMFPACK_refinement_step","0");
  
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
  Pset_default("is_Aij_on_at_AH?","1");
  
  /* the tolerance of BH mass while adjusting the AH excision radius to
  // meet the target BH mass */
  Pset_default("BH_mass_tolerance","1e-2");
  
  /* the weight we use to update enthalpy */
  Pset_default("NS_enthalpy_update_weight","1e-3");
  
  /* the weight we use to update NS surface */
  Pset_default("NS_surface_update_weight","1e-3");
  
  /* the weight we use to update BH excision radius */
  Pset_default("BH_r_excision_update_weight","1e-2");
  
  /* how to extrapolate fluid_fields outside the NS 
  // options:
  // 1. "phi:exp_exp_continuity,enthalpy:Ylm" -->
  // 	for phi requires to be exponentially decreasing and be C^2
  //    and for enthalpy it requires to be sum{Clm * r^-(l+1) * Ylm}.
  // 2. "slop_method" -->
  //    required to have C^1 field across the boundary.
  */
  Pset_default("extrapolate_fluid_fields_method","slop_method");
  
  /* how to extrapolate fields inside the BH 
  // options:
  // 1. Ylm   : using Ylm expansion
  // 2. linear: simply using a linear function.
  //    this is used during the elliptic solve.
  // 3. WTGR  : using the method developed by Wolfgang and George. 
  //    it's C^1 continues.
  // 4. ChebTn_Ylm: using 3rd order ChebTn in r and Ylm in 
  //    directions. it's C^2 continues.
  // 5. 4th_Poly_Ylm: 4th order Polynomial in r and Ylm in 
  //    angular directions. it's C^2 continues. */
  // 6. C2_EllEq_Brown: solving ell. eqs. inside the BH (Brown's method).
  //    it's C^2 continues. */
  Pset_default("extrapolate_inside_BH_method","linear");
  
  /* max allowed enthalpy L2 norm residual; if root finder of NS surface 
  // gets L2 norm residual below this number => it's been failed. */
  Pset_default("NS_enthalpy_allowed_residual","1e-4");
  
  /* weight for update of BH's spin */
  Pset_default("BH_spin_update_weight","0.0");
  
  /* relative tolerance(error) for BH's mass */
  Pset_default("BH_mass_tolerance","0.1");
  
  /* type of free data for BH
  // options:
  // 1. boosted_KerrSchild_metric
  // 2. conformally_flat_metric => this option is not complete.
  */
  Pset_default("BH_NS_free_data_metric","boosted_KerrSchild_metric");
  
  /* test elliptic equations converge [yes/no] */
  Pset_default("Elliptic_Convergence_Test","NO");
  
  /* update weight for Euler constant */
  Pset_default("NS_Euler_const_update_weight","1.0");
  
  /* root finder method for NS surface
  // options:
  // 1. Bisect_Single: using bisect method.
  // 2. Steepest_Descent: using steepest descent method. */
  Pset_default("NS_surface_root_finder","Bisect_Single");
 
}


