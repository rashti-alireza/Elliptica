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
  
}
