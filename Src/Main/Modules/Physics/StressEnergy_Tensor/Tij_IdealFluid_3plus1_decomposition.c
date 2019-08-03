/*
// Alireza Rashti
// August 2019
*/

#include "Tij_IdealFluid_3plus1_decomposition.h"

/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid. thus, 
// Tij_IF means stress enerfy of ideal fluid. 
// given all of the fields needed it builds 
// the first component of fluid four velocity, i.e. 
// u_mu = (u_U0,u_U1,u_U2,u_U3). */
void Tij_IF_u0(Patch_T *const patch)
{
}



/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid. thus, 
// Tij_IF means stress enerfy of ideal fluid.
// given all of the fields needed for J = -gamma(i,-mu)*T(mu,nu)*n(-nu) 
// it builds "momentum current * psi^6", where psi is conformal factor,
// and puts it to _J_U?.
// note: if patch does not contain fluid, it does nothing. */
void Tij_IF_build_psi6J_Ui(Patch_T *const patch)
{
}




/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid. thus, 
// Tij_IF means stress enerfy of ideal fluid.
// given all of the fields needed for E = T(mu,nu)*n(-mu)*n(-nu)
// it builds "total energy density * psi^6", where psi is conformal factor,
// and puts it to _E.
// note: if patch does not contain fluid, it does nothing. */
void Tij_IF_build_psi6E(Patch_T *const patch)
{
}


/* Note: Tij stands for Stress Energy tensor.
// IF stands for ideal fluid. thus, 
// Tij_IF means stress enerfy of ideal fluid.
// given all of the fields needed for 
// S = T(mu,nu)*gamma(i,j)*gamma(-i,-mu)*gamma(-j,-nu)
// it builds "total energy density * psi^6", 
// where psi is conformal factor, and puts it to _S. 
// note: if patch does not contain fluid, it does nothing. */
void Tij_IF_build_psi6S(Patch_T *const patch)
{
}
