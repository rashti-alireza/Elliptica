/*
// Alireza Rashti
// August 2019
*/

#include "Tij_NS_ifluid_xcts_gConf_update.h"


/* updating all matter related and their derivatives. */
void Tij_NS_idealfluid_XCTS_gConf_update(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid  = mygrid(phys,Ftype("NS"));
  Patch_T *patch      = 0;
  const int  neat = strstr_i(Gets("enthalpy_neat"),"yes");
  const double W  = (neat == 0 ? Getd("enthalpy_update_weight") : 0.);
  const double Euler_const = Getd("Euler_equation_constant");
  const double NS_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  double X[3] = {0};
                                 
  printf(Pretty0"update weight  = %e\n"
         Pretty0"neat it?       = %d\n"
         Pretty0"Euler constant = %e\n",
         W,neat,Euler_const);
  fflush(stdout);
  
  EoS_T *eos     = init_EoS(phys);
  FOR_ALL_p(grid->np)
  {
    patch = grid->patch[p];
    
    /* update values */
    if (!EQL(W,0.))
      RELAX_UPDATE_FUNC(Tij_NS_IF_XCTS_gConf_enthalpy(patch,Euler_const),
                        patch,enthalpy,W);
    
    if (neat)
      Tij_NS_neat_enthalpy(patch);
    
    Tij_NS_eos_update_rho0(patch,eos);
    Tij_NS_IF_XCTS_gConf_u0(patch);
    
    /* update derivatives */
    if (1)
    {
      dField_di(denthalpy_D0);
      dField_di(denthalpy_D1);
      dField_di(denthalpy_D2);
      dField_di(du0_D0);
      dField_di(du0_D1);
      dField_di(du0_D2);
      /* NOTE: rho0 for pwp eos is C^0 continuous so it is prone to Gibbs phenomenon
      // in particular where h ~ 1, i.e., close to NS's surface. hence, calculating 
      // drho0 using analytical values and thus chain rule is preferable. */
      READ_v(enthalpy)
      READ_v(denthalpy_D0)
      READ_v(denthalpy_D1)
      READ_v(denthalpy_D2)
      
      REALLOC_v_WRITE_v(drho0_D0)
      REALLOC_v_WRITE_v(drho0_D1)
      REALLOC_v_WRITE_v(drho0_D2)
      
      FOR_ALL_ijk
      {
        eos->h          = enthalpy[ijk];
        double drho0_dh = eos->drho0_dh(eos);
        drho0_D0[ijk] = drho0_dh*denthalpy_D0[ijk];
        drho0_D1[ijk] = drho0_dh*denthalpy_D1[ijk];
        drho0_D2[ijk] = drho0_dh*denthalpy_D2[ijk];
      }
    }
    else/* the following not using chain rule for rho0 */
    {
      Tij_NS_IF_XCTS_gConf_derives(patch);
    }
    
    /* sources */
    Tij_NS_IF_XCTS_gConf_psi6J_Ui(patch,eos);
    Tij_NS_IF_XCTS_gConf_psi6E(patch,eos);
    Tij_NS_IF_XCTS_gConf_psi6S(patch,eos);
  }
  
  /* update central quantities */
  Interpolation_T *interp_s = init_interpolation();
  patch = x_in_which_patch(NS_center,grid->patch,grid->np);
  assert(patch);
  assert(X_of_x(X,NS_center,patch));
  interp_s->field = patch->fields[Ind("enthalpy")];
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  eos->h = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  Setd("rho0_center",eos->rest_mass_density(eos));
  Setd("pressure_center",eos->pressure(eos));
  Setd("energy_density_center",eos->energy_density(eos));
  
  free_EoS(eos);
  
  FUNC_TOC
}

/* updating rho0 using eos */
void Tij_NS_eos_update_rho0(Patch_T *const patch,EoS_T *const eos)
{
  READ_v(enthalpy)
  REALLOC_v_WRITE_v(rho0)
  const Uint nn = patch->nn;
  Uint ijk;

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
}

/* after finding new NS surface, root finder might find h ~ 1
// so put it h = 1, to prevent nan in matter fields  */
void Tij_NS_neat_enthalpy(Patch_T *const patch)
{
  WRITE_v(enthalpy)
  
  const Uint nn = patch->nn;
  Uint ijk;
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    if (LSSEQL(enthalpy[ijk],1))
      enthalpy[ijk] = 1;
  }

}

