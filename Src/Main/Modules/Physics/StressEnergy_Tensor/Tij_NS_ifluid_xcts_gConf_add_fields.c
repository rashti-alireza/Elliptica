/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "Tij_header.h"


void Tij_NS_idealfluid_XCTS_gConf_add_fields(Grid_T *const grid)
{
 Uint p;
 FOR_ALL_PATCHES(p,grid)
 {
 Patch_T *patch = grid->patch[p];


  /* declaring: */
  ADD_FIELD(enthalpy)
  ADD_FIELD(denthalpy_D2)
  ADD_FIELD(denthalpy_D0)
  ADD_FIELD(denthalpy_D1)
  ADD_FIELD(u0)
  ADD_FIELD(du0_D0)
  ADD_FIELD(du0_D1)
  ADD_FIELD(du0_D2)
  ADD_FIELD(rho0)
  ADD_FIELD(drho0_D2)
  ADD_FIELD(drho0_D0)
  ADD_FIELD(drho0_D1)
  ADD_AND_ALLOC_FIELD(EConf)
  ADD_AND_ALLOC_FIELD(JConf_U0)
  ADD_AND_ALLOC_FIELD(JConf_U1)
  ADD_AND_ALLOC_FIELD(JConf_U2)
  ADD_AND_ALLOC_FIELD(SConf)


 }
}