/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "bbn_headers.h"
#include "maths_calculus_lib.h"

#define add_field_and_get_field(name) \
 ADD_FIELD(name); DECLARE_FIELD(name);

void bbn_1st_2nd_derivatives_conformal_metric(Patch_T *const patch);
void bbn_1st_2nd_derivatives_conformal_metric(Patch_T *const patch)
{

  /* declaring: */
  DECLARE_FIELD(_gamma_D2D2)
  DECLARE_FIELD(_gamma_D0D2)
  DECLARE_FIELD(_gamma_D1D2)
  DECLARE_FIELD(_gamma_D0D1)
  DECLARE_FIELD(_gamma_D1D1)
  DECLARE_FIELD(_gamma_D0D0)
  add_field_and_get_field(_dgamma_D0D2D0)
  add_field_and_get_field(_dgamma_D1D2D2)
  add_field_and_get_field(_dgamma_D0D0D1)
  add_field_and_get_field(_dgamma_D1D2D0)
  add_field_and_get_field(_dgamma_D0D0D0)
  add_field_and_get_field(_dgamma_D0D2D1)
  add_field_and_get_field(_dgamma_D0D2D2)
  add_field_and_get_field(_dgamma_D2D2D1)
  add_field_and_get_field(_dgamma_D0D0D2)
  add_field_and_get_field(_dgamma_D2D2D0)
  add_field_and_get_field(_dgamma_D0D1D2)
  add_field_and_get_field(_dgamma_D1D2D1)
  add_field_and_get_field(_dgamma_D1D1D0)
  add_field_and_get_field(_dgamma_D2D2D2)
  add_field_and_get_field(_dgamma_D0D1D0)
  add_field_and_get_field(_dgamma_D0D1D1)
  add_field_and_get_field(_dgamma_D1D1D2)
  add_field_and_get_field(_dgamma_D1D1D1)
  add_field_and_get_field(_ddgamma_D2D2D2D2)
  add_field_and_get_field(_ddgamma_D0D0D2D2)
  add_field_and_get_field(_ddgamma_D0D1D1D1)
  add_field_and_get_field(_ddgamma_D1D1D0D2)
  add_field_and_get_field(_ddgamma_D1D2D2D2)
  add_field_and_get_field(_ddgamma_D1D1D1D0)
  add_field_and_get_field(_ddgamma_D1D2D1D2)
  add_field_and_get_field(_ddgamma_D1D2D0D2)
  add_field_and_get_field(_ddgamma_D2D2D1D1)
  add_field_and_get_field(_ddgamma_D1D1D2D0)
  add_field_and_get_field(_ddgamma_D1D2D0D1)
  add_field_and_get_field(_ddgamma_D0D1D1D0)
  add_field_and_get_field(_ddgamma_D0D2D2D1)
  add_field_and_get_field(_ddgamma_D0D0D2D1)
  add_field_and_get_field(_ddgamma_D0D0D1D2)
  add_field_and_get_field(_ddgamma_D0D1D0D0)
  add_field_and_get_field(_ddgamma_D0D1D1D2)
  add_field_and_get_field(_ddgamma_D0D2D0D0)
  add_field_and_get_field(_ddgamma_D0D2D0D1)
  add_field_and_get_field(_ddgamma_D0D1D2D1)
  add_field_and_get_field(_ddgamma_D0D1D0D2)
  add_field_and_get_field(_ddgamma_D0D1D2D2)
  add_field_and_get_field(_ddgamma_D1D2D0D0)
  add_field_and_get_field(_ddgamma_D0D2D1D2)
  add_field_and_get_field(_ddgamma_D0D1D2D0)
  add_field_and_get_field(_ddgamma_D1D2D2D1)
  add_field_and_get_field(_ddgamma_D2D2D2D0)
  add_field_and_get_field(_ddgamma_D1D1D0D1)
  add_field_and_get_field(_ddgamma_D1D2D1D1)
  add_field_and_get_field(_ddgamma_D0D0D2D0)
  add_field_and_get_field(_ddgamma_D1D1D2D2)
  add_field_and_get_field(_ddgamma_D2D2D2D1)
  add_field_and_get_field(_ddgamma_D2D2D0D2)
  add_field_and_get_field(_ddgamma_D0D2D1D1)
  add_field_and_get_field(_ddgamma_D2D2D0D1)
  add_field_and_get_field(_ddgamma_D1D1D1D1)
  add_field_and_get_field(_ddgamma_D0D2D0D2)
  add_field_and_get_field(_ddgamma_D1D1D1D2)
  add_field_and_get_field(_ddgamma_D2D2D1D2)
  add_field_and_get_field(_ddgamma_D0D2D2D2)
  add_field_and_get_field(_ddgamma_D1D2D1D0)
  add_field_and_get_field(_ddgamma_D0D0D0D1)
  add_field_and_get_field(_ddgamma_D0D1D0D1)
  add_field_and_get_field(_ddgamma_D0D2D2D0)
  add_field_and_get_field(_ddgamma_D1D1D0D0)
  add_field_and_get_field(_ddgamma_D0D0D1D1)
  add_field_and_get_field(_ddgamma_D1D2D2D0)
  add_field_and_get_field(_ddgamma_D0D2D1D0)
  add_field_and_get_field(_ddgamma_D1D1D2D1)
  add_field_and_get_field(_ddgamma_D0D0D1D0)
  add_field_and_get_field(_ddgamma_D2D2D0D0)
  add_field_and_get_field(_ddgamma_D0D0D0D0)
  add_field_and_get_field(_ddgamma_D0D0D0D2)
  add_field_and_get_field(_ddgamma_D2D2D1D0)



  _dgamma_D1D1D2->v=Partial_Derivative(_gamma_D1D1,"z");
  _dgamma_D2D2D0->v=Partial_Derivative(_gamma_D2D2,"x");
  _dgamma_D0D2D0->v=Partial_Derivative(_gamma_D0D2,"x");
  _dgamma_D1D1D1->v=Partial_Derivative(_gamma_D1D1,"y");
  _dgamma_D0D1D2->v=Partial_Derivative(_gamma_D0D1,"z");
  _dgamma_D0D1D1->v=Partial_Derivative(_gamma_D0D1,"y");
  _dgamma_D0D0D1->v=Partial_Derivative(_gamma_D0D0,"y");
  _dgamma_D0D1D0->v=Partial_Derivative(_gamma_D0D1,"x");
  _dgamma_D1D2D1->v=Partial_Derivative(_gamma_D1D2,"y");
  _dgamma_D2D2D1->v=Partial_Derivative(_gamma_D2D2,"y");
  _dgamma_D2D2D2->v=Partial_Derivative(_gamma_D2D2,"z");
  _dgamma_D0D0D2->v=Partial_Derivative(_gamma_D0D0,"z");
  _dgamma_D0D2D1->v=Partial_Derivative(_gamma_D0D2,"y");
  _dgamma_D1D2D0->v=Partial_Derivative(_gamma_D1D2,"x");
  _dgamma_D0D0D0->v=Partial_Derivative(_gamma_D0D0,"x");
  _dgamma_D0D2D2->v=Partial_Derivative(_gamma_D0D2,"z");
  _dgamma_D1D1D0->v=Partial_Derivative(_gamma_D1D1,"x");
  _dgamma_D1D2D2->v=Partial_Derivative(_gamma_D1D2,"z");
  _ddgamma_D0D0D0D2->v=Partial_Derivative(_dgamma_D0D0D0,"z");
  _ddgamma_D0D1D1D1->v=Partial_Derivative(_dgamma_D0D1D1,"y");
  _ddgamma_D2D2D2D2->v=Partial_Derivative(_dgamma_D2D2D2,"z");
  _ddgamma_D0D0D0D1->v=Partial_Derivative(_dgamma_D0D0D0,"y");
  _ddgamma_D0D1D0D0->v=Partial_Derivative(_dgamma_D0D1D0,"x");
  _ddgamma_D0D2D2D2->v=Partial_Derivative(_dgamma_D0D2D2,"z");
  _ddgamma_D0D2D1D2->v=Partial_Derivative(_dgamma_D0D2D1,"z");
  _ddgamma_D1D1D0D2->v=Partial_Derivative(_dgamma_D1D1D0,"z");
  _ddgamma_D2D2D0D0->v=Partial_Derivative(_dgamma_D2D2D0,"x");
  _ddgamma_D1D1D0D1->v=Partial_Derivative(_dgamma_D1D1D0,"y");
  _ddgamma_D1D2D2D0->v=Partial_Derivative(_dgamma_D1D2D2,"x");
  _ddgamma_D2D2D2D0->v=Partial_Derivative(_dgamma_D2D2D2,"x");
  _ddgamma_D0D0D1D2->v=Partial_Derivative(_dgamma_D0D0D1,"z");
  _ddgamma_D0D2D1D0->v=Partial_Derivative(_dgamma_D0D2D1,"x");
  _ddgamma_D0D1D2D2->v=Partial_Derivative(_dgamma_D0D1D2,"z");
  _ddgamma_D0D0D2D2->v=Partial_Derivative(_dgamma_D0D0D2,"z");
  _ddgamma_D0D0D1D0->v=Partial_Derivative(_dgamma_D0D0D1,"x");
  _ddgamma_D1D2D1D0->v=Partial_Derivative(_dgamma_D1D2D1,"x");
  _ddgamma_D1D1D2D0->v=Partial_Derivative(_dgamma_D1D1D2,"x");
  _ddgamma_D0D1D1D0->v=Partial_Derivative(_dgamma_D0D1D1,"x");
  _ddgamma_D2D2D1D1->v=Partial_Derivative(_dgamma_D2D2D1,"y");
  _ddgamma_D2D2D0D2->v=Partial_Derivative(_dgamma_D2D2D0,"z");
  _ddgamma_D2D2D1D2->v=Partial_Derivative(_dgamma_D2D2D1,"z");
  _ddgamma_D0D1D0D1->v=Partial_Derivative(_dgamma_D0D1D0,"y");
  _ddgamma_D0D1D2D1->v=Partial_Derivative(_dgamma_D0D1D2,"y");
  _ddgamma_D0D0D2D0->v=Partial_Derivative(_dgamma_D0D0D2,"x");
  _ddgamma_D1D1D2D2->v=Partial_Derivative(_dgamma_D1D1D2,"z");
  _ddgamma_D1D2D2D2->v=Partial_Derivative(_dgamma_D1D2D2,"z");
  _ddgamma_D2D2D0D1->v=Partial_Derivative(_dgamma_D2D2D0,"y");
  _ddgamma_D0D2D2D1->v=Partial_Derivative(_dgamma_D0D2D2,"y");
  _ddgamma_D0D0D1D1->v=Partial_Derivative(_dgamma_D0D0D1,"y");
  _ddgamma_D1D1D1D0->v=Partial_Derivative(_dgamma_D1D1D1,"x");
  _ddgamma_D1D1D1D1->v=Partial_Derivative(_dgamma_D1D1D1,"y");
  _ddgamma_D0D2D2D0->v=Partial_Derivative(_dgamma_D0D2D2,"x");
  _ddgamma_D2D2D1D0->v=Partial_Derivative(_dgamma_D2D2D1,"x");
  _ddgamma_D1D2D1D1->v=Partial_Derivative(_dgamma_D1D2D1,"y");
  _ddgamma_D0D1D2D0->v=Partial_Derivative(_dgamma_D0D1D2,"x");
  _ddgamma_D0D2D0D2->v=Partial_Derivative(_dgamma_D0D2D0,"z");
  _ddgamma_D1D2D2D1->v=Partial_Derivative(_dgamma_D1D2D2,"y");
  _ddgamma_D0D2D1D1->v=Partial_Derivative(_dgamma_D0D2D1,"y");
  _ddgamma_D0D0D0D0->v=Partial_Derivative(_dgamma_D0D0D0,"x");
  _ddgamma_D0D2D0D0->v=Partial_Derivative(_dgamma_D0D2D0,"x");
  _ddgamma_D0D1D1D2->v=Partial_Derivative(_dgamma_D0D1D1,"z");
  _ddgamma_D1D2D0D0->v=Partial_Derivative(_dgamma_D1D2D0,"x");
  _ddgamma_D0D2D0D1->v=Partial_Derivative(_dgamma_D0D2D0,"y");
  _ddgamma_D0D0D2D1->v=Partial_Derivative(_dgamma_D0D0D2,"y");
  _ddgamma_D1D1D2D1->v=Partial_Derivative(_dgamma_D1D1D2,"y");
  _ddgamma_D1D1D0D0->v=Partial_Derivative(_dgamma_D1D1D0,"x");
  _ddgamma_D1D1D1D2->v=Partial_Derivative(_dgamma_D1D1D1,"z");
  _ddgamma_D1D2D0D1->v=Partial_Derivative(_dgamma_D1D2D0,"y");
  _ddgamma_D1D2D0D2->v=Partial_Derivative(_dgamma_D1D2D0,"z");
  _ddgamma_D0D1D0D2->v=Partial_Derivative(_dgamma_D0D1D0,"z");
  _ddgamma_D1D2D1D2->v=Partial_Derivative(_dgamma_D1D2D1,"z");
  _ddgamma_D2D2D2D1->v=Partial_Derivative(_dgamma_D2D2D2,"y");

}