/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "bbn_headers.h"
#include "maths_equation_solvings_lib.h"

void *bbn_bhf_eq_Brown(void *vp1,void *vp2);
void *bbn_bhf_eq_Brown(void *vp1,void *vp2)
{
  DDM_SCHUR_EQ_DECLARE
  unsigned ijk;/* node index */
  const Solving_Man_T *const sol = patch->solving_man;
  const char *const fld_name     = sol->field_name[sol->cf];
  char eq_fld_name[100] = {'\0'};


  /* declaring: */


sprintf(eq_fld_name,"src_%s",fld_name);
const double *const B = patch->pool[Ind(eq_fld_name)]->v;
sprintf(eq_fld_name,"dd%s_D0D0",fld_name);
const double *const ddB_xx = patch->pool[Ind(eq_fld_name)]->v;
sprintf(eq_fld_name,"dd%s_D1D1",fld_name);
const double *const ddB_yy = patch->pool[Ind(eq_fld_name)]->v;
sprintf(eq_fld_name,"dd%s_D2D2",fld_name);
const double *const ddB_zz = patch->pool[Ind(eq_fld_name)]->v;
  DDM_SCHUR_EQ_OPEN

  double F_eq = 
-B[ijk] + ddB_xx[ijk] + ddB_yy[ijk] + ddB_zz[ijk];

  F[n] = F_eq;

  DDM_SCHUR_EQ_CLOSE

  return 0;
}