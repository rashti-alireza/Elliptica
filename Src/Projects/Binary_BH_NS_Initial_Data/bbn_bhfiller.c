/*
// Alireza Rashti
// September 2020
*/
/* filling the inside of BH for evolution purposes */

#include "bbn_bhfiller.h"

/* these fields to be extrapolated  */
static const char *const fields_name[] = {
  "psi","eta","K",
  "Beta_U0","Beta_U1","Beta_U2",
  "_gamma_D2D2","_gamma_D0D2",
  "_gamma_D0D0","_gamma_D0D1",
  "_gamma_D1D2","_gamma_D1D1",0};

/* initialize the bhfiller struct */
struct BHFiller_S* 
bhf_init
  (
  Grid_T *const grid/* the whole grid */,
  const char *const method/* the method to be used for extrapolating */
  )
{
  struct BHFiller_S *const bhf = calloc(1,sizeof(*bhf));IsNull(bhf);
  const double EPS   = 1;
  const double Ma    = Pgetd("BH_irreducible_mass");
  const unsigned npi = 7;/* number of patches inside BH */
  const unsigned npo = 6;/* number of patches outside BH */
  const double fr0_Beta_U0     = 0;
  const double fr0_Beta_U1     = 0;
  const double fr0_Beta_U2     = 0;
  const double fr0_gamma_D0D0  = 1;
  const double fr0_gamma_D0D1  = 0;
  const double fr0_gamma_D0D2  = 0;
  const double fr0_gamma_D1D1  = 1;
  const double fr0_gamma_D1D2  = 0;
  const double fr0_gamma_D2D2  = 1;
  const double fr0_K           = 0;
  const double fr0_alpha       = 0.1;
  const double fr0_psi = 2+Ma/(2*EPS);
  const double fr0_eta = fr0_alpha*fr0_psi;
  const char *s = 0;
  unsigned f,nf,i,j,p;
  
  bhf->npo = npo;
  bhf->npi = npi;
  nf = 0;/* number of fields */
  while(fields_name[nf]) ++nf;
  bhf->nf = nf;
  /* using Cheb_Tn and Ylm as the bases 
  // and extrapolate demanding C2 continuity */
  if (strcmp_i(method,"TnYlm_C2"))
  {
    /* if par doesn't exist */
    int lmax_par = PgetiEZ("bbn_bhfiller_lmax");
    if (lmax_par == INT_MAX)
      lmax_par = 10;
    const unsigned lmax   = (unsigned)lmax_par;
    const unsigned Ntheta = 2*lmax+1;
    const unsigned Nphi   = 2*lmax+1;
    const unsigned N      = Ntheta*Nphi;
    
    /* initialize tables */
    init_Legendre_root_function();
  
    bhf->lmax   = lmax;
    bhf->Ntheta = Ntheta;
    bhf->Nphi   = Nphi;
    bhf->fld    = calloc(nf,sizeof(*bhf->fld));IsNull(bhf->fld);
    sprintf(bhf->method,"%s",method);
    for (f = 0; f < nf; ++f)
    {
      bhf->fld[f] = calloc(1,sizeof(*bhf->fld[f]));IsNull(bhf->fld[f]);
      /* names of fields and its derivatives */
      sprintf(bhf->fld[f]->f,"%s",fields_name[f]);
      if (fields_name[f][0] == '_')/* => _dgamma */
      {
        s = bhf->fld[f]->f+1;
        /* if it is indexed */
        if (regex_search(".+_(U|D)[[:digit:]]",s))
        {
          for (i = 0; i < 3; ++i)
          {
            sprintf(bhf->fld[f]->df[i],"_d%sD%u",s,i);
            for (j = i; j < 3; ++j)
              sprintf(bhf->fld[f]->ddf[IJsymm3(i,j)],"_dd%sD%uD%u",s,i,j);
          }
        }
        else/* not indexed */
        {
          for (i = 0; i < 3; ++i)
          {
            sprintf(bhf->fld[f]->df[0],"_d%s_D%u",s,i);
            for (j = i; j < 3; ++j)
              sprintf(bhf->fld[f]->ddf[IJsymm3(i,j)],"_dd%s_D%uD%u",s,i,j);
          }
        }
      }
      else/* if no _ at the beginning */
      {
        s = fields_name[f];
        /* if it is indexed */
        if (regex_search(".+_(U|D)[[:digit:]]",s))
        {
          for (i = 0; i < 3; ++i)
          {
            sprintf(bhf->fld[f]->df[i],"d%sD%u",s,i);
            for (j = i; j < 3; ++j)
              sprintf(bhf->fld[f]->ddf[IJsymm3(i,j)],"dd%sD%uD%u",s,i,j);
          }
        }
        else/* not indexed */
        {
          for (i = 0; i < 3; ++i)
          {
            sprintf(bhf->fld[f]->df[i],"d%s_D%u",s,i);
            for (j = i; j < 3; ++j)
              sprintf(bhf->fld[f]->ddf[IJsymm3(i,j)],"dd%s_D%uD%u",s,i,j);
          }
        }
      }
      
      /* alloc ChebTn_coeffs */
      for (i = 0 ; i < 4; ++i)
      {
        bhf->fld[f]->ChebTn_coeffs[i]  = alloc_double(N);
        bhf->fld[f]->realYlm_coeffs[i] = alloc_ClmYlm(lmax);
        bhf->fld[f]->imagYlm_coeffs[i] = alloc_ClmYlm(lmax);
      }
    }/* for (f = 0; f < nf; ++f) */
  }/* if (strcmp_i(method,"TnYlm_C2")) */
  else
    Error0(NO_OPTION);
  
  /* quick test for names */
  if (0)
  {
    /* show contents */
    for (f = 0; f < nf; ++f)
    {
      pr_line();
      printf("fld[%u] = %s\n",f,bhf->fld[f]->f);
      printf("d[%s] = (%s,%s,%s)\n",
          bhf->fld[f]->f,bhf->fld[f]->df[0],
          bhf->fld[f]->df[1],bhf->fld[f]->df[2]);
      for (i = 0; i < 3; ++i)
      {
        for (j = i; j < 3; ++j)
        {
          printf("dd[%s](%u,%u) = %s\n",
            bhf->fld[f]->f,i,j,bhf->fld[f]->ddf[IJsymm3(i,j)]);
        }
      }
    }
  }

  /* set values of field at r=0 */
  for (f = 0; f < nf; ++f)
  {
    if (strcmp_i(fields_name[f],"psi"))
    {
      bhf->fld[f]->f_r0 = fr0_psi;
    }
    else if (strcmp_i(fields_name[f],"eta"))
    {
      bhf->fld[f]->f_r0 = fr0_eta;
    }
    else if (strcmp_i(fields_name[f],"K"))
    {
      bhf->fld[f]->f_r0 = fr0_K;
    }
    else if (strcmp_i(fields_name[f],"Beta_U0"))
    {
      bhf->fld[f]->f_r0 = fr0_Beta_U0;
    }
    else if (strcmp_i(fields_name[f],"Beta_U1"))
    {
      bhf->fld[f]->f_r0 = fr0_Beta_U1;
    }
    else if (strcmp_i(fields_name[f],"Beta_U2"))
    {
      bhf->fld[f]->f_r0 = fr0_Beta_U2;
    }
    else if (strcmp_i(fields_name[f],"_gamma_D2D2"))
    {
      bhf->fld[f]->f_r0 = fr0_gamma_D2D2;
    }
    else if (strcmp_i(fields_name[f],"_gamma_D0D2"))
    {
      bhf->fld[f]->f_r0 = fr0_gamma_D0D2;
    }
    else if (strcmp_i(fields_name[f],"_gamma_D0D0"))
    {
      bhf->fld[f]->f_r0 = fr0_gamma_D0D0;
    }
    else if (strcmp_i(fields_name[f],"_gamma_D0D1"))
    {
      bhf->fld[f]->f_r0 = fr0_gamma_D0D1;
    }
    else if (strcmp_i(fields_name[f],"_gamma_D1D2"))
    {
      bhf->fld[f]->f_r0 = fr0_gamma_D1D2;
    }
    else if (strcmp_i(fields_name[f],"_gamma_D1D1"))
    {
      bhf->fld[f]->f_r0 = fr0_gamma_D1D1;
    }
    else
      Error0(NO_OPTION);
    
  }/* for (f = 0; f < nf ++f) */
  
  /* patches outside the BH */
  bhf->patches_outBH = calloc(npo,sizeof(*bhf->patches_outBH));
  IsNull(bhf->patches_outBH);
  /* patches inside the BH */
  bhf->patches_inBH = calloc(npi,sizeof(*bhf->patches_inBH));
  IsNull(bhf->patches_inBH);
  
  i = j = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if (IsItHorizonPatch(patch))
      bhf->patches_outBH[i++] = patch;
    else if (IsItInsideBHPatch(patch))
     bhf->patches_inBH[j++] = patch;
  }
  assert(i == npo);
  assert(j == npi);
  
  return bhf;
}

/* free bhfiller struct */
static void bhf_free(struct BHFiller_S *const bhf)
{
  unsigned i,f;
  
  if(!bhf)
    return;
  
  for (f = 0; f < bhf->nf; ++f)
  {
    for (i = 0; i < 4; ++i)
    {
      _free(bhf->fld[f]->ChebTn_coeffs[i]);
      _free(bhf->fld[f]->realYlm_coeffs[i]);
      _free(bhf->fld[f]->imagYlm_coeffs[i]);
    }
  }
  free_2d_mem(bhf->fld,bhf->nf);
  _free(bhf->patches_outBH);
  _free(bhf->patches_inBH);
  free(bhf);
}

/* ->: EXIT_SUCCESS if succeeds. 
// extrapolating inside the BH */
int 
bbn_bhfiller
  (
  Grid_T *const grid/* the whole grid */,
  const char *const method/* the method to be used for extrapolating */
  )
{
  /* check if it is perfect sphere */
  if (!Pcmps("BH_R_type","PerfectSphere"))
    Error0(NO_OPTION);

  struct BHFiller_S *const bhf = bhf_init(grid,method);
  const double EPS   = 1E-12;
  const unsigned npo = bhf->npo;
  const unsigned npi = bhf->npi;
  const unsigned nf  = bhf->nf;/* numebr of fields */
  const unsigned lmax   = bhf->lmax;
  const unsigned Ntheta = bhf->Ntheta;
  const unsigned Nphi   = bhf->Nphi;
  const double r_fill = Pgetd("r_excision");
  const double r_fill3= pow(r_fill,3);
  unsigned p,fld;
  
  /* update all coeffs to avoid race condition */
  printf("--> Updating coefficients ...\n");
  fflush(stdout);
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = bhf->patches_outBH[p];
    unsigned f;

    bbn_1st_2nd_derivatives_conformal_metric(patch);
    bbn_add_and_take_2nd_derivatives_K(patch);
    Field_T *R1_f  = patch->CoordSysInfo->CubedSphericalCoord->R1_f;
    Field_T *R2_f  = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
    if (R1_f)
      make_coeffs_2d(R1_f,0,1);/* X and Y direction */
    if (R2_f)
      make_coeffs_2d(R2_f,0,1);/* X and Y direction */

    /* make coeffs for all fields inside this patch */
    for (f = 0; f < patch->nfld; ++f)
    {
      if (patch->pool[f]->v      &&
          patch->pool[f] != R1_f && 
          patch->pool[f] != R2_f    )
        make_coeffs_3d(patch->pool[f]);
    }
  }
  
  /* populating f, df/dr, d^2f/dr^2 at each (th,ph) points */
  printf("--> Populating extrapolation coefficients ...\n");
  fflush(stdout);
  OpenMP_1d_Pragma(omp parallel for)
  for (fld = 0; fld < nf; ++fld)
  {
    unsigned i,j,_i,_j;
    for (i = 0; i < Ntheta; ++i)
    {
      double theta = acos(-Legendre_root_function(i,Ntheta));
      for (j = 0; j < Nphi; ++j)
      {
        double phi = j*2*M_PI/Nphi;
        Patch_T *patch    = 0;
        double KD[2]      = {0,1};
        double df_dx[3]   = {0};
        double ddf_ddx[6] = {0};
        double ddfddr = 0,dfdr = 0,f_r1 = 0,f_r0 = 0;
        double a[4] = {0};
        double _ddfddr[3] = {0,0,0};
        unsigned ij = IJ(i,j,Nphi);
        unsigned d1,d2;/* derivative */
        double X[3],x[3],_x[3],N[3];
        /* find patch for the given theta and phi */
        find_XYZ_and_patch_of_theta_phi_BH_CS(X,&patch,theta,phi,grid);

        /* r = r_fill(sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^) */
        _x[0] = r_fill*sin(theta)*cos(phi);
        _x[1] = r_fill*sin(theta)*sin(phi);
        _x[2] = r_fill*cos(theta);
        x[0]  = _x[0] + patch->c[0];
        x[1]  = _x[1] + patch->c[1];
        x[2]  = _x[2] + patch->c[2];
        assert(X_of_x(X,x,patch));
        
        /* normal vector */
        N[0]  = sin(theta)*cos(phi);
        N[1]  = sin(theta)*sin(phi);
        N[2]  = cos(theta);
        
        /* interpolate on the surface */
        Interpolation_T *interp_s = init_interpolation();
        interp_s->XY_dir_flag = 1;
        interp_s->X = X[0];
        interp_s->Y = X[1];
        interp_s->K = 0;
        /* f value at r = r0 and r = 0*/
        interp_s->field = patch->pool[Ind(bhf->fld[fld]->f)];
        plan_interpolation(interp_s);
        f_r1 = execute_interpolation(interp_s);
        f_r0 = bhf->fld[fld]->f_r0;
        
        /* df/dx value */
        for (d1 = 0; d1 < 3; d1++)
        {
          interp_s->field = patch->pool[Ind(bhf->fld[fld]->df[d1])];
          plan_interpolation(interp_s);
          df_dx[d1] = execute_interpolation(interp_s);
        }
        /* d^2f/dx^2 value */
        for (d1 = 0; d1 < 3; d1++)
        {
          for (d2 = d1; d2 < 3; d2++)
          {
            interp_s->field = 
              patch->pool[Ind(bhf->fld[fld]->ddf[IJsymm3(d1,d2)])];
            plan_interpolation(interp_s);
            ddf_ddx[IJsymm3(d1,d2)] = execute_interpolation(interp_s);
          }
        }
        free_interpolation(interp_s);
        
        /* df/dr */
        dfdr = (N[0]*df_dx[0]+N[1]*df_dx[1]+N[2]*df_dx[2]);
        
        ddfddr = 0;
        _ddfddr[0] = _ddfddr[1] = _ddfddr[2] = 0;
        /* d^2f/dr^2 */
        for (_i = 0; _i < 3; ++_i)
        {
          for (_j = 0; _j < 3; ++_j)
          {
            _ddfddr[_i] += (KD[_i==_j]/r_fill - _x[_i]*_x[_j]/r_fill3)*df_dx[_j];
            _ddfddr[_i] += N[_j]*ddf_ddx[IJsymm3(_i,_j)];
          }
          ddfddr += _ddfddr[_i]*N[_i]; \
        }
        
        a[0] = (10*f_r0 + 22*f_r1 + 2*ddfddr*r_fill - 6*dfdr*r_fill)/32.;
        a[1] = (-15*f_r0 + 15*f_r1 - ddfddr*r_fill + dfdr*r_fill)/32.;
        a[2] = (6*f_r0 - 6*f_r1 - 2*ddfddr*r_fill + 6*dfdr*r_fill)/32.;
        a[3] = (-f_r0 + f_r1 + ddfddr*r_fill - dfdr*r_fill)/32.;
        
        for (_i = 0; _i < 4; _i++)
          bhf->fld[fld]->ChebTn_coeffs[_i][ij] = a[_i];
        
      }
    }/* for (i = 0; i < Ntheta; ++i) */
    /* now populate the Ylm coeffs */
    for ( i = 0 ; i < 4; ++i)
    {
      double *rC = bhf->fld[fld]->realYlm_coeffs[i];
      double *iC = bhf->fld[fld]->imagYlm_coeffs[i];
      double *v  = bhf->fld[fld]->ChebTn_coeffs[i];
      get_Ylm_coeffs(rC,iC,v,Ntheta,Nphi,lmax);
    }
  }/* for (fld = 0; fld < nf ++fld) */
  
  /* now fill the BH */
  printf("--> Fill the holes ...\n");
  fflush(stdout);
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npi; p++)
  {
    Patch_T *patch = bhf->patches_inBH[p];
    unsigned nn = patch->nn;
    double theta = 0,phi = 0, t = 0;
    unsigned ijk;
    unsigned f,i;

    bbn_add_fields_in_patch(patch);
    
    /* loop over all fields to be extrapolated */
    for (f = 0; f < nf; ++f)
    {
      Field_T *u = patch->pool[Ind(bhf->fld[f]->f)];
      empty_field(u);
      u->v      = alloc_double(patch->nn);
      double *v = u->v;
      
      for (ijk = 0; ijk < nn; ++ijk)
      {
        DEF_RELATIVE_x
        DEF_RELATIVE_y
        DEF_RELATIVE_z
        DEF_RELATIVE_r
        if (!EQL(r,0))
        {
          theta = acos(z/r);
        }
        else
        {
          r     = EPS;
          theta = 0;
        }
        t   = 2*r/r_fill-1;
        phi = arctan(y,x);
        for (i = 0; i < 4; ++i)
        {
          double *rC = bhf->fld[f]->realYlm_coeffs[i];
          double *iC = bhf->fld[f]->imagYlm_coeffs[i];
          v[ijk] += 
            interpolation_Ylm(rC,iC,lmax,theta,phi)*Cheb_Tn((int)i,t);
        }
      }/* for (ijk = 0; ijk < nn; ++ijk) */
    }/* for (f = 0; f < nf ++f) */
    
    /* _gammaI */
    REALLOC_v_WRITE_v(_gammaI_U2U2)
    REALLOC_v_WRITE_v(_gammaI_U0U2)
    REALLOC_v_WRITE_v(_gammaI_U0U0)
    REALLOC_v_WRITE_v(_gammaI_U0U1)
    REALLOC_v_WRITE_v(_gammaI_U1U2)
    REALLOC_v_WRITE_v(_gammaI_U1U1)
    READ_v(_gamma_D2D2)
    READ_v(_gamma_D0D2)
    READ_v(_gamma_D0D0)
    READ_v(_gamma_D0D1)
    READ_v(_gamma_D1D2)
    READ_v(_gamma_D1D1)
    READ_v(Beta_U0)
    READ_v(Beta_U1)
    READ_v(Beta_U2)
    REALLOC_v_WRITE_v(B0_U0)
    REALLOC_v_WRITE_v(B0_U1)
    REALLOC_v_WRITE_v(B0_U2)
    bbn_update_B1_U012(patch);
    READ_v(B1_U0)
    READ_v(B1_U1)
    READ_v(B1_U2)
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      B0_U0[ijk] = Beta_U0[ijk]-B1_U0[ijk];
      B0_U1[ijk] = Beta_U1[ijk]-B1_U1[ijk];
      B0_U2[ijk] = Beta_U2[ijk]-B1_U2[ijk];
      /* _gammaI =  _gamma inverse */
      COMPUTE_gammaI(_gamma_D0D0[ijk],_gamma_D0D1[ijk],_gamma_D0D2[ijk],
                     _gamma_D0D1[ijk],_gamma_D1D1[ijk],_gamma_D1D2[ijk],
                     _gamma_D0D2[ijk],_gamma_D1D2[ijk],_gamma_D2D2[ijk])
                     
      /* quick test check _gamma * _gammaI = delta */
      if (0)
      {
          double delta_U0D0 = 
        _gammaI_U0U0[ijk]*_gamma_D0D0[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D0D2[ijk];

          double delta_U0D1 = 
        _gammaI_U0U0[ijk]*_gamma_D0D1[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D1D2[ijk];

          double delta_U0D2 = 
        _gammaI_U0U0[ijk]*_gamma_D0D2[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U0U2[ijk]*_gamma_D2D2[ijk];

          double delta_U1D2 = 
        _gammaI_U0U1[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U1U2[ijk]*_gamma_D2D2[ijk];

          double delta_U1D0 = 
        _gammaI_U0U1[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D0D2[ijk];

         double delta_U1D1 = 
        _gammaI_U0U1[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D1D2[ijk];

          double delta_U2D2 = 
        _gammaI_U0U2[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U2U2[ijk]*_gamma_D2D2[ijk];

          double delta_U2D0 = 
        _gammaI_U0U2[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D0D2[ijk];

          double delta_U2D1 = 
        _gammaI_U0U2[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D1D2[ijk];

        if(!EQL(delta_U1D1,1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D1,0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D2,0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U1D2,0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D0,1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D1,0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D2,1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D0,0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U1D0,0))  Error0("_gammaI is not correct!\n");
      }
    }/* for (ijk = 0; ijk < nn; ++ijk) */
    /* update derivatives */
    bbn_update_derivative_Beta_U0(patch);
    bbn_update_derivative_Beta_U1(patch);
    bbn_update_derivative_Beta_U2(patch);
    bbn_update_derivative_psi(patch);
    bbn_update_derivative_eta(patch);
    
    /* for K_{ij} inside BH patches */
    bbn_1st_derivatives_conformal_metric(patch);
    bbn_free_data_Gamma_patch(patch);
    bbn_rm_1st_derivatives_conformal_metric(patch);
  }/* for (p = 0; p < npi; p++) */
  
  /* free */
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_up",grid));
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_down",grid));
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_left",grid));
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_right",grid));
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_back",grid));
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_front",grid));
  
  bhf_free(bhf);
  return EXIT_SUCCESS;
}


/* given theta, phi and knowing the fact that they are on BH surface, 
// it finds the corresponding patch and X,Y,Z coordinate. */
static void find_XYZ_and_patch_of_theta_phi_BH_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid)
{
  const double tan_phi    = tan(phi);
  const double cos_theta  = cos(theta);
  const double tan_phi2   = Pow2(tan_phi);
  const double cos_theta2 = Pow2(cos_theta);
  Flag_T found_flg = NO;
  unsigned p;
  
  X[2] = 0;/* since we are on BH surface from BH surrounding side */
  
  /* check all of BH patches in which (x,y,z) and 
  // (X,Y,Z) and (theta,phi) are consistent */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItHorizonPatch(patch))
      continue;

    Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
    const double *c = patch->c;
    double a = 0, b = 0;
    double a_sign = 0,b_sign = 0,c_sign = 0;
    double x[3],phi2,theta2,r;
    
    /* we know that theta = 0 or Pi occures only at UP or DOWN patches
    // so don't bother to follow algorithm all the way down.
    // furthermore, this prevent 0 in denominator of unrelated patches. */
    if (EQL(theta,0) || EQL(theta,M_PI))
    {
      if (side == LEFT || side == RIGHT || 
          side == BACK || side == FRONT   )
        continue;
    }
    
    /* first calculate the magnetitude of a and b 
    // which are related to X[0] and X[1] with a sign */
    switch (side)
    {
      case UP:
        a = Sqrt((1 - cos_theta2)/(cos_theta2 + cos_theta2*tan_phi2));
        b = tan_phi*Sqrt((1 - cos_theta2)/(cos_theta2*(1 + tan_phi2)));
      break;
      case DOWN:
        b = Sqrt((1 - cos_theta2)/(cos_theta2 + cos_theta2*tan_phi2));
        a = tan_phi*Sqrt((1 - cos_theta2)/(cos_theta2*(1 + tan_phi2)));
      break;
      case LEFT:
        a = 1/tan_phi;
        b = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/((1. - cos_theta2)*tan_phi2));
      break;
      case RIGHT:
        b = 1/tan_phi;
        a = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/((1. - cos_theta2)*tan_phi2));
      break;
      case BACK:
        a = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/(1 - cos_theta2));
        b = tan_phi;
      break;
      case FRONT:
        b = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/(1 - cos_theta2));
        a = tan_phi;
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* having found the magnitude of a and b, we need to find out the sign of them.
    // this is done by paying attention to side, signum(cos_theta) and range of tanphi */
    switch (side)
    {
      case UP:
        arctan_argument_signum(&b_sign,&a_sign,phi);
      break;
      case DOWN:
        arctan_argument_signum(&a_sign,&b_sign,phi);
      break;
      case LEFT:
        arctan_argument_signum(&c_sign,&a_sign,phi);
        if (cos_theta > 0) b_sign = 1;
        else		   b_sign = -1;
      break;
      case RIGHT:
        arctan_argument_signum(&c_sign,&b_sign,phi);
        if (cos_theta > 0) a_sign = 1;
        else		   a_sign = -1;
      break;
      case BACK:
        arctan_argument_signum(&b_sign,&c_sign,phi);
        if (cos_theta > 0) a_sign = 1;
        else		   a_sign = -1;
      break;
      case FRONT:
        arctan_argument_signum(&a_sign,&c_sign,phi);
        if (cos_theta > 0) b_sign = 1;
        else		   b_sign = -1;
      break;
      default:
        Error0(NO_OPTION);
    }
    
    X[0] = fabs(a)*a_sign;
    X[1] = fabs(b)*b_sign;
    
    /* check if x of X really gives you the correct angles */
    x_of_X(x,X,patch);
    x[0] -= c[0];
    x[1] -= c[1];
    x[2] -= c[2];
    r = root_square(3,x,0);
    theta2 = acos(x[2]/r);
    phi2   = arctan(x[1],x[0]);
    if (EQL(theta2,theta) && EQL(phi2,phi))
    {
      found_flg = YES;
      *ppatch = patch;
      break;
    }
  }
  if (found_flg == NO)
    Error0("(X,Y,Z) or patch could not be found.\n");
}





