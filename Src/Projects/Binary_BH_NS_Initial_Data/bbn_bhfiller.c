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

/* ->: EXIT_SUCCESS if succeeds, otherwise an error code
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

  int ret = 0;
  
  /* initialize */
  struct BHFiller_S *const bhf = bhf_init(grid,method);
  
  /* call bh-filler */
  ret = bhf->bhfiller(bhf);
  
  /* free */
  bhf_free(bhf);
  
  return ret;
}


/* initialize the bhfiller struct */
struct BHFiller_S* 
bhf_init
  (
  Grid_T *const grid/* the whole grid */,
  const char *const method/* the method to be used for extrapolating */
  )
{
  struct BHFiller_S *const bhf = calloc(1,sizeof(*bhf));IsNull(bhf);
  /* grid */
  bhf->grid = grid;
  /* method */
  sprintf(bhf->method,"%s",method);
  
  /* using Cheb_Tn and Ylm as the bases 
  // and extrapolate demanding C2 continuity */
  if (strcmp_i(method,"3rd_ChebTn_Ylm"))
  {
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
    unsigned f,nf,i,j,p;
    int lmax_par = PgetiEZ("bbn_bhfiller_lmax");
    if (lmax_par == INT_MAX)/* if par doesn't exist */
      lmax_par = 10;
    const unsigned lmax   = (unsigned)lmax_par;
    const unsigned Ntheta = 2*lmax+1;
    const unsigned Nphi   = 2*lmax+1;
    const unsigned N      = Ntheta*Nphi;
    
    nf = 0;/* number of fields */
    while(fields_name[nf]) ++nf;
    
    /* alloc */
    bhf->nf = nf;
    bhf->fld= calloc(nf,sizeof(*bhf->fld));IsNull(bhf->fld);
    
    /* set the method function */
    bhf->bhfiller = bhf_3rd_ChebTn_Ylm;
    
    /* collect names */
    collect_names(bhf,nf);
    
    /* initialize tables */
    init_Legendre_root_function();
    bhf->lmax   = lmax;
    bhf->Ntheta = Ntheta;
    bhf->Nphi   = Nphi;
    /* alloc radial_coeffs */
    for (f = 0; f < nf; ++f)
    {
      for (i = 0 ; i < 4; ++i)
      {
        bhf->fld[f]->radial_coeffs[i]  = alloc_double(N);
        bhf->fld[f]->realYlm_coeffs[i] = alloc_ClmYlm(lmax);
        bhf->fld[f]->imagYlm_coeffs[i] = alloc_ClmYlm(lmax);
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
    bhf->npo = npo;
    bhf->npi = npi;
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
  }/* if (strcmp_i(method,"3rd_ChebTn_Ylm")) */
  /* using George's thesis method */
  else if (strcmp_i(method,"WTGR"))
  {
    /* set the method function */
    bhf->bhfiller = bhf_WTGR;
  }
  else if (strcmp_i(method,"4th_Poly_Ylm"))
  {
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
    unsigned f,nf,i,j,p;
    int lmax_par = PgetiEZ("bbn_bhfiller_lmax");
    if (lmax_par == INT_MAX)/* if par doesn't exist */
      lmax_par = 10;
    const unsigned lmax   = (unsigned)lmax_par;
    const unsigned Ntheta = 2*lmax+1;
    const unsigned Nphi   = 2*lmax+1;
    const unsigned N      = Ntheta*Nphi;
    
    nf = 0;/* number of fields */
    while(fields_name[nf]) ++nf;
    
    /* alloc */
    bhf->nf = nf;
    bhf->fld= calloc(nf,sizeof(*bhf->fld));IsNull(bhf->fld);
    
    /* set the method function */
    bhf->bhfiller = bhf_4th_Poly_Ylm;
    
    /* collect names */
    collect_names(bhf,nf);
    
    /* initialize tables */
    init_Legendre_root_function();
    bhf->lmax   = lmax;
    bhf->Ntheta = Ntheta;
    bhf->Nphi   = Nphi;
    /* alloc radial_coeffs */
    for (f = 0; f < nf; ++f)
    {
      for (i = 0 ; i < 5; ++i)
      {
        bhf->fld[f]->radial_coeffs[i]  = alloc_double(N);
        bhf->fld[f]->realYlm_coeffs[i] = alloc_ClmYlm(lmax);
        bhf->fld[f]->imagYlm_coeffs[i] = alloc_ClmYlm(lmax);
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
    bhf->npo = npo;
    bhf->npi = npi;
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
  }
  /* solving elliptic equations in the hole */
  else if (strcmp_i(method,"C2_EllEq_Brown"))
  {
    const unsigned npi = 7;/* number of patches inside BH */
    const unsigned npo = 6;/* number of patches outside BH */
    unsigned nf,i,j,p;
    
    nf = 0;/* number of fields */
    while(fields_name[nf]) ++nf;
    
    /* alloc */
    bhf->nf = nf;
    bhf->fld= calloc(nf,sizeof(*bhf->fld));IsNull(bhf->fld);
    
    /* collect names */
    collect_names(bhf,nf);
        
    /* set the method function */
    bhf->bhfiller = bhf_ell_Brown;
    
    /* patches outside the BH */
    bhf->npo = npo;
    bhf->npi = npi;
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
  }
  else
    Error0(NO_OPTION);
  
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
    for (i = 0; i < MAX_COEFFS; ++i)
    {
      _free(bhf->fld[f]->radial_coeffs[i]);
      _free(bhf->fld[f]->realYlm_coeffs[i]);
      _free(bhf->fld[f]->imagYlm_coeffs[i]);
    }
  }
  free_2d_mem(bhf->fld,bhf->nf);
  _free(bhf->patches_outBH);
  _free(bhf->patches_inBH);
  free(bhf);
}

/* ->: EXIT_SUCESS if succeeds, otherwise an error code.
// method to fill BH is 3rd_ChebTn_Ylm:
// ===============================
//
// f(r,th,ph) = a0(th,ph)Tn(0,t)+a1(th,ph)Tn(1,t)+
//              a2(th,ph)Tn(2,t)+a3(th,ph)Tn(3,t)
// where, t = 2*r/rfill-1.
// the coeffs a's are determinded by demaning the C2 continuity
// across the AH and the value of the function at r = 0.
// */
static int bhf_3rd_ChebTn_Ylm(struct BHFiller_S *const bhf)
{
  printf("|--> BH-filler method = 3rd_ChebTn_Ylm.\n");
  fflush(stdout);
  Grid_T *const grid = bhf->grid;
  //const double EPS   = 1E-12;
  const unsigned NCoeffs = 4;
  const unsigned npo = bhf->npo;
  const unsigned npi = bhf->npi;
  const unsigned nf  = bhf->nf;/* numebr of fields */
  const unsigned lmax   = bhf->lmax;
  const unsigned Ntheta = bhf->Ntheta;
  const unsigned Nphi   = bhf->Nphi;
  const double rfill = Pgetd("r_excision");
  const double rfill3= pow(rfill,3);
  const double att_a = 0.8;
  const double att_b = -0.6;
  unsigned p,fld;
  
  /* update all coeffs to avoid race condition */
  printf("|--> Updating coefficients ...\n");
  fflush(stdout);
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = bhf->patches_outBH[p];
    unsigned f;

    bbn_1st_2nd_derivatives_conformal_metric(patch);
    bbn_add_and_take_2nd_derivatives_K(patch);
    Field_T *R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f;
    Field_T *R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
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
        make_coeffs_2d(patch->pool[f],0,1);/* X and Y direction */
    }
  }
  
  /* populating f, df/dr, d^2f/dr^2 at each (th,ph) points */
  printf("|--> Populating extrapolation coefficients ...\n");
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
        double ddfddr = 0,dfdr = 0,fr1 = 0,fr0 = 0;
        double a[NCoeffs];
        double _ddfddr[3] = {0,0,0};
        unsigned ij = IJ(i,j,Nphi);
        unsigned d1,d2;/* derivative */
        double X[3],x[3],_x[3],N[3];
        /* find patch for the given theta and phi */
        find_XYZ_and_patch_of_theta_phi_BH_CS(X,&patch,theta,phi,grid);

        /* r = rfill(sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^) */
        _x[0] = rfill*sin(theta)*cos(phi);
        _x[1] = rfill*sin(theta)*sin(phi);
        _x[2] = rfill*cos(theta);
        x[0]  = _x[0] + patch->c[0];
        x[1]  = _x[1] + patch->c[1];
        x[2]  = _x[2] + patch->c[2];
        assert(X_of_x(X,x,patch));
        
        /* normal vector */
        N[0]  = sin(theta)*cos(phi);
        N[1]  = sin(theta)*sin(phi);
        N[2]  = cos(theta);
        
        /* 2d interpolate on the surface */
        Interpolation_T *interp_s = init_interpolation();
        interp_s->XY_dir_flag = 1;
        interp_s->X = X[0];
        interp_s->Y = X[1];
        interp_s->K = 0;
        /* f value at r = r0 and r = 0*/
        interp_s->field = patch->pool[Ind(bhf->fld[fld]->f)];
        plan_interpolation(interp_s);
        fr1 = execute_interpolation(interp_s);
        fr0 = bhf->fld[fld]->f_r0;
        
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
            _ddfddr[_i] += (KD[_i==_j]/rfill - _x[_i]*_x[_j]/rfill3)*df_dx[_j];
            _ddfddr[_i] += N[_j]*ddf_ddx[IJsymm3(_i,_j)];
          }
          ddfddr += _ddfddr[_i]*N[_i]; \
        }
        
        a[0] = (10*fr0 + 22*fr1 + rfill*(-6*dfdr + ddfddr*rfill))/32.;
        a[1] = (-30*fr0 + 30*fr1 + rfill*(2*dfdr - ddfddr*rfill))/64.;
        a[2] = (6*fr0 - 6*fr1 + rfill*(6*dfdr - ddfddr*rfill))/32.;
        a[3] = (-2*fr0 + 2*fr1 + rfill*(-2*dfdr + ddfddr*rfill))/64.;
        
        for (_i = 0; _i < NCoeffs; _i++)
          bhf->fld[fld]->radial_coeffs[_i][ij] = a[_i];
        
      }
    }/* for (i = 0; i < Ntheta; ++i) */
    /* now populate the Ylm coeffs */
    for ( i = 0 ; i < NCoeffs; ++i)
    {
      double *rC       = bhf->fld[fld]->realYlm_coeffs[i];
      double *iC       = bhf->fld[fld]->imagYlm_coeffs[i];
      const double *v  = bhf->fld[fld]->radial_coeffs[i];
      get_Ylm_coeffs(rC,iC,v,Ntheta,Nphi,lmax);
    }
  }/* for (fld = 0; fld < nf ++fld) */
  
  /* now fill the BH */
  printf("|--> Fill the holes ...\n");
  fflush(stdout);
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npi; p++)
  {
    Patch_T *patch = bhf->patches_inBH[p];
    unsigned nn = patch->nn;
    double theta = 0,phi = 0, t = 0;
    double Y;
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
        
        if (r > rfill)
          r = rfill;
        
        theta = acos(z/r);
        phi = arctan(y,x);
        t   = 2*r/rfill-1;
        for (i = 0; i < NCoeffs; ++i)
        {
          const double *rC = bhf->fld[f]->realYlm_coeffs[i];
          const double *iC = bhf->fld[f]->imagYlm_coeffs[i];
          v[ijk] += 
            interpolation_Ylm(rC,iC,lmax,theta,phi)*Cheb_Tn((int)i,t);
        }
        Y = 0.5*(1+tanh(att_a*(rfill/(rfill-r)+att_b*(rfill/r))));
        int p = Pgetd("p");
        double Yn = pow(Y,p);
        v[ijk] = Yn*v[ijk]+(1-Yn)*bhf->fld[f]->f_r0;
        
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


/* for BAM initial data reader,
// extrapolate the fields Beta,eta, psi, _gamma's
// inside the BH, using the method developed by Wolfgang and Geroge,
// more info "http://fau.digital.flvc.org/islandora/object/fau%3A4224". */

/* handy macros for extrapolating inside BH: */
#define WTGR_StringIt(x)  #x

#define WTGR_EXTRAPOLATE_scalar(x)   \
        double x##_onAH       = interpolate_from_patch_prim(WTGR_StringIt(x)        ,X_on_BHsurf,BHsurf_patch); \
        double d##x##_D0_onAH = interpolate_from_patch_prim(WTGR_StringIt(d##x##_D0),X_on_BHsurf,BHsurf_patch); \
        double d##x##_D1_onAH = interpolate_from_patch_prim(WTGR_StringIt(d##x##_D1),X_on_BHsurf,BHsurf_patch); \
        double d##x##_D2_onAH = interpolate_from_patch_prim(WTGR_StringIt(d##x##_D2),X_on_BHsurf,BHsurf_patch); \
        double dur_##x        = (N[0]*d##x##_D0_onAH+N[1]*d##x##_D1_onAH+N[2]*d##x##_D2_onAH); \
        double ur_##x         = x##_onAH + dur_##x*dr; \
        x[ijk]                = ur_##x*Y + u0_##x*(1-Y);

#define WTGR_EXTRAPOLATE_Beta(x)   \
        double x##_onAH      = interpolate_from_patch_prim(WTGR_StringIt(x)       ,X_on_BHsurf,BHsurf_patch); \
        double d##x##D0_onAH = interpolate_from_patch_prim(WTGR_StringIt(d##x##D0),X_on_BHsurf,BHsurf_patch); \
        double d##x##D1_onAH = interpolate_from_patch_prim(WTGR_StringIt(d##x##D1),X_on_BHsurf,BHsurf_patch); \
        double d##x##D2_onAH = interpolate_from_patch_prim(WTGR_StringIt(d##x##D2),X_on_BHsurf,BHsurf_patch); \
        double dur_##x       = (N[0]*d##x##D0_onAH+N[1]*d##x##D1_onAH+N[2]*d##x##D2_onAH); \
        double ur_##x        = x##_onAH + dur_##x*dr; \
        x[ijk]               = ur_##x*Y + u0_##x*(1-Y);
        
#define WTGR_EXTRAPOLATE_gammabar(x)   \
        double x##_onAH      = interpolate_from_patch_prim(WTGR_StringIt(_##x)     ,X_on_BHsurf,BHsurf_patch); \
        double d##x##D0_onAH = interpolate_from_patch_prim(WTGR_StringIt(_d##x##D0),X_on_BHsurf,BHsurf_patch); \
        double d##x##D1_onAH = interpolate_from_patch_prim(WTGR_StringIt(_d##x##D1),X_on_BHsurf,BHsurf_patch); \
        double d##x##D2_onAH = interpolate_from_patch_prim(WTGR_StringIt(_d##x##D2),X_on_BHsurf,BHsurf_patch); \
        double dur_##x       = (N[0]*d##x##D0_onAH+N[1]*d##x##D1_onAH+N[2]*d##x##D2_onAH); \
        double ur_##x        = x##_onAH + dur_##x*dr; \
        _##x[ijk]            = ur_##x*Y + u0__##x*(1-Y);


static int bhf_WTGR(struct BHFiller_S *const bhf)
{
  printf("|--> BH-filler method = WTGR.\n");
  fflush(stdout);
  Grid_T *const grid          = bhf->grid;
  const double EPS            = 1E-12;/* to avoid division by zero */
  const double EPS2           = 1E-6;/* to increase r_fill radius a bit */
  const double r_fill         = Pgetd("BH_R_size")*(1+EPS2);
  const double Ma             = Pgetd("BH_irreducible_mass");
  const double att_a          = 0.8;
  const double att_b          = -0.6;
  const double u0_Beta_U0     = 0;
  const double u0_Beta_U1     = 0;
  const double u0_Beta_U2     = 0;
  const double u0__gamma_D0D0 = 1;
  const double u0__gamma_D0D1 = 0;
  const double u0__gamma_D0D2 = 0;
  const double u0__gamma_D1D1 = 1;
  const double u0__gamma_D1D2 = 0;
  const double u0__gamma_D2D2 = 1;
  const double u0_K           = 0;
  const double u0_alpha       = 0.1;
  Needle_T *patch_numbers = 0;
  const unsigned npi = 7;/* number of patches inside BH */
  const unsigned npo = 6;/* number of patches outside BH */
  unsigned p;
  
  /* check if it is perfect sphere */
  if (!Pcmps("BH_R_type","PerfectSphere"))
    bbn_bam_error("This function is used when "
        "the BH surface is a perfect sphere!",__FILE__,__LINE__);
  
  /* update coeffs to avoid race condition */
  patch_numbers       = alloc_needle();
  patch_numbers->grid = grid;
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_up",grid));
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_down",grid));
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_left",grid));
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_right",grid));
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_back",grid));
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_front",grid));
  if (patch_numbers->Nin != npo)
    bbn_bam_error("Wrong patch number",__FILE__,__LINE__);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = grid->patch[patch_numbers->in[p]];
    unsigned f;
    
    bbn_1st_derivatives_conformal_metric(patch);
    
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
  free_needle(patch_numbers);
  
  /* extrapolate */
  patch_numbers       = alloc_needle();
  patch_numbers->grid = grid;
  needle_in(patch_numbers,GetPatch("right_BH_up",grid));
  needle_in(patch_numbers,GetPatch("right_BH_down",grid));
  needle_in(patch_numbers,GetPatch("right_BH_left",grid));
  needle_in(patch_numbers,GetPatch("right_BH_right",grid));
  needle_in(patch_numbers,GetPatch("right_BH_back",grid));
  needle_in(patch_numbers,GetPatch("right_BH_front",grid));
  needle_in(patch_numbers,GetPatch("right_central_box",grid));
  if(patch_numbers->Nin != npi)
    bbn_bam_error("Wrong patch number",__FILE__,__LINE__);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npi; p++)
  {
    Patch_T *patch = grid->patch[patch_numbers->in[p]];
    unsigned nn = patch->nn;
    double Y;
    double theta = 0,phi = 0,dr = 0;
    double x_on_BHsurf[3]={0},
           X_on_BHsurf[3]={0},
           N[3] = {0};
    unsigned ijk;
    
    /* fill needle */
    Needle_T *needle = alloc_needle();
    needle->grid = grid;
    needle_in(needle,GetPatch("right_BH_surrounding_up",grid));
    needle_in(needle,GetPatch("right_BH_surrounding_down",grid));
    needle_in(needle,GetPatch("right_BH_surrounding_left",grid));
    needle_in(needle,GetPatch("right_BH_surrounding_right",grid));
    needle_in(needle,GetPatch("right_BH_surrounding_back",grid));
    needle_in(needle,GetPatch("right_BH_surrounding_front",grid));
    
    bbn_add_fields_in_patch(patch);
    REALLOC_v_WRITE_v(Beta_U0)
    REALLOC_v_WRITE_v(Beta_U1)
    REALLOC_v_WRITE_v(Beta_U2)
    REALLOC_v_WRITE_v(B0_U0)
    REALLOC_v_WRITE_v(B0_U1)
    REALLOC_v_WRITE_v(B0_U2)
    REALLOC_v_WRITE_v(psi)
    REALLOC_v_WRITE_v(eta)
    REALLOC_v_WRITE_v(K)
    REALLOC_v_WRITE_v(_gamma_D2D2)
    REALLOC_v_WRITE_v(_gamma_D0D2)
    REALLOC_v_WRITE_v(_gamma_D0D0)
    REALLOC_v_WRITE_v(_gamma_D0D1)
    REALLOC_v_WRITE_v(_gamma_D1D2)
    REALLOC_v_WRITE_v(_gamma_D1D1)
    
    REALLOC_v_WRITE_v(_gammaI_U2U2)
    REALLOC_v_WRITE_v(_gammaI_U0U2)
    REALLOC_v_WRITE_v(_gammaI_U0U0)
    REALLOC_v_WRITE_v(_gammaI_U0U1)
    REALLOC_v_WRITE_v(_gammaI_U1U2)
    REALLOC_v_WRITE_v(_gammaI_U1U1)
    
    /* making B1 */
    bbn_update_B1_U012(patch);
    READ_v(B1_U0)
    READ_v(B1_U1)
    READ_v(B1_U2)
    for (ijk = 0; ijk < nn; ++ijk)
    {
      DEF_RELATIVE_x
      DEF_RELATIVE_y
      DEF_RELATIVE_z
      DEF_RELATIVE_r
      Patch_T *BHsurf_patch = 0;
      if (!EQL(r,0))
      {
        N[0]  = x/r;
        N[1]  = y/r;
        N[2]  = z/r;
        theta = acos(z/r);
      }
      else
      {
        N[0]  = 0;
        N[1]  = 0;
        N[2]  = 0;
        r     = EPS;
        theta = 0;
      }
      dr    = r - r_fill;
      phi   = arctan(y,x);
      Y     = 0.5*(1+tanh(att_a*(r_fill/(r_fill-r)+att_b*(r_fill/r))));
      if(!isfinite(Y))
        bbn_bam_error("BH filler Y goes wrong.",__FILE__,__LINE__);
      
      x_on_BHsurf[0] = r_fill*sin(theta)*cos(phi)+patch->c[0];
      x_on_BHsurf[1] = r_fill*sin(theta)*sin(phi)+patch->c[1];
      x_on_BHsurf[2] = r_fill*cos(theta)         +patch->c[2];
      
      /* find the patch and X which has this point */
      needle->x = x_on_BHsurf;
      point_finder(needle);
      if (!needle->Nans)
        bbn_bam_error("Could not find the given point!\n",__FILE__,__LINE__);
      BHsurf_patch = grid->patch[needle->ans[0]];
      if(!X_of_x(X_on_BHsurf,x_on_BHsurf,BHsurf_patch))
        bbn_bam_error("X is wrong.",__FILE__,__LINE__);
      
      _free(needle->ans);
      needle->ans  = 0;
      needle->Nans = 0;
      
      /* extrapolate */
      double u0_psi = 2+Ma/(2*r);
      double u0_eta = u0_alpha*u0_psi;
      
      WTGR_EXTRAPOLATE_scalar(psi)
      WTGR_EXTRAPOLATE_scalar(eta)
      WTGR_EXTRAPOLATE_scalar(K)

      WTGR_EXTRAPOLATE_Beta(Beta_U0)
      WTGR_EXTRAPOLATE_Beta(Beta_U1)
      WTGR_EXTRAPOLATE_Beta(Beta_U2)
      
      B0_U0[ijk] = Beta_U0[ijk]-B1_U0[ijk];
      B0_U1[ijk] = Beta_U1[ijk]-B1_U1[ijk];
      B0_U2[ijk] = Beta_U2[ijk]-B1_U2[ijk];
      
      WTGR_EXTRAPOLATE_gammabar(gamma_D2D2)
      WTGR_EXTRAPOLATE_gammabar(gamma_D0D2)
      WTGR_EXTRAPOLATE_gammabar(gamma_D0D0)
      WTGR_EXTRAPOLATE_gammabar(gamma_D0D1)
      WTGR_EXTRAPOLATE_gammabar(gamma_D1D2)
      WTGR_EXTRAPOLATE_gammabar(gamma_D1D1)
      
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

        if(!EQL(delta_U1D1,1))  bbn_bam_error("_gammaI is not correct!\n",__FILE__,__LINE__);
        if(!EQL(delta_U0D1,0))  bbn_bam_error("_gammaI is not correct!\n",__FILE__,__LINE__);
        if(!EQL(delta_U0D2,0))  bbn_bam_error("_gammaI is not correct!\n",__FILE__,__LINE__);
        if(!EQL(delta_U1D2,0))  bbn_bam_error("_gammaI is not correct!\n",__FILE__,__LINE__);
        if(!EQL(delta_U0D0,1))  bbn_bam_error("_gammaI is not correct!\n",__FILE__,__LINE__);
        if(!EQL(delta_U2D1,0))  bbn_bam_error("_gammaI is not correct!\n",__FILE__,__LINE__);
        if(!EQL(delta_U2D2,1))  bbn_bam_error("_gammaI is not correct!\n",__FILE__,__LINE__);
        if(!EQL(delta_U2D0,0))  bbn_bam_error("_gammaI is not correct!\n",__FILE__,__LINE__);
        if(!EQL(delta_U1D0,0))  bbn_bam_error("_gammaI is not correct!\n",__FILE__,__LINE__);
      }
    }
    /* free */
    free_needle(needle);
    
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
    /* bbn_update_psi10A_UiUj(patch); */
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  /* free */
  free_needle(patch_numbers);
  bbn_rm_1st_derivatives_conformal_metric(GetPatch("right_BH_surrounding_up",grid));
  bbn_rm_1st_derivatives_conformal_metric(GetPatch("right_BH_surrounding_down",grid));
  bbn_rm_1st_derivatives_conformal_metric(GetPatch("right_BH_surrounding_left",grid));
  bbn_rm_1st_derivatives_conformal_metric(GetPatch("right_BH_surrounding_right",grid));
  bbn_rm_1st_derivatives_conformal_metric(GetPatch("right_BH_surrounding_back",grid));
  bbn_rm_1st_derivatives_conformal_metric(GetPatch("right_BH_surrounding_front",grid));
  /* free all coeffs */
  for (p = 0; p < grid->np; p++)
  {
    Patch_T *patch = grid->patch[p];
    unsigned f;
    for (f = 0; f < patch->nfld; ++f)
    {
      free_coeffs(patch->pool[f]);
    }
  }
  printf("|--> memory usege     = %0.2f(Gb)\n",how_much_memory("gb"));
  fflush(stdout);
  
  return EXIT_SUCCESS;
}

/* undef the macros */
#ifdef WTGR_StringIt
#undef WTGR_StringIt
#endif

#ifdef WTGR_EXTRAPOLATE_scalar
#undef WTGR_EXTRAPOLATE_scalar
#endif

#ifdef WTGR_EXTRAPOLATE_Beta
#undef WTGR_EXTRAPOLATE_Beta
#endif

#ifdef WTGR_EXTRAPOLATE_gammabar
#undef WTGR_EXTRAPOLATE_gammabar
#endif

/* given field name, X and patch, finds the value of the field in X  
// using interpolation.
// ->return value: f(X) */
static double interpolate_from_patch_prim(const char *const field,const double *const X,Patch_T *const patch)
{
  double interp;
  Interpolation_T *interp_s = init_interpolation();
  Field_T *const F_field    = patch->pool[Ind(field)];
  
  interp_s->field = F_field;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}

/* collect names of the fields and their derivatives */
static void collect_names(struct BHFiller_S *const bhf,const unsigned nf)
{
  const char *s = 0;
  unsigned f,i,j;
    
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
  }/* for (f = 0; f < nf; ++f) */
  
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
}

/* ->: EXIT_SUCESS if succeeds, otherwise an error code.
// method to fill BH is 4th_Poly_Ylm:
// ===============================
//
// f(r,th,ph) = a0(th,ph)+a1(th,ph)*r + a2(th,ph)*r^2 +
//               + a3(th,ph)*r^3 + a4(th,ph)*r^4.
// the coeffs a's are determinded by demaning the C2 continuity
// across the AH and the value of the function at r = 0 and a point
// near horizon from outside.
// */
static int bhf_4th_Poly_Ylm(struct BHFiller_S *const bhf)
{
  printf("|--> BH-filler method = 4th_Poly_Ylm.\n");
  fflush(stdout);
  Grid_T *const grid = bhf->grid;
  const double EPS   = 1E-12;
  const unsigned NCoeffs = 5;
  const unsigned npo = bhf->npo;
  const unsigned npi = bhf->npi;
  const unsigned nf  = bhf->nf;/* numebr of fields */
  const unsigned lmax   = bhf->lmax;
  const unsigned Ntheta = bhf->Ntheta;
  const unsigned Nphi   = bhf->Nphi;
  const double rfill = Pgetd("r_excision");
  const double rfill3= pow(rfill,3);
  const double r_2   = rfill*(1+0.01);/* a bit farther from AH */
  unsigned p,fld;
  
  /* update all coeffs to avoid race condition */
  printf("|--> Updating coefficients ...\n");
  fflush(stdout);
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = bhf->patches_outBH[p];
    unsigned f;

    bbn_1st_2nd_derivatives_conformal_metric(patch);
    bbn_add_and_take_2nd_derivatives_K(patch);
    Field_T *R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f;
    Field_T *R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
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
  printf("|--> Populating extrapolation coefficients ...\n");
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
        Patch_T *patch = 0;
        double KD[2]   = {0,1};
        unsigned ij    = IJ(i,j,Nphi);
        double df_dx[3]   = {0};
        double ddf_ddx[6] = {0};
        double ddfddr = 0,dfdr = 0,
               fr_0 = 0,fr_1 = 0,fr_2 = 0;
        double a[NCoeffs];
        double _ddfddr[3] = {0,0,0};
        unsigned d1,d2;/* derivative */
        double X[3],x[3],_x[3],N[3],Xr_2[3];
        /* find patch for the given theta and phi */
        find_XYZ_and_patch_of_theta_phi_BH_CS(X,&patch,theta,phi,grid);
        
        /* r = r_2 */
        _x[0] = r_2*sin(theta)*cos(phi);
        _x[1] = r_2*sin(theta)*sin(phi);
        _x[2] = r_2*cos(theta);
        x[0]  = _x[0] + patch->c[0];
        x[1]  = _x[1] + patch->c[1];
        x[2]  = _x[2] + patch->c[2];
        assert(X_of_x(Xr_2,x,patch));
        
        /* r = rfill(sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^) */
        _x[0] = rfill*sin(theta)*cos(phi);
        _x[1] = rfill*sin(theta)*sin(phi);
        _x[2] = rfill*cos(theta);
        x[0]  = _x[0] + patch->c[0];
        x[1]  = _x[1] + patch->c[1];
        x[2]  = _x[2] + patch->c[2];
        assert(X_of_x(X,x,patch));
        
        /* normal vector */
        N[0]  = sin(theta)*cos(phi);
        N[1]  = sin(theta)*sin(phi);
        N[2]  = cos(theta);
        
        /* 3d interpolate on the surface */
        Interpolation_T *interp_s = init_interpolation();
        interp_s->XYZ_dir_flag = 1;
        
        /* f value at r = r_2 */
        interp_s->X = Xr_2[0];
        interp_s->Y = Xr_2[1];
        interp_s->Z = Xr_2[2];
        interp_s->field = patch->pool[Ind(bhf->fld[fld]->f)];
        plan_interpolation(interp_s);
        fr_2 = execute_interpolation(interp_s);
        
        /* f value at r = r0 and r = 0 */
        interp_s->X = X[0];
        interp_s->Y = X[1];
        interp_s->Z = X[2];
        interp_s->field = patch->pool[Ind(bhf->fld[fld]->f)];
        plan_interpolation(interp_s);
        fr_1 = execute_interpolation(interp_s);
        /* r = 0 */
        fr_0 = bhf->fld[fld]->f_r0;
        
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
            _ddfddr[_i] += (KD[_i==_j]/rfill - _x[_i]*_x[_j]/rfill3)*df_dx[_j];
            _ddfddr[_i] += N[_j]*ddf_ddx[IJsymm3(_i,_j)];
          }
          ddfddr += _ddfddr[_i]*N[_i]; \
        }
        
        a[0] = fr_0;
        
        a[1] = (-2*fr_0*Power(r_2 - rfill,3)*(3*r_2 + rfill) + 
     2*fr_1*Power(r_2,2)*(3*Power(r_2,2) - 8*r_2*rfill + 6*Power(rfill,2)) + 
     rfill*(-2*dfdr*Power(r_2,2)*(2*Power(r_2,2) - 5*r_2*rfill + 3*Power(rfill,2)) + 
        rfill*(ddfddr*Power(r_2,2)*Power(r_2 - rfill,2) - 2*fr_2*Power(rfill,2))))/
   (2.*r_2*Power(r_2 - rfill,3)*rfill);
   
        a[2] = (6*fr_0*Power(r_2 - rfill,3)*(r_2 + rfill) - 
     6*fr_1*(Power(r_2,4) - 2*Power(r_2,3)*rfill + 2*r_2*Power(rfill,3)) + 
     rfill*(6*fr_2*Power(rfill,3) - 
        ddfddr*r_2*Power(r_2 - rfill,2)*rfill*(2*r_2 + rfill) + 
        6*dfdr*r_2*(Power(r_2,3) - 2*Power(r_2,2)*rfill + Power(rfill,3))))/
   (2.*r_2*Power(r_2 - rfill,3)*Power(rfill,2));
   
        a[3] = (-2*fr_0*Power(r_2 - rfill,3)*(r_2 + 3*rfill) + 
     2*fr_1*r_2*(Power(r_2,3) - 6*r_2*Power(rfill,2) + 8*Power(rfill,3)) + 
     rfill*(-6*fr_2*Power(rfill,3) + 
        ddfddr*r_2*Power(r_2 - rfill,2)*rfill*(r_2 + 2*rfill) - 
        2*dfdr*r_2*(Power(r_2,3) - 6*r_2*Power(rfill,2) + 5*Power(rfill,3))))/
   (2.*r_2*Power(r_2 - rfill,3)*Power(rfill,3));
   
        a[4] = (2*fr_0*Power(r_2 - rfill,3) - 2*fr_1*r_2*
      (Power(r_2,2) - 3*r_2*rfill + 3*Power(rfill,2)) + 
     rfill*(-(rfill*(ddfddr*r_2*Power(r_2 - rfill,2) - 2*fr_2*rfill)) + 
        2*dfdr*r_2*(Power(r_2,2) - 3*r_2*rfill + 2*Power(rfill,2))))/
   (2.*r_2*Power(r_2 - rfill,3)*Power(rfill,3));
        
        for (_i = 0; _i < NCoeffs; _i++)
          bhf->fld[fld]->radial_coeffs[_i][ij] = a[_i];
        
      }
    }/* for (i = 0; i < Ntheta; ++i) */
    /* now populate the Ylm coeffs */
    for ( i = 0 ; i < NCoeffs; ++i)
    {
      double *rC       = bhf->fld[fld]->realYlm_coeffs[i];
      double *iC       = bhf->fld[fld]->imagYlm_coeffs[i];
      const double *v  = bhf->fld[fld]->radial_coeffs[i];
      get_Ylm_coeffs(rC,iC,v,Ntheta,Nphi,lmax);
    }
  }/* for (fld = 0; fld < nf ++fld) */
  
  /* now fill the BH */
  printf("|--> Fill the holes ...\n");
  fflush(stdout);
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npi; p++)
  {
    Patch_T *patch = bhf->patches_inBH[p];
    unsigned nn = patch->nn;
    double theta = 0,phi = 0;
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
        phi = arctan(y,x);
        for (i = 0; i < NCoeffs; ++i)
        {
          const double *rC = bhf->fld[f]->realYlm_coeffs[i];
          const double *iC = bhf->fld[f]->imagYlm_coeffs[i];
          v[ijk] += 
            interpolation_Ylm(rC,iC,lmax,theta,phi)*pow(r,i);
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
    
  return EXIT_SUCCESS;
}

/* ->: EXIT_SUCCESS.
// filling the hole using elliptic equations (Brown's method)
// this is elliptic eq. of 3rd order i.e (Laplace^2)^3 f = 0.
// ref:https://arxiv.org/pdf/0809.3533.pdf 
// NOTE: this override solving parameters thus MUST be called
// at the very end when the ID constructor is done.
// NOTE: one must use finite center difference to calculate
// Laplace^2 f in order to achieve C^2 cos the values outside and inside
// of horizon get interwoven; but, in spectral method derivatives are
// computed in an individual patch thus one must achieve C^2 by imposing
// b.c on the horizon. BC for spectral method are not trivial though! */
static int bhf_ell_Brown(struct BHFiller_S *const bhf)
{
  printf("|--> BH-filler method = C2_EllEq_Brown.\n");
  fflush(stdout);
  const unsigned Nf  = bhf->nf;
  const unsigned Npi = bhf->npi;
  const unsigned Npo = bhf->npo;
  const double rfill = Pgetd("r_excision");
  const double rfill3= pow(rfill,3);
  char s_solve_eq[MAX_STR_LARGE] = {'\0'};
  Grid_T *inbh_grid  = calloc(1,sizeof(*inbh_grid));IsNull(inbh_grid);
  /* init equations */
  Solve_Equations_T *SolveEqs = init_solve_equations(bhf->grid);
  sEquation_T **field_eq = init_eq()/* field equation */,
            **bc_eq = init_eq()/* B.C. for the field */,
            **jacobian_field_eq = init_eq()/* jacobian for field equation */,
            **jacobian_bc_eq = init_eq()/* jacobian for B.C. */;
  unsigned f,p;
  
  /* updating coeffs */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < Npo; p++)
  {
    Patch_T *patch = bhf->patches_outBH[p];
    unsigned fn;

    bbn_1st_2nd_derivatives_conformal_metric(patch);
    bbn_add_and_take_2nd_derivatives_K(patch);
    Field_T *R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f;
    Field_T *R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
    if (R1_f)
      make_coeffs_2d(R1_f,0,1);/* X and Y direction */
    if (R2_f)
      make_coeffs_2d(R2_f,0,1);/* X and Y direction */

    /* make coeffs for all fields inside this patch */
    for (fn = 0; fn < patch->nfld; ++fn)
    {
      if (patch->pool[fn]->v      &&
          patch->pool[fn] != R1_f && 
          patch->pool[fn] != R2_f    )
        make_coeffs_2d(patch->pool[fn],0,1);/* X and Y direction */
    }
  }
  
  /* adding field */
  for (p = 0; p < Npi; ++p)
  {
    Patch_T *patch = bhf->patches_inBH[p];
    bbn_add_fields_in_patch(patch);
  }
  
  /* setup fields and ell. eq. and B.C. .
  // note: each field has 3 ell. eqs => 3 b.c. */
  for (f = 0; f < Nf; ++f)
  {
    char s[MAX_STR2] = {'\0'};
    
    sprintf(s,"%s.2,",bhf->fld[f]->f);
    strcat(s_solve_eq,s);
    sprintf(s,"%s.1,",bhf->fld[f]->f);
    strcat(s_solve_eq,s);
    sprintf(s,"%s.0,",bhf->fld[f]->f);
    strcat(s_solve_eq,s);
    
    /* adding fields */
    for (p = 0; p < Npi; ++p)
    {
      Patch_T *patch = bhf->patches_inBH[p];
    
      /* add 3 fields and their derivatives 
      // with suffix 0,1,2 higher order fields */
      /* field */
      sprintf(s,"%s.0",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"%s.1",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"%s.2",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      
      /* dfield/d? */
      sprintf(s,"d%s.0_D0",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"d%s.0_D1",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"d%s.0_D2",bhf->fld[f]->f);
      
      add_field(s,0,patch,YES);
      sprintf(s,"d%s.1_D0",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"d%s.1_D1",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"d%s.1_D2",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      
      sprintf(s,"d%s.2_D0",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"d%s.2_D1",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"d%s.2_D2",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      
      /* ddfield/d^2? */
      sprintf(s,"dd%s.0_D0D0",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"dd%s.0_D1D1",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"dd%s.0_D2D2",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      
      sprintf(s,"dd%s.1_D0D0",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"dd%s.1_D1D1",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"dd%s.1_D2D2",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      
      sprintf(s,"dd%s.2_D0D0",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"dd%s.2_D1D1",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"dd%s.2_D2D2",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      
      /* b.c field */
      sprintf(s,"bc_%s.0",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"bc_%s.1",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"bc_%s.2",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      
      /* source fields */
      sprintf(s,"src_%s.0",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"src_%s.1",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
      sprintf(s,"src_%s.2",bhf->fld[f]->f);
      add_field(s,0,patch,YES);
    }  
    /* adding eqs. (each field has 3 eqs.) */  
    sprintf(s,"eq_%s.0",bhf->fld[f]->f);
    add_eq(&field_eq,bbn_bhf_eq_Brown,s);
    sprintf(s,"eq_%s.1",bhf->fld[f]->f);
    add_eq(&field_eq,bbn_bhf_eq_Brown,s);
    sprintf(s,"eq_%s.2",bhf->fld[f]->f);
    add_eq(&field_eq,bbn_bhf_eq_Brown,s);
    
    sprintf(s,"bc_%s.0",bhf->fld[f]->f);
    add_eq(&bc_eq,bbn_bhf_bc_Brown,s);
    sprintf(s,"bc_%s.1",bhf->fld[f]->f);
    add_eq(&bc_eq,bbn_bhf_bc_Brown,s);
    sprintf(s,"bc_%s.2",bhf->fld[f]->f);
    add_eq(&bc_eq,bbn_bhf_bc_Brown,s);
    
    sprintf(s,"jacobian_eq_%s.0",bhf->fld[f]->f);
    add_eq(&jacobian_field_eq,bbn_bhf_jacobian_eq_Brown,s);
    sprintf(s,"jacobian_eq_%s.1",bhf->fld[f]->f);
    add_eq(&jacobian_field_eq,bbn_bhf_jacobian_eq_Brown,s);
    sprintf(s,"jacobian_eq_%s.2",bhf->fld[f]->f);
    add_eq(&jacobian_field_eq,bbn_bhf_jacobian_eq_Brown,s);
    
    sprintf(s,"jacobian_bc_%s.0",bhf->fld[f]->f);
    add_eq(&jacobian_bc_eq ,bbn_bhf_jacobian_bc_Brown,s);
    sprintf(s,"jacobian_bc_%s.1",bhf->fld[f]->f);
    add_eq(&jacobian_bc_eq ,bbn_bhf_jacobian_bc_Brown,s);
    sprintf(s,"jacobian_bc_%s.2",bhf->fld[f]->f);
    add_eq(&jacobian_bc_eq ,bbn_bhf_jacobian_bc_Brown,s);
  }/* for (f = 0; f < Nf; ++f) */
  
  /* populating bc values */
  OpenMP_Patch_Pragma(omp parallel for)
  for (f = 0; f < Nf; ++f)
  {
    char s[MAX_STR2] = {'\0'};
    unsigned pp;
    
    for (pp = 0; pp < Npi; ++pp)
    {
      Patch_T *patch = bhf->patches_inBH[pp];
      
      /* this patch does not have BC */
      if (strstr(patch->name,"right_central_box"))
        continue;
      
      const unsigned *n = patch->n;
      unsigned i,j,k,ijk,_i,_j;
      
      sprintf(s,"bc_%s.0",bhf->fld[f]->f);
      double *bc0 = patch->pool[Ind(s)]->v;
      
      sprintf(s,"bc_%s.1",bhf->fld[f]->f);
      double *bc1 = patch->pool[Ind(s)]->v;
      
      sprintf(s,"bc_%s.2",bhf->fld[f]->f);
      double *bc2 = patch->pool[Ind(s)]->v;
      
      for (i = 0; i < n[0]; ++i)
      {
        for (j = 0; j < n[1]; ++j)
        {
          ijk = L(n,i,j,n[2]-1);/* on the boundary */
          double x[3] = { patch->node[ijk]->x[0],
                          patch->node[ijk]->x[1],
                          patch->node[ijk]->x[2] };
          double Xp[3] = {0};/* (X,Y,Z)(x,y,z) in grid_prev */
          Patch_T *patchp = 0;/* patch in grid_prev contains (x,y,z) */
          double KD[2]      = {0,1};
          double df_dx[3]   = {0};
          double ddf_ddx[6] = {0};
          double ddfddr = 0,dfdr = 0,fr1;
          double _ddfddr[3] = {0,0,0};
          unsigned d1,d2;/* derivative */
          double _x[3],N[3],r;
          
          /* finding X and patch outside BH associated to x */
          assert(find_X_and_patch_outside_BH(x,"",bhf->grid,Xp,&patchp));
          
          _x[0] = x[0]-patch->c[0];
          _x[1] = x[1]-patch->c[1];
          _x[2] = x[2]-patch->c[2];
          r     = root_square(3,_x,0);
          
          /* normal vector */
          N[0]  = _x[0]/r;
          N[1]  = _x[1]/r;
          N[2]  = _x[2]/r;
          
          /* 2d interpolate on the surface */
          Interpolation_T *interp_s = init_interpolation();
          interp_s->XY_dir_flag = 1;
          interp_s->X = Xp[0];
          interp_s->Y = Xp[1];
          interp_s->K = 0;
          /* f value at r = r1 */
          interp_s->field = 
            patchp->pool[LookUpField_E(bhf->fld[f]->f,patchp)];
          plan_interpolation(interp_s);
          fr1 = execute_interpolation(interp_s);
          
          /* df/dx value */
          for (d1 = 0; d1 < 3; d1++)
          {
            interp_s->field = 
              patchp->pool[LookUpField_E(bhf->fld[f]->df[d1],patchp)];
            plan_interpolation(interp_s);
            df_dx[d1] = execute_interpolation(interp_s);
          }
          /* d^2f/dx^2 value */
          for (d1 = 0; d1 < 3; d1++)
          {
            for (d2 = d1; d2 < 3; d2++)
            {
              interp_s->field = 
                patchp->pool[LookUpField_E(bhf->fld[f]->ddf[IJsymm3(d1,d2)],patchp)];
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
              _ddfddr[_i] += (KD[_i==_j]/rfill - _x[_i]*_x[_j]/rfill3)*df_dx[_j];
              _ddfddr[_i] += N[_j]*ddf_ddx[IJsymm3(_i,_j)];
            }
            ddfddr += _ddfddr[_i]*N[_i]; \
          }
        
          /* now fill the bc values */
          for (k = 0; k < n[2]; ++k)
          {
            ijk = L(n,i,j,k);
            bc0[ijk] = fr1;
            bc1[ijk] = dfdr;
            bc2[ijk] = ddfddr;
          }
        
        }/* for (j = 0; j < n[1]; ++j) */
      }/* for (i = 0; i < n[0]; ++i) */
    }
  }/* for (f = 0; f < Nf; ++f) */
  
  /* initialize solving parameters */
  Psets("Solving_Order",s_solve_eq);
  Psetd("Solving_Residual",1E-10);
  Pseti("Solving_Max_Number_of_Newton_Step",3);
  Psetd("Solving_Newton_update_weight",1.0);
  initialize_solving_man(bhf->grid,field_eq,bc_eq,
                         jacobian_field_eq,jacobian_bc_eq);
  /* select inside bh patch to solve the field there */
  inbh_grid->kind = bhf->grid->kind;
  inbh_grid->gn   = bhf->grid->gn;
  inbh_grid->np   = Npi;
  inbh_grid->patch = calloc(Npi+1,sizeof(*inbh_grid->patch));
  IsNull(inbh_grid->patch);
  for (p = 0; p < Npi; ++p)
  {
    Patch_T *patch = bhf->patches_inBH[p];
    inbh_grid->patch[p] = calloc(1,sizeof(*inbh_grid->patch[p]));
    IsNull(inbh_grid->patch[p]);
    /* note, for the following all of the pointers inside the structures 
    // will be equal, since this is not a deep copy. */
    inbh_grid->patch[p][0]    = patch[0];
    inbh_grid->patch[p]->pn   = p;
    inbh_grid->patch[p]->grid = inbh_grid;
    inbh_grid->nn            += patch->nn;
    /* the following needs to be constructed from scratch */
    inbh_grid->patch[p]->interface   = 0;
    inbh_grid->patch[p]->solving_man = 
      calloc(1,sizeof(*inbh_grid->patch[p]->solving_man));
    IsNull(inbh_grid->patch[p]->solving_man);
    inbh_grid->patch[p]->solving_man[0]     = patch->solving_man[0];
    inbh_grid->patch[p]->solving_man->patch = inbh_grid->patch[p];
    inbh_grid->patch[p]->solving_man->jacobian = 0;
    inbh_grid->patch[p]->solving_man->nj       = 0;
    inbh_grid->patch[p]->solving_man->method->Schur_Complement = 0;
    inbh_grid->patch[p]->solving_man->method->SchurC = 0;
  }
  /* now let's fill up inbh_grid->patch[?]->interface */
  realize_geometry(inbh_grid);
  
  /* add special grids */
  for (f = 0; f < Nf; ++f)
  {
    char s[MAX_STR2] = {'\0'};
    /* adding special grid */
    sprintf(s,"%s.0",bhf->fld[f]->f);
    add_special_grid_solve_equations(inbh_grid,s,SolveEqs);
    sprintf(s,"%s.1",bhf->fld[f]->f);
    add_special_grid_solve_equations(inbh_grid,s,SolveEqs);
    sprintf(s,"%s.2",bhf->fld[f]->f);
    add_special_grid_solve_equations(inbh_grid,s,SolveEqs);
  }/* for (f = 0; f < Nf; ++f) */
    
  /* populating solution managment */
  SolveEqs->solving_order = Pgets("Solving_Order");
  SolveEqs->FieldUpdate   = bbn_bhf_ell_Brown_field_update;
  SolveEqs->SourceUpdate  = bbn_bhf_ell_Brown_source_update;
  SolveEqs->umfpack_refine= PgetdEZ("Solving_UMFPACK_refinement_step");
  SolveEqs->umfpack_size  = PgetiEZ("Solving_UMFPACK_size");
  
  /* solve equations */
  solve_eqs(SolveEqs);
  
  /* now fill the BH */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < Npi; p++)
  {
    Patch_T *patch = bhf->patches_inBH[p];
    unsigned nn = patch->nn;
    char s[MAX_STR2] = {'\0'};
    unsigned ijk;
    unsigned fn;
    
    /* copy the solutions */
    for (fn = 0; fn < Nf; ++fn)
    {
      sprintf(s,"%s.0",bhf->fld[fn]->f);
      const double *sol = patch->pool[Ind(s)]->v;
      Field_T *u = patch->pool[Ind(bhf->fld[fn]->f)];
      empty_field(u);
      u->v      = alloc_double(patch->nn);
      double *v = u->v;
      
      for (ijk = 0; ijk < nn; ++ijk)
        v[ijk] = sol[ijk];
    }/* for (f = 0; f < Nf ++f) */
    
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
  }/* for (p = 0; p < Npi; p++) */
  
  /* free */
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_up",bhf->grid));
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_down",bhf->grid));
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_left",bhf->grid));
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_right",bhf->grid));
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_back",bhf->grid));
  bbn_rm_1st_2nd_derivatives_conformal_metric
    (GetPatch("right_BH_surrounding_front",bhf->grid));
  
  /* free SolveEqs and special grid */
  free_solve_equations(SolveEqs);
  /* free data base of equations */
  free_db_eqs(field_eq);
  free_db_eqs(bc_eq);
  free_db_eqs(jacobian_field_eq);
  free_db_eqs(jacobian_bc_eq);     
  /* free grid and patches */
  if (inbh_grid)
  {
    FOR_ALL_PATCHES(p,inbh_grid)
    {
      Patch_T *patch = inbh_grid->patch[p];

      free_patch_interface(patch);
      free_patch_SolMan_method_Schur(patch);
      free(patch->solving_man);
    }
    free_2d_mem(inbh_grid->patch,inbh_grid->np);
    free(inbh_grid);
  }
  
  return EXIT_SUCCESS;                  
}

/* given a cartesian point x on the grid, it finds the corresponding X and patch 
// on which this x takes place. note: it only finds in patches covered outside of BH.
// hint, is the name of the patch that potentially has the given x 
// -> return value: for success 1, otherwise 0 */
static int find_X_and_patch_outside_BH(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch)
{
  Needle_T *needle = alloc_needle();
  const double LOW_RES_ERR = 1E-9;
  unsigned *found;
  unsigned p;
  int ret = 1;
  
  needle->grid = grid;
  needle->x    = x;
  
  /* find this point everywhere */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if (IsItHorizonPatch(patch))
      needle_in(needle,patch);
      
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  point_finder(needle);
  found = needle->ans;
  
  /* if it could not find X in neither the given hint patch nor its neighbors */
  if (needle->Nans)
  {
    *ppatch = grid->patch[found[0]];
    X_of_x(X,x,*ppatch);
  }
  else/* if no patch found let's find it in the other patches */
  {
    needle->ex   = needle->in;
    needle->Nex  = needle->Nin;
    needle->in   = 0;
    needle->Nin  = 0;
    
    point_finder(needle);
    found = needle->ans;
    
    /* if not found */
    if (!needle->Nans)
    {
      /* at the horizon of BH at low resulotion some points 
      // might not be found, let's fix this by hand. */
      if (strstr(hint,"right_BH"))
      {
        *ppatch = GetPatch(hint,grid);
        X_of_x(X,x,*ppatch);
        
        /* NOTE: This is assumed Cubed Spherical coords */
        if (LSS(fabs(X[2]),LOW_RES_ERR))
        {
          X[2] = 0;/* the culprit is at low resolution, 
                   // X[2] won't be found very close to 0! */
          /* make sure the X falls in the interval */
          assert(LSSEQL(X[0],1) && GRTEQL(X[0],-1));
          assert(LSSEQL(X[1],1) && GRTEQL(X[1],-1));
        }
        else
        {
          fprintf(stderr,"The point (%g,%g,%g) could not be found!\n",x[0],x[1],x[2]);
          fflush(stdout);
          fflush(stderr);
          ret = 0;
          /* Error0("Point not found!\n"); */
        }
      }
      else
      {
        fprintf(stderr,"The point (%g,%g,%g) could not be found!\n",x[0],x[1],x[2]);
        fflush(stdout);
        fflush(stderr);
        ret = 0;
        /* Error0("Point not found!\n"); */
      }
    }
    else
    {  
      *ppatch = grid->patch[found[0]];
      X_of_x(X,x,*ppatch);
    }
  }
  free_needle(needle);
  
  return ret;
}


/* field update: updating the derivatives, this is used in 
// Brown's approach ell eq. */
static void 
bbn_bhf_ell_Brown_field_update
  (
  Patch_T *const patch,
  const char *const name
  )
{
  char s[MAX_STR2] = {'\0'};
  
  /* field */
  Field_T *const f = patch->pool[Ind(name)];
  
  /* dfield/d? */
  sprintf(s,"d%s_D0",name);
  Field_T *const df_dx = patch->pool[Ind(s)];
  empty_field(df_dx);
  df_dx->v = Partial_Derivative(f,"x");
  
  sprintf(s,"d%s_D1",name);
  Field_T *const df_dy = patch->pool[Ind(s)];
  empty_field(df_dy);
  df_dy->v = Partial_Derivative(f,"y"); 
  
  sprintf(s,"d%s_D2",name);
  Field_T *const df_dz = patch->pool[Ind(s)];
  empty_field(df_dz);
  df_dz->v = Partial_Derivative(f,"z"); 
  
  /* d^2field/d?^2 */
  sprintf(s,"dd%s_D0D0",name);
  Field_T *const ddf_ddx = patch->pool[Ind(s)];
  empty_field(ddf_ddx);
  ddf_ddx->v = Partial_Derivative(df_dx,"x");
  
  sprintf(s,"dd%s_D1D1",name);
  Field_T *const ddf_ddy = patch->pool[Ind(s)];
  empty_field(ddf_ddy);
  ddf_ddy->v = Partial_Derivative(df_dy,"y");
  
  sprintf(s,"dd%s_D2D2",name);
  Field_T *const ddf_ddz = patch->pool[Ind(s)];
  empty_field(ddf_ddz);
  ddf_ddz->v = Partial_Derivative(df_dz,"z");
}

/* source update: updating sources of ell. eqs., this is used in 
// Brown's approach ell eqs. 
// note: we have:
// Laplac^2 psi_{n}- psi_{n+1} = 0, here psi_{n+1} is the source. */
static void 
bbn_bhf_ell_Brown_source_update
  (
  Grid_T *const grid,
  const char *const name
  )
{
  char s1[MAX_STR2] = {'\0'};
  char s2[MAX_STR2] = {'\0'};
  char s[MAX_STR] = {'\0'};
  const char *s_last;
  int f_index;
  unsigned p;
  
  sprintf(s,"%s",name);
  s_last  = name+strlen(name)-1;/* => 0 in psi.0 */
  f_index = atoi(s_last);/* => 0 in psi.0 */
  
  if (f_index > 0) /* -> psi.1 => psi.0 */
  {
    sprintf(&s[strlen(name)-1],"%d",f_index-1);
    sprintf(s1,"src_%s",s);
    sprintf(s2,"%s",name);
  }
    
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Field_T *src1,*src2;
    unsigned nn = patch->nn;
    unsigned ijk;
    
    if (f_index == 0)
      break;
      
    /* src1 */
    src1 = patch->pool[Ind(s1)];
    empty_field(src1);
    src1->v = alloc_double(patch->nn);/* set to zero */
    /* src2 */
    src2 = patch->pool[Ind(s2)];
    
    /* copy */
    for (ijk = 0; ijk < nn; ++ijk)
      src1->v[ijk] = src2->v[ijk];
  }
  
}



