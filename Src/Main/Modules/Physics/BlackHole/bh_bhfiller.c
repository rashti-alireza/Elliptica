/*
// Alireza Rashti
// September 2020
*/
/* filling black hole inside. */

#include "bh_bhfiller.h"


/* extrapolating inside the BH */
int bh_fill_inside_black_hole(Physics_T *const phys)
{
  FUNC_TIC
  
  int ret = -1;
  
  IF_sval("BHfiller_method","ChebTn_Ylm")
  {
    /* these fields to be extrapolated  */
    const char *fields_name[] = 
    {
      "psi","eta","K",
      "Beta_U0","Beta_U1","Beta_U2",
      "K_DiDj_D0D0","K_DiDj_D0D1",
      "K_DiDj_D0D2","K_DiDj_D1D1",
      "K_DiDj_D1D2","K_DiDj_D2D2",
      0/* end */
    };
    ret = bh_bhfiller(phys,fields_name,"ChebTn_Ylm");
  }
  else
  {
    Error0(NO_OPTION);
  }

  FUNC_TOC;
  
  return ret;
}

int 
bh_bhfiller
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  )
{
  if (phys->ctype != BH)
   Error0("Wrong physics!");

  printf(Pretty0"method = %s\n",method);
  
  int ret = -1;
  
  /* initialize */
  struct BHFiller_S *const bhf = bhf_init(phys,fields_name,method);
  
  /* call bh-filler */
  ret = bhf->bhfiller(bhf);
  
  /* free */
  bhf_free(bhf);
  
  return ret;
}


/* initialize the bhfiller struct */
static struct BHFiller_S* 
bhf_init
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  )
{
  struct BHFiller_S *const bhf = calloc(1,sizeof(*bhf));IsNull(bhf);
  /* grid */
  Grid_T *const grid = phys->grid;
  bhf->grid = grid;
  /* method */
  sprintf(bhf->method,"%s",method);
  
  /* using Cheb_Tn and Ylm as the bases 
  // and extrapolate demanding C2 continuity.
  // NOTE: among the other methods this is the best. */
  if (strcmp_i(method,"ChebTn_Ylm"))
  {
    const unsigned NCoeffs = 10;/* number of coeffs in ChebTn expansion */
    const unsigned lmax   = (unsigned)Geti("bhfiller_Ylm_expansion_lmax");
    const unsigned Ntheta = Ntheta_Ylm(lmax);
    const unsigned Nphi   = Nphi_Ylm(lmax);
    const unsigned N      = Ntotal_Ylm(lmax);
    unsigned npi;/* number of patches inside BH */
    unsigned npo;/* number of patches outside BH */
    /* values of extrapolant function at the center of BH f(r=0) */
    const double fr0_Beta_U0     = 0;
    const double fr0_Beta_U1     = 0;
    const double fr0_Beta_U2     = 0;
    const double fr0_gConf_D0D0  = 1;
    const double fr0_gConf_D0D1  = 0;
    const double fr0_gConf_D0D2  = 0;
    const double fr0_gConf_D1D1  = 1;
    const double fr0_gConf_D1D2  = 0;
    const double fr0_gConf_D2D2  = 1;
    const double fr0_K_DiDj_D0D0 = 0;
    const double fr0_K_DiDj_D0D1 = 0;
    const double fr0_K_DiDj_D0D2 = 0;
    const double fr0_K_DiDj_D1D1 = 0;
    const double fr0_K_DiDj_D1D2 = 0;
    const double fr0_K_DiDj_D2D2 = 0;
    const double fr0_K           = 0;
    const double fr0_alpha       = 0.2;
    const double fr0_psi         = 2;/* big enough */
    const double fr0_eta         = fr0_alpha*fr0_psi;
    unsigned f,nf,i;
    
    nf = 0;/* number of fields */
    while(fields_name[nf]) ++nf;
    
    /* alloc */
    bhf->nf = nf;
    bhf->fld= calloc(nf,sizeof(*bhf->fld));IsNull(bhf->fld);
    
    if (grid->kind == Grid_SplitCubedSpherical_NSNS ||
        grid->kind == Grid_SplitCubedSpherical_BHBH ||
        grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH  ||
        grid->kind == Grid_SplitCubedSpherical_SNS  ||
        grid->kind == Grid_CubedSpherical_BHNS      ||
        grid->kind == Grid_CubedSpherical_NSNS
       )
    {
     /* set the method function */
     IF_sval("surface_type","perfect_s2")
      bhf->bhfiller = bhf_ChebTn_Ylm_pefect_S2_CS;
     else
      Error0(NO_OPTION);
    }
    else
     Error0(NO_OPTION);
    
    /* collect names */
    collect_names(bhf,fields_name,nf);
    
    /* initialize tables */
    init_Legendre_root_function();
    bhf->lmax   = lmax;
    bhf->Ntheta = Ntheta;
    bhf->Nphi   = Nphi;
    bhf->NCoeffs= NCoeffs;
    
    /* alloc radial_coeffs */
    for (f = 0; f < nf; ++f)
    {
      for (i = 0 ; i < NCoeffs; ++i)
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
        bhf->fld[f]->f_r0    = fr0_psi;
        bhf->fld[f]->func_r0 = punc_psi;
      }
      else if (strcmp_i(fields_name[f],"eta"))
      {
        bhf->fld[f]->f_r0    = fr0_eta;
        bhf->fld[f]->func_r0 = punc_eta;
      }
      else if (strcmp_i(fields_name[f],"K"))
      {
        bhf->fld[f]->f_r0    = fr0_K;
        bhf->fld[f]->func_r0 = punc_K;
      }
      else if (strcmp_i(fields_name[f],"Beta_U0"))
      {
        bhf->fld[f]->f_r0    = fr0_Beta_U0;
        bhf->fld[f]->func_r0 = punc_Beta_U0;
      }
      else if (strcmp_i(fields_name[f],"Beta_U1"))
      {
        bhf->fld[f]->f_r0    = fr0_Beta_U1;
        bhf->fld[f]->func_r0 = punc_Beta_U1;
      }
      else if (strcmp_i(fields_name[f],"Beta_U2"))
      {
        bhf->fld[f]->f_r0    = fr0_Beta_U2;
        bhf->fld[f]->func_r0 = punc_Beta_U2;
      }
      else if (strcmp_i(fields_name[f],"gConf_D2D2"))
      {
        bhf->fld[f]->f_r0    = fr0_gConf_D2D2;
        bhf->fld[f]->func_r0 = punc_gConf_D2D2;
      }
      else if (strcmp_i(fields_name[f],"gConf_D0D2"))
      {
        bhf->fld[f]->f_r0    = fr0_gConf_D0D2;
        bhf->fld[f]->func_r0 = punc_gConf_D0D2;
      }
      else if (strcmp_i(fields_name[f],"gConf_D0D0"))
      {
        bhf->fld[f]->f_r0    = fr0_gConf_D0D0;
        bhf->fld[f]->func_r0 = punc_gConf_D0D0;
      }
      else if (strcmp_i(fields_name[f],"gConf_D0D1"))
      {
        bhf->fld[f]->f_r0    = fr0_gConf_D0D1;
        bhf->fld[f]->func_r0 = punc_gConf_D0D1;
      }
      else if (strcmp_i(fields_name[f],"gConf_D1D2"))
      {
        bhf->fld[f]->f_r0    = fr0_gConf_D1D2;
        bhf->fld[f]->func_r0 = punc_gConf_D1D2;
      }
      else if (strcmp_i(fields_name[f],"gConf_D1D1"))
      {
        bhf->fld[f]->f_r0    = fr0_gConf_D1D1;
        bhf->fld[f]->func_r0 = punc_gConf_D1D1;
      }
      else if (strcmp_i(fields_name[f],"K_DiDj_D2D2"))
      {
        bhf->fld[f]->f_r0    = fr0_K_DiDj_D2D2;
        bhf->fld[f]->func_r0 = 0;
      }
      else if (strcmp_i(fields_name[f],"K_DiDj_D0D2"))
      {
        bhf->fld[f]->f_r0    = fr0_K_DiDj_D0D2;
        bhf->fld[f]->func_r0 = 0;
      }
      else if (strcmp_i(fields_name[f],"K_DiDj_D0D0"))
      {
        bhf->fld[f]->f_r0    = fr0_K_DiDj_D0D0;
        bhf->fld[f]->func_r0 = 0;
      }
      else if (strcmp_i(fields_name[f],"K_DiDj_D0D1"))
      {
        bhf->fld[f]->f_r0    = fr0_K_DiDj_D0D1;
        bhf->fld[f]->func_r0 = 0;
      }
      else if (strcmp_i(fields_name[f],"K_DiDj_D1D2"))
      {
        bhf->fld[f]->f_r0    = fr0_K_DiDj_D1D2;
        bhf->fld[f]->func_r0 = 0;
      }
      else if (strcmp_i(fields_name[f],"K_DiDj_D1D1"))
      {
        bhf->fld[f]->f_r0    = fr0_K_DiDj_D1D1;
        bhf->fld[f]->func_r0 = 0;
      }
      else
        Error0(NO_OPTION);
      
    }/* for (f = 0; f < nf ++f) */
    
    /* patches outside the BH */
    bhf->patches_outBH = 
     collect_patches(phys->grid,Ftype("BH_around_IB"),&npo);
    IsNull(bhf->patches_outBH);
    bhf->npo = npo;
    /* patches inside the BH */
    bhf->patches_inBH = 
     collect_patches(phys->grid,Ftype("BH"),&npi);
    IsNull(bhf->patches_inBH);
    bhf->npi = npi;
    
  }/* if (strcmp_i(method,"ChebTn_Ylm")) */
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
// method to fill BH is ChebTn_Ylm with the following extrapolant:
// ===============================================================
//
// f(r(t),th,ph) = C_{ilm}*ChebT_i(t)*Y_{lm}(th,ph)
//              
// where, t = 2*r/rfill-1.
// the coeffs a's are determinded by demaning the C2 continuity
// across the AH and the value of the function at r = 0.
// note: it has some assumptions which are only true in cubed spherical
// and only for perfect sphere. */
static int bhf_ChebTn_Ylm_pefect_S2_CS(struct BHFiller_S *const bhf)
{
  const unsigned NCoeffs = bhf->NCoeffs;
  const unsigned npo     = bhf->npo;
  const unsigned npi     = bhf->npi;
  const unsigned nf     = bhf->nf;/* numebr of fields */
  const unsigned lmax   = bhf->lmax;
  const unsigned Ntheta = bhf->Ntheta;
  const unsigned Nphi   = bhf->Nphi;
  const double rfill = Pgetd("perfect_S2_radius");
  const double rfill3= pow(rfill,3);
  unsigned p,fld;

  /* update all coeffs to avoid race condition */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = bhf->patches_outBH[p];
    unsigned f = 0;

    /* make coeffs in  X and Y direction inside this patch */
    for (f = 0; f < nf; ++f)
    {
      int ii;
      
      /* must have the field */
      make_coeffs_2d(patch->fields[Ind(bhf->fld[f]->f)],0,1);
      
      /* must have dfield */
      for (ii = 0; ii < 3; ++ii)
      {
        int indxf = _Ind(bhf->fld[f]->df[ii]);
        if (indxf >= 0)
        {
         make_coeffs_2d(patch->fields[indxf],0,1);
        }
        else
        {
        printf(Pretty0"could not find %s in %s "
                       "=> calculate it numerically ...\n",
                       bhf->fld[f]->df[ii],patch->name);
         bhf->fld[f]->did_add_df = 1;
         Field_T *df = add_field(bhf->fld[f]->df[ii],0,patch,NO);
         partial_derivative(df);
        }
      }
      
      /* must have ddfield */
      for (ii = 0; ii < 6; ++ii)
      {
        int indxf = _Ind(bhf->fld[f]->ddf[ii]);
        if (indxf >= 0)
        {
          make_coeffs_2d(patch->fields[indxf],0,1);
        }
        else
        {
         printf(Pretty0"could not find %s in %s "
                       "=> calculate it numerically ...\n",
                       bhf->fld[f]->ddf[ii],patch->name);
         bhf->fld[f]->did_add_ddf = 1;
         Field_T *ddf = add_field(bhf->fld[f]->ddf[ii],0,patch,NO);
         partial_derivative(ddf);
        }
      }
      
      ++f;
    }
  }

  /* populating f, df/dr, d^2f/dr^2 at each (th,ph) points */
  printf(Pretty0"Populating extrapolation coefficients ...\n");
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
        unsigned ij = IJ_Ylm(i,j,Nphi);
        unsigned d1,d2;/* derivative */
        double X[3],x[3],_x[3],N[3];
        
        /* find patch for the given theta and phi */
        X[2] = 0.;
        find_XYZ_and_patch_of_theta_phi_CS
         (X,&patch,theta,phi,bhf->patches_outBH,bhf->npo);
         
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
        interp_s->field = patch->fields[Ind(bhf->fld[fld]->f)];
        plan_interpolation(interp_s);
        fr1 = execute_interpolation(interp_s);
        fr0 = bhf->fld[fld]->f_r0;
        
        /* df/dx value */
        for (d1 = 0; d1 < 3; d1++)
        {
          interp_s->field = patch->fields[Ind(bhf->fld[fld]->df[d1])];
          plan_interpolation(interp_s);
          df_dx[d1] = execute_interpolation(interp_s);
        }
        /* d^2f/dx^2 value */
        for (d1 = 0; d1 < 3; d1++)
        {
          for (d2 = d1; d2 < 3; d2++)
          {
            interp_s->field = 
              patch->fields[Ind(bhf->fld[fld]->ddf[IJsymm3(d1,d2)])];
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
        
        /* find and set the coeffs */
        bh_bhf_ChebTn_extrapolate
        (a,fr0,fr1,dfdr,ddfddr,rfill,NCoeffs);
        
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
  printf(Pretty0"Fill the hole ...\n");
  fflush(stdout);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npi; p++)
  {
    Patch_T *patch = bhf->patches_inBH[p];
    unsigned nn    = patch->nn;
    unsigned f,ijk;
    
    /* loop over all fields to be extrapolated */
    for (f = 0; f < nf; ++f)
    {
      Field_T *u = patch->fields[Ind(bhf->fld[f]->f)];
      empty_field(u);
      u->v      = alloc_double(patch->nn);
      double *v = u->v;
      double theta,phi,t;
      unsigned i;
      
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
        
      }/* for (ijk = 0; ijk < nn; ++ijk) */
    }/* for (f = 0; f < nf ++f) */
  }
    
  return EXIT_SUCCESS;
}

/* collect names of the fields and their derivatives */
static void collect_names(struct BHFiller_S *const bhf,const char **const fields_name,const unsigned nf)
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

/* puncture behavior function */
static double punc_psi(void *const params)
{
  struct Param_S *const par = params;
  double r     = par->r;
  //double rfill = par->rfill;
  double M     = par->M;
  double eps   = par->eps;
  
  if (r < eps)
    r = eps;
  
  return 2+M/(2*r);
}

/* puncture behavior function */
static double punc_eta(void *const params)
{
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //double eps   = par->eps;
  
  return 0.1*punc_psi(params);
}

/* puncture behavior function */
static double punc_K(void *const params)
{
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //double eps   = par->eps;
  UNUSED(params);
  return 0;
}

/* puncture behavior function */
static double punc_Beta_U0(void *const params)
{
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //double eps   = par->eps;
  UNUSED(params);
  return 0;
}

/* puncture behavior function */
static double punc_Beta_U1(void *const params)
{
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //const double eps   = par->eps;
  UNUSED(params);
  return 0;
}

/* puncture behavior function */
static double punc_Beta_U2(void *const params)
{
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //double eps   = par->eps;
  UNUSED(params);
  return 0;
}

/* puncture behavior function */
static double punc_gConf_D0D0(void *const params)
{
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //double eps   = par->eps;
  UNUSED(params);
  return 1;
}

/* puncture behavior function */
static double punc_gConf_D0D1(void *const params)
{
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //double eps   = par->eps;
  UNUSED(params);
  return 0;
}

/* puncture behavior function */
static double punc_gConf_D0D2(void *const params)
{
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //double eps   = par->eps;
  UNUSED(params);
  return 0;
}

/* puncture behavior function */
static double punc_gConf_D1D1(void *const params)
{
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //double eps   = par->eps;
  UNUSED(params);
  return 1;
}

/* puncture behavior function */
static double punc_gConf_D1D2(void *const params)
{	
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //double eps   = par->eps;
  UNUSED(params);
  return 0;
}

/* puncture behavior function */
static double punc_gConf_D2D2(void *const params)
{
  //struct Param_S *const par = params;
  //double r     = par->r;
  //double rfill = par->rfill;
  //double M     = par->M;
  //double eps   = par->eps;
  UNUSED(params);
  return 1;
}

/* smoothing function for inside of the BH */
double bh_bhf_smoother(const double r, const double rmax,const double rmin)
{
  double ret = 0;
  
  if (r >= rmax)
    return 1;
  else if (r <= rmin)
    return 0;
  else
  {
    return bh_bhf_poly_smoother(r,rmax,rmin);
    
    if (0)
      return polynomial5(r,rmax,rmin);
    if (0)
      return polynomial7(r,rmax,rmin);
  }
  
  return ret;
}

/* polynomial of order 5 with coeffs a, this polynomial is:
// P(rmin) = P'(rmin) = P"(rmin) = 0 and
// P(rmax) = 1 and P'(rmax) = P"(rmax) = 0. */
static double polynomial5(const double r, const double rmax,const double rmin)
{
  double a[6] = {0};
  double ret;
  
  a[0] = (Power(rmin,3)*(10*Power(rmax,2) - 5*rmax*rmin + Power(rmin,2)))/
   Power(-rmax + rmin,5);
   
  a[1] = (30*Power(rmax,2)*Power(rmin,2))/Power(rmax - rmin,5);
  
  a[2] = (-30*rmax*rmin*(rmax + rmin))/Power(rmax - rmin,5);
  
  a[3] = (10*(Power(rmax,2) + 4*rmax*rmin + Power(rmin,2)))/
   Power(rmax - rmin,5);
   
  a[4] = (-15*(rmax + rmin))/Power(rmax - rmin,5);
  
  a[5] = 6./Power(rmax - rmin,5);
 
  ret  = a[0] + a[1]*r   + a[2]*Power(r,2) + 
         a[3]*Power(r,3) + a[4]*Power(r,4) + 
         a[5]*Power(r,5);
          
  return ret;
}

/* polynomial of order 7 with coeffs a, this polynomial is:
// P(rmin) = P'(rmin) = P''(rmin) = P'''(rmin) = 0 and
// P(rmax) = 1 and P'(rmax) = P''(rmax) = P'''(rmax) = 0. */
static double polynomial7(const double r, const double rmax,const double rmin)
{
  double a[8] = {0};
  double ret;
  
  a[0] = (Power(rmin,4)*(35*Power(rmax,3) - 21*Power(rmax,2)*rmin + 
       7*rmax*Power(rmin,2) - Power(rmin,3)))/Power(rmax - rmin,7);
  
  a[1] = (-140*Power(rmax,3)*Power(rmin,3))/Power(rmax - rmin,7);
 
  a[2] = (210*Power(rmax,2)*Power(rmin,2)*(rmax + rmin))/
   Power(rmax - rmin,7);
  
  a[3] = (-140*rmax*rmin*(Power(rmax,2) + 3*rmax*rmin + Power(rmin,2)))/
   Power(rmax - rmin,7);
  
  a[4] = (35*(Power(rmax,3) + 9*Power(rmax,2)*rmin + 9*rmax*Power(rmin,2) + 
       Power(rmin,3)))/Power(rmax - rmin,7);
  
  a[5] = (-84*(Power(rmax,2) + 3*rmax*rmin + Power(rmin,2)))/
   Power(rmax - rmin,7);
  
  a[6] = (70*(rmax + rmin))/Power(rmax - rmin,7);
  
  a[7] = -20/Power(rmax - rmin,7);
  
  ret  = a[0] + a[1]*r   + a[2]*Power(r,2) + 
         a[3]*Power(r,3) + a[4]*Power(r,4) + 
         a[5]*Power(r,5) + a[6]*Power(r,6) +
         a[7]*Power(r,7);
          
  return ret;
}


