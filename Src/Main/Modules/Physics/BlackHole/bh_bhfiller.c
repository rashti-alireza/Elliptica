/*
// Alireza Rashti
// September 2020
*/
/* filling black hole inside. */

#include "bh_bhfiller.h"

/* make extrapolation weaker or stronger */
static const double Attenuate_Factor0 = 1;/* came from experiment */

/* extrapolating inside the BH */
int bh_fill_inside_black_hole(Physics_T *const phys)
{
  FUNC_TIC
  
  IF_sval("filler_method","none")
  {
    printf(Pretty0"No black hole filling requested.\n");
    FUNC_TOC
    return EXIT_SUCCESS;
  }
  
  Verbose = Geti("filler_verbose");
  
  int ret = EXIT_SUCCESS;
  /* these fields to be extrapolated  */
  char **fields_name = 
         read_separated_items_in_string(Gets("filler_fields"),',');
  Grid_T *grid;
  Uint f,p;
   
  printf(Pretty0"fields = %s\n",Gets("filler_fields"));
   
  /* first add patches */
  bh_add_patch_inside_black_hole(phys,Ftype("BH"));
  
  /* add fields if they not exists */
  grid = mygrid(phys,Ftype("BH"));
  f = 0;
  while (fields_name[f])
  {
    FOR_ALL_PATCHES(p,grid)
    {
     Patch_T *patch = grid->patch[p];
     
     /* if not there */
     if (_Ind(fields_name[f]) < 0)
      add_field(fields_name[f],0,patch,NO);
    }
    ++f;
  }
  
  /* now fill */
  ret = bh_bhfiller(phys,fields_name,Gets("filler_method"));

  free_2d(fields_name);
  
  FUNC_TOC;
  return ret;
}

int 
bh_bhfiller
  (
  Physics_T *const phys/* physics of interest */,
  char **const fields_name/* ends determined by 0 */,
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
  char **const fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  )
{
  struct BHFiller_S *const bhf = calloc(1,sizeof(*bhf));IsNull(bhf);
  /* physics */
  bhf->phys = phys;
  /* grid */
  Grid_T *const grid = phys->grid;
  bhf->grid = grid;
  /* method */
  sprintf(bhf->method,"%s",method);
  
  /* using Cheb_Tn and Ylm as the bases 
  // and extrapolate demanding C2 continuity. */
  if (strcmp_i(method,"ChebTn_Ylm_perfect_s2"))
  {
    const Uint NCoeffs = 6;/* number of coeffs in ChebTn expansion */
    const Uint lmax   = (Uint)Geti("filler_Ylm_expansion_lmax");
    const Uint Ntheta = Ntheta_Ylm(lmax);
    const Uint Nphi   = Nphi_Ylm(lmax);
    const Uint N      = Ntotal_Ylm(lmax);
    Uint npi;/* number of patches inside BH */
    Uint npo;/* number of patches outside BH */
    /* values of extrapolant function at the center of BH f(r=0) */
    char par[MAX_STR0];
    Uint f,nf,i;
    
    nf = 0;/* number of fields */
    while(fields_name[nf]) ++nf;
    
    /* alloc */
    bhf->nf = nf;
    bhf->fld= calloc(nf,sizeof(*bhf->fld));IsNull(bhf->fld);
    
    if (grid->kind == Grid_SplitCubedSpherical_BHBH ||
        grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH  ||
        grid->kind == Grid_CubedSpherical_BHNS
       )
    {
     /* set the method function */
     IF_sval("surface_type","perfect_s2")
       bhf->bhfiller = bhf_ChebTn_Ylm_pefect_S2_CS;
     else
       bhf->bhfiller = bhf_ChebTn_general_S2_CS;

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
      sprintf(par,"filler_r0_%s",fields_name[f]);
      bhf->fld[f]->f_r0 = Getd(par);
    }/* for (f = 0; f < nf ++f) */
    
  }/* if (strcmp_i(method,"ChebTn_Ylm_perfect_s2")) */
  else if (strcmp_i(method,"ChebTn_general_s2"))
  {
    const Uint NCoeffs = 6;/* number of coeffs in ChebTn expansion */
    Uint npi;/* number of patches inside BH */
    Uint npo;/* number of patches outside BH */
    /* values of extrapolant function at the center of BH f(r=0) */
    char par[MAX_STR0];
    Uint f,nf;
    
    nf = 0;/* number of fields */
    while(fields_name[nf]) ++nf;
    
    bhf->NCoeffs = NCoeffs;
    /* alloc */
    bhf->nf = nf;
    bhf->fld= calloc(nf,sizeof(*bhf->fld));IsNull(bhf->fld);
    
    if (grid->kind == Grid_SplitCubedSpherical_BHBH ||
        grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH  ||
        grid->kind == Grid_CubedSpherical_BHNS
       )
    {
     /* set the method function */
     bhf->bhfiller = bhf_ChebTn_general_S2_CS;
     
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
    }
    else
     Error0(NO_OPTION);
    
    /* collect names */
    collect_names(bhf,fields_name,nf);
    
    /* set values of field at r=0 */
    for (f = 0; f < nf; ++f)
    {
      sprintf(par,"filler_r0_%s",fields_name[f]);
      bhf->fld[f]->f_r0 = Getd(par);
    }/* for (f = 0; f < nf ++f) */
    
  }/* else if (strcmp_i(method,"ChebTn_general_s2")) */
  else if (strcmp_i(method,"expmr_C0_perfect_s2") ||
           strcmp_i(method,"r_expmr_C1_perfect_s2"))
  {
    Uint npi;/* number of patches inside BH */
    Uint npo;/* number of patches outside BH */
    Uint nf;
    
    nf = 0;/* number of fields */
    while(fields_name[nf]) ++nf;
    
    /* alloc */
    bhf->nf = nf;
    bhf->fld= calloc(nf,sizeof(*bhf->fld));IsNull(bhf->fld);
    
    if (grid->kind == Grid_SplitCubedSpherical_BHBH ||
        grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH  ||
        grid->kind == Grid_CubedSpherical_BHNS
       )
    {
     /* set the method function */
     bhf->bhfiller = bhf_f_df_ddf_perfect_s2_CS;
     
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
    }
    else
     Error0(NO_OPTION);
    
    if (strcmp_i(method,"expmr_C0_perfect_s2"))
    {
      bhf->extrap   = approx_expmr_C0;
      bhf->C2       = 0;
      bhf->C1       = 0;
    }
    else if (strcmp_i(method,"r_expmr_C1_perfect_s2")) 
    {
      bhf->extrap   = approx_r_expmr_C1;
      bhf->C2       = 0;
      bhf->C1       = 1;
    }
    else
      Error0(NO_OPTION);
      
    /* collect names */
    collect_names(bhf,fields_name,nf);
  }
  else
    Error0(NO_OPTION);
  
  return bhf;
}

/* free bhfiller struct */
static void bhf_free(struct BHFiller_S *const bhf)
{
  Uint i,f;
  
  if(!bhf)
    return;
  
  for (f = 0; f < bhf->nf; ++f)
  {
    for (i = 0; i < MAX_COEFFS; ++i)
    {
      Free(bhf->fld[f]->radial_coeffs[i]);
      Free(bhf->fld[f]->realYlm_coeffs[i]);
      Free(bhf->fld[f]->imagYlm_coeffs[i]);
    }
  }
  free_2d_mem(bhf->fld,bhf->nf);
  Free(bhf->patches_outBH);
  Free(bhf->patches_inBH);
  free(bhf);
}

/* ->: EXIT_SUCCESS if succeeds, otherwise an error code.
// method to fill BH is ChebTn_Ylm_perfect_s2 with the following extrapolant:
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
  Physics_T *const phys  = bhf->phys;
  const Uint NCoeffs = bhf->NCoeffs;
  const Uint npo     = bhf->npo;
  const Uint npi     = bhf->npi;
  const Uint nf     = bhf->nf;/* numebr of fields */
  const Uint lmax   = bhf->lmax;
  const Uint Ntheta = bhf->Ntheta;
  const Uint Nphi   = bhf->Nphi;
  const double rfill = Getd("perfect_S2_radius");
  const double rfill3= pow(rfill,3);
  const double BH_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};

  Uint p,fld;

  /* update all coeffs to avoid race condition */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = bhf->patches_outBH[p];
    Uint f = 0;

    /* make coeffs in  X and Y direction inside this patch */
    for (f = 0; f < nf; ++f)
    {
      int ii;
      
      /* add dfield if does not exist, field must exist already. */
      for (ii = 0; ii < 3; ++ii)
      {
        if (_Ind(bhf->fld[f]->df[ii]) < 0)
        {
          if(Verbose)
            printf(Pretty0"compute %s in %s\n",
                       bhf->fld[f]->df[ii],patch->name);
          bhf->fld[f]->did_add_df = 1;
          Field_T *df = add_field(bhf->fld[f]->df[ii],0,patch,NO);
          partial_derivative(df);
        }
      }
      /* add ddfield if does not exist, dfield must exist already. */
      for (ii = 0; ii < 6; ++ii)
      {
        if (_Ind(bhf->fld[f]->ddf[ii]) < 0)
        {
          if(Verbose)
            printf(Pretty0"compute %s in %s\n",
                       bhf->fld[f]->ddf[ii],patch->name);
          bhf->fld[f]->did_add_ddf = 1;
          Field_T *ddf = add_field(bhf->fld[f]->ddf[ii],0,patch,NO);
          partial_derivative(ddf);
        }
      }
      /* Note: partial derivatives modify coeffs thus it is made here */
      /* populate coeffs */
      make_coeffs_2d(patch->fields[Ind(bhf->fld[f]->f)],0,1);
      
      for (ii = 0; ii < 3; ++ii)
        make_coeffs_2d(patch->fields[Ind(bhf->fld[f]->df[ii])],0,1);
        
      for (ii = 0; ii < 6; ++ii)
        make_coeffs_2d(patch->fields[Ind(bhf->fld[f]->ddf[ii])],0,1);
    }
  }

  /* populating f, df/dr, d^2f/dr^2 at each (th,ph) points */
  printf(Pretty0"Populating extrapolation coefficients ...\n");
  fflush(stdout);
  OpenMP_1d_Pragma(omp parallel for)
  for (fld = 0; fld < nf; ++fld)
  {
    Uint i,j,_i,_j;
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
        Uint ij = IJ_Ylm(i,j,Nphi);
        Uint d1,d2;/* derivative */
        double X[3],x[3],_x[3],N[3];
        
        /* find patch for the given theta and phi */
        X[2] = 0.;
        find_XYZ_and_patch_of_theta_phi_CS
         (X,&patch,BH_center,theta,phi,bhf->patches_outBH,bhf->npo);
         
        /* r = rfill(sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^) */
        _x[0] = rfill*sin(theta)*cos(phi);
        _x[1] = rfill*sin(theta)*sin(phi);
        _x[2] = rfill*cos(theta);
        x[0]  = _x[0] + BH_center[0];
        x[1]  = _x[1] + BH_center[1];
        x[2]  = _x[2] + BH_center[2];
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
        /* f value at r = r1 and r = 0 */
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
    Uint nn    = patch->nn;
    Uint f,ijk;
    if(Verbose)
      printf(Pretty0"%s\n",patch->name),fflush(stdout);
    
    /* loop over all fields to be extrapolated */
    for (f = 0; f < nf; ++f)
    {
      Field_T *u = patch->fields[Ind(bhf->fld[f]->f)];
      empty_field(u);
      u->v      = alloc_double(patch->nn);
      double *v = u->v;
      double theta,phi,t;
      Uint i;
      
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

/* ->: EXIT_SUCCESS if succeeds, otherwise an error code.
// a various kind of extrapolations inside a perfect s2 BH.
// note: it has some assumptions which are only true in cubed spherical
// and only for perfect sphere. */
static int bhf_f_df_ddf_perfect_s2_CS(struct BHFiller_S *const bhf)
{
  Physics_T *const phys  = bhf->phys;
  const Uint npo = bhf->npo;
  const Uint npi = bhf->npi;
  const Uint nf  = bhf->nf;/* numebr of fields */
  const double rSurf = Getd("perfect_S2_radius");
  const double rSurf3= pow(rSurf,3);
  const double BH_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  Uint p;

  /* update all coeffs to avoid race condition */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = bhf->patches_outBH[p];
    Uint f = 0;

    /* make coeffs in  X and Y direction inside this patch */
    for (f = 0; f < nf; ++f)
    {
      int ii;
      
      /* add dfield if does not exist, field must exist already. */
      if (bhf->C1)
      for (ii = 0; ii < 3; ++ii)
      {
        if (_Ind(bhf->fld[f]->df[ii]) < 0)
        {
          if(Verbose)
            printf(Pretty0"compute %s in %s\n",
                       bhf->fld[f]->df[ii],patch->name);
          bhf->fld[f]->did_add_df = 1;
          Field_T *df = add_field(bhf->fld[f]->df[ii],0,patch,NO);
          partial_derivative(df);
        }
      }
      /* add ddfield if does not exist, dfield must exist already. */
      if (bhf->C2)
      for (ii = 0; ii < 6; ++ii)
      {
        if (_Ind(bhf->fld[f]->ddf[ii]) < 0)
        {
          if(Verbose)
            printf(Pretty0"compute %s in %s\n",
                       bhf->fld[f]->ddf[ii],patch->name);
          bhf->fld[f]->did_add_ddf = 1;
          Field_T *ddf = add_field(bhf->fld[f]->ddf[ii],0,patch,NO);
          partial_derivative(ddf);
        }
      }
      /* Note: partial derivatives modify coeffs thus it is made here */
      /* populate coeffs */
      make_coeffs_2d(patch->fields[Ind(bhf->fld[f]->f)],0,1);
      
      if (bhf->C1)
      for (ii = 0; ii < 3; ++ii)
        make_coeffs_2d(patch->fields[Ind(bhf->fld[f]->df[ii])],0,1);
      
      if (bhf->C2)
      for (ii = 0; ii < 6; ++ii)
        make_coeffs_2d(patch->fields[Ind(bhf->fld[f]->ddf[ii])],0,1);
    }
  }

  /* extrapolating inside */
  OpenMP_1d_Pragma(omp parallel for)
  for (p = 0; p < npi; p++)
  {
    Patch_T *ipatch = bhf->patches_inBH[p];
    Uint f = 0;
    
    for (f = 0; f < nf; ++f)
    {
     struct Demand_S demand[1] = {0};
     
     Field_T *field = 
       ipatch->fields[LookUpField_E(bhf->fld[f]->f,ipatch)];
     empty_field(field);
     field->v = alloc_double(ipatch->nn);
     
     for(Uint ijk = 0; ijk < ipatch->nn; ++ijk)
     {
      Patch_T *patch = 0;/* patch outside BH to be used for f,df,ddf */
      double th = 0,ph = 0;
      double X[3] = {0}, N[3] = {0}, x[3] = {0};
      double KD[2]      = {0,1};
      double df_dx[3]   = {0};
      double ddf_ddx[6] = {0};
      double ddfddr = 0,dfdr = 0,fr0 = 0;
      double _ddfddr[3] = {0,0,0};
      double r;
      Uint d1,d2;/* derivative */

      /* inside r,th,ph */
      x[0]= ipatch->node[ijk]->x[0]-BH_center[0];
      x[1]= ipatch->node[ijk]->x[1]-BH_center[1];
      x[2]= ipatch->node[ijk]->x[2]-BH_center[2];
      r   = sqrt(Pow2(x[0])+Pow2(x[1])+Pow2(x[2]));
      th = acos(x[2]/r);
      ph = arctan(x[1],x[0]);
      
      /* normal vector */
      N[0]  = sin(th)*cos(ph);
      N[1]  = sin(th)*sin(ph);
      N[2]  = cos(th); 

      X[2] = 0.;
      find_XYZ_and_patch_of_theta_phi_CS
        (X,&patch,BH_center,th,ph,bhf->patches_outBH,npo);
      
      /* find f(r0),df(r0),ddf(r0): */
      /* 2d interpolate on the surface */
      Interpolation_T *interp_s = init_interpolation();
      interp_s->XY_dir_flag = 1;
      interp_s->X = X[0];
      interp_s->Y = X[1];
      
      if (patch->coordsys == CubedSpherical)
      {
        interp_s->K = 0;
      }
      else 
        Error0(NO_OPTION);
       
      /* f value at r = r0 = rSurf*/
      interp_s->field = patch->fields[Ind(bhf->fld[f]->f)];
      plan_interpolation(interp_s);
      fr0 = execute_interpolation(interp_s);
       
      /* df/dx value */
      dfdr = 0.;
      if (bhf->C1)
      {
        for (d1 = 0; d1 < 3; d1++)
        {
          interp_s->field = patch->fields[Ind(bhf->fld[f]->df[d1])];
          plan_interpolation(interp_s);
          df_dx[d1] = execute_interpolation(interp_s);
        }
        /* df/dr */
        dfdr = (N[0]*df_dx[0]+N[1]*df_dx[1]+N[2]*df_dx[2]);
      }
      
      /* d^2f/dx^2 value */
      ddfddr = 0;
      if (bhf->C2)
      {
        for (d1 = 0; d1 < 3; d1++)
        {
          for (d2 = d1; d2 < 3; d2++)
          {
            interp_s->field = 
              patch->fields[Ind(bhf->fld[f]->ddf[IJsymm3(d1,d2)])];
            plan_interpolation(interp_s);
            ddf_ddx[IJsymm3(d1,d2)] = execute_interpolation(interp_s);
          }
        }
        _ddfddr[0] = _ddfddr[1] = _ddfddr[2] = 0;
        /* d^2f/dr^2 */
        for (int _i = 0; _i < 3; ++_i)
        {
          for (int _j = 0; _j < 3; ++_j)
          {
            _ddfddr[_i] += (KD[_i==_j]/rSurf - x[_i]*x[_j]/rSurf3)*df_dx[_j];
            _ddfddr[_i] += N[_j]*ddf_ddx[IJsymm3(_i,_j)];
          }
          ddfddr += _ddfddr[_i]*N[_i];
        }
      }
      
      /* free */
      free_interpolation(interp_s);
      
      /* bhf */
      demand->r     = r;
      demand->r0    = rSurf;
      demand->fr0   = fr0;
      demand->dfr0  = dfdr;
      demand->ddfr0 = ddfddr;
      
      field->v[ijk] = bhf->extrap(demand);
     }
    }
  }
  
  /* remove automatically added fields */
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = bhf->patches_outBH[p];
    Uint f = 0;

    for (f = 0; f < nf; ++f)
    {
      int ii;

      if (bhf->fld[f]->did_add_df)
      for (ii = 0; ii < 3; ++ii)
      {
         Field_T *df = patch->fields[Ind(bhf->fld[f]->df[ii])];
         REMOVE_FIELD(df);
      }

      if (bhf->fld[f]->did_add_ddf)
      for (ii = 0; ii < 6; ++ii)
      {
        Field_T *ddf = patch->fields[Ind(bhf->fld[f]->ddf[ii])];
        REMOVE_FIELD(ddf);
      }
    }
  }
 
  return EXIT_SUCCESS;
}

/* ->: EXIT_SUCCESS if succeeds, otherwise an error code.
// method to fill BH is ChebTn_general_s2 with the following extrapolant:
// ===============================================================
//
// f(r(t),th,ph) = C_{i}*ChebT_i(t)
//              
// where, t = 2*r/rfill-1.
// the coeffs a's are determinded by demaning the C2 continuity
// across the AH and the value of the function at r = 0.
// note: it has some assumptions which are only true in cubed spherical */
static int bhf_ChebTn_general_S2_CS(struct BHFiller_S *const bhf)
{
  Physics_T *const phys  = bhf->phys;
  const Uint NCoeffs = bhf->NCoeffs;
  const Uint npo     = bhf->npo;
  const Uint npi     = bhf->npi;
  const Uint nf      = bhf->nf;/* numebr of fields */
  const double Rmin  = Getd("filler_Rmin_cutoff");
  const double BH_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  Uint p,fld;

  /* update all coeffs to avoid race condition */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = bhf->patches_outBH[p];
    Uint f = 0;

    /* make coeffs in  X and Y direction inside this patch */
    for (f = 0; f < nf; ++f)
    {
      int ii;
      
      /* add dfield if does not exist, field must exist already. */
      for (ii = 0; ii < 3; ++ii)
      {
        if (_Ind(bhf->fld[f]->df[ii]) < 0)
        {
          if(Verbose)
            printf(Pretty0"compute %s in %s\n",
                       bhf->fld[f]->df[ii],patch->name);
          bhf->fld[f]->did_add_df = 1;
          Field_T *df = add_field(bhf->fld[f]->df[ii],0,patch,NO);
          partial_derivative(df);
        }
      }
      /* add ddfield if does not exist, dfield must exist already. */
      for (ii = 0; ii < 6; ++ii)
      {
        if (_Ind(bhf->fld[f]->ddf[ii]) < 0)
        {
          if(Verbose)
            printf(Pretty0"compute %s in %s\n",
                       bhf->fld[f]->ddf[ii],patch->name);
          bhf->fld[f]->did_add_ddf = 1;
          Field_T *ddf = add_field(bhf->fld[f]->ddf[ii],0,patch,NO);
          partial_derivative(ddf);
        }
      }
      /* Note: partial derivatives modify coeffs thus it is made here */
      /* populate coeffs */
      make_coeffs_2d(patch->fields[Ind(bhf->fld[f]->f)],0,1);
      
      for (ii = 0; ii < 3; ++ii)
        make_coeffs_2d(patch->fields[Ind(bhf->fld[f]->df[ii])],0,1);
        
      for (ii = 0; ii < 6; ++ii)
        make_coeffs_2d(patch->fields[Ind(bhf->fld[f]->ddf[ii])],0,1);
    }
  }

  /* fill */
  printf(Pretty0"Fill the hole ...\n");
  fflush(stdout);
  OpenMP_1d_Pragma(omp parallel for)
  for (fld = 0; fld < nf; ++fld)
  {
    const double KD[2] = {0,1};
    for (Uint ip = 0; ip < npi; ip++)
    {
      /* i for inside patch */
      Patch_T *ipatch = bhf->patches_inBH[ip];
      Field_T *u = ipatch->fields[LookUpField_E(bhf->fld[fld]->f,ipatch)];
      empty_field(u);
      u->v      = alloc_double(ipatch->nn);
      double *v = u->v;
      
      for(Uint ijk = 0; ijk < ipatch->nn; ++ijk)
      {
        /* o prefix stands for outside */
        Patch_T *patch = 0;/* this is the outside of BH patch */
        double theta,phi;
        double t,rSurf,rSurf3;
        double df_dx[3]   = {0};
        double ddf_ddx[6] = {0};
        double ddfddr = 0,dfdr = 0,fr1 = 0,fr0 = 0;
        double a[NCoeffs];
        double _ddfddr[3] = {0,0,0};
        double oX[3],ox[3],N[3];
        Uint d1,d2;/* derivative */
        Uint _i,_j;
        
        /* inside */
        double x=ipatch->node[ijk]->x[0]-BH_center[0];
        double y=ipatch->node[ijk]->x[1]-BH_center[1];
        double z=ipatch->node[ijk]->x[2]-BH_center[2];
        double r=sqrt(Pow2(x)+Pow2(y)+Pow2(z));
        
        /* if we don't want r smaller than Rmin */
        if (r < Rmin) continue;
        
        /* find outside patch for the given theta and phi */
        theta = acos(z/r);
        phi = arctan(y,x);
        oX[2] = 0.;
        find_XYZ_and_patch_of_theta_phi_CS
         (oX,&patch,BH_center,theta,phi,bhf->patches_outBH,bhf->npo);
        
        /* find r surface */
        assert(x_of_X(ox,oX,patch));
        ox[0] -= BH_center[0];
        ox[1] -= BH_center[1];
        ox[2] -= BH_center[2];
        rSurf  = sqrt(Pow2(ox[0])+Pow2(ox[1])+Pow2(ox[2]));
        rSurf3 = rSurf*Pow2(rSurf);
        
        /* r = rSurf(sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^) */
        ox[0] = rSurf*sin(theta)*cos(phi);
        ox[1] = rSurf*sin(theta)*sin(phi);
        ox[2] = rSurf*cos(theta);
        
        /* normal vector */
        N[0]  = sin(theta)*cos(phi);
        N[1]  = sin(theta)*sin(phi);
        N[2]  = cos(theta);
        
        /* 2d interpolate on the surface */
        Interpolation_T *interp_s = init_interpolation();
        interp_s->XY_dir_flag = 1;
        interp_s->X = oX[0];
        interp_s->Y = oX[1];
        interp_s->K = 0;
        /* f value at r = r1 and r = 0 */
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
            _ddfddr[_i] += (KD[_i==_j]/rSurf - ox[_i]*ox[_j]/rSurf3)*df_dx[_j];
            _ddfddr[_i] += N[_j]*ddf_ddx[IJsymm3(_i,_j)];
          }
          ddfddr += _ddfddr[_i]*N[_i];
        }
        
        /* find and set the coeffs */
        bh_bhf_ChebTn_extrapolate
        (a,fr0,fr1,dfdr,ddfddr,rSurf,NCoeffs);
        
        t = 2*r/rSurf-1;
        for (_i = 0; _i < NCoeffs; _i++)
          v[ijk] += a[_i]*Cheb_Tn((int)_i,t);
      }
    }
  }/* for (fld = 0; fld < nf ++fld) */
  
  return EXIT_SUCCESS;
}

/* collect names of the fields and their derivatives */
static void collect_names(struct BHFiller_S *const bhf,char **const fields_name,const Uint nf)
{
  const char *s = 0;
  Uint f,i,j;
    
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

/* ->: 1 if adds any patches inside the black hole, 0 otherwise.
// if the black hole excised from grid, this function add relevant
// patches to cover the inside and if already covered, it does nothing. 
// note: is also allocates and adds nodes. */
int bh_add_patch_inside_black_hole(Physics_T *const phys,
                                   const char *const region)
{
  Grid_T *const grid = phys->grid;
  assert(grid);
  
  /* check if there is already patches covering inside */
  for (Uint p = 0; p < grid->np; ++p) 
  {
    if (IsItCovering(grid->patch[p],region))
    {
      return 0;
    }
  }
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_BHBH ||
      grid->kind == Grid_SplitCubedSpherical_SBH
     )
  {
    populate_CS_patch_SplitCS(grid,region,phys->pos);
    populate_box_patch_SplitCS(grid,"central_box",phys->pos,region);
  }
  else
    Error0(NO_OPTION);
  
  alloc_nodes(mygrid(phys,region));
  make_nodes(mygrid(phys,region));
  make_JacobianT(mygrid(phys,region));
  
  return 1;
}


/* printing the specified quantities along a line for quality check */
void 
bh_interpolating_fields_on_a_line
  (
  Physics_T *const phys/* physics of interest */,
  const char *const sfields_name/* comma separated fields */,
  const char *const dir/* output directory */,
  const char *const stem_g/* if stem of a metric given => test det(g) > 0 */
  )
{
  FUNC_TIC
  
  AssureType(phys->ctype == BH);
  
  const double SmallDet = 1E-2;
  Grid_T *const grid = mygrid(phys,".*");
  char **const fields_name = 
    read_separated_items_in_string(sfields_name,','); 
  /* strcut for point where interpolate taken place */
  struct interpolation_points
  {
    double *x,*y,*z;/* (x,y,z) coords */
    double *X,*Y,*Z;/* (X,Y,Z) coords */
    Uint *patchn;/* patch number for each coord */
    Uint npoints;/* number of coords */
    int **f_index;/* field index for each patch and for each field
                  // ex: f_index[p][f] = for patch p and field f. */
  }pnt[1] = {0};
  /* the line eq.:
  // x = x_0 + t*mx
  // y = y_0 + t*my
  // z = z_0 + t*mz. */
  const double mx    = Getd("filler_test_print_1d_slop_x");
  const double my    = Getd("filler_test_print_1d_slop_y");
  const double mz    = Getd("filler_test_print_1d_slop_z");
  const double x_0   = Getd("filler_test_print_1d_init_x");
  const double y_0   = Getd("filler_test_print_1d_init_y");
  const double z_0   = Getd("filler_test_print_1d_init_z");
  const double Len   = Getd("filler_test_print_1d_length");
  const Uint npoints = (Uint)Geti("filler_test_print_1d_points");
  const double t     = Len/npoints;
  char fname[MAX_STR_LARGE] = {'\0'};
  double *interp_v = 0;
  FILE *file;
  Uint count_f;
  Uint p,f;
  
  if (npoints == 0)
    Error0("Bad parameter: number of points are zero.\n");
    
  /* populate points along y-axis since the objects are there */
  pnt->npoints = npoints;
  pnt->x       = alloc_double(npoints);
  pnt->y       = alloc_double(npoints);
  pnt->z       = alloc_double(npoints);
  pnt->X       = alloc_double(npoints);
  pnt->Y       = alloc_double(npoints);
  pnt->Z       = alloc_double(npoints);
  pnt->patchn  = calloc(npoints,sizeof(*pnt->patchn));
  IsNull(pnt->patchn);
  
  /* fill coords along the line. */
  for (p = 0; p < npoints; ++p)
  {
    pnt->x[p] = x_0+p*t*mx;
    pnt->y[p] = y_0+p*t*my;
    pnt->z[p] = z_0+p*t*mz;
  }
  
  /* find the corresponding X and patch */
  OpenMP_1d_Pragma(omp parallel for)
  for (p = 0; p < npoints; ++p)
  {
    Patch_T *patch = 0;
    double x[3],X[3];
    
    x[0] = pnt->x[p];
    x[1] = pnt->y[p];
    x[2] = pnt->z[p];
    
    patch = x_in_which_patch(x,grid->patch,grid->np);
    if(!patch || !X_of_x(X,x,patch))
      Error0("It could not find X(x,y,z) in any patch!\n");

    pnt->X[p]      = X[0];
    pnt->Y[p]      = X[1];
    pnt->Z[p]      = X[2];
    pnt->patchn[p] = patch->pn;
  }
  
  /* to avoid race condition between threads write all coeffs */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Uint fn = 0;
    
    while (fields_name[fn])
    {
      make_coeffs_3d(patch->fields[Ind(fields_name[fn])]);
      fn++;
    }
  }
  
  /* set f_index, note: it must be set right before interpolation
  // to make sure all fields are added already. */
  pnt->f_index = calloc(grid->np,sizeof(*pnt->f_index)); 
  IsNull(pnt->f_index);
  /* count f */
  count_f = 0;
  while(fields_name[count_f])
    ++count_f;
  
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    assert(patch->pn == p);
    
    pnt->f_index[p] = calloc(count_f,sizeof(*pnt->f_index[p]));
    IsNull(pnt->f_index[p]);
    
    f = 0;
    while(fields_name[f])
    {
      pnt->f_index[p][f] = Ind(fields_name[f]);
      ++f;
    }
  }
  
  interp_v = alloc_double(npoints);
  f = 0;
  while(fields_name[f])
  {
    /* interpolating each fields at the all given points */
    OpenMP_1d_Pragma(omp parallel for)
    for (p = 0; p < npoints; ++p)
    {
      Patch_T *patch  = grid->patch[pnt->patchn[p]];
      Interpolation_T *interp_s = init_interpolation();
      interp_s->field = patch->fields[pnt->f_index[patch->pn][f]];
      interp_s->XYZ_dir_flag = 1;
      interp_s->X = pnt->X[p];
      interp_s->Y = pnt->Y[p];
      interp_s->Z = pnt->Z[p];
      plan_interpolation(interp_s);
      interp_v[p] = execute_interpolation(interp_s);
      free_interpolation(interp_s);
    }
    
    /* write */
    sprintf(fname,"%s/%s_on_line_%0.1f_%0.1f_%0.1f.txt",dir,fields_name[f],mx,my,mz);
    file = Fopen(fname,"w");
    fprintf(file,"# fields value along the line:\n"
                 "# {\n"
                 "#   x = %g + t*%g,\n"
                 "#   y = %g + t*%g,\n"
                 "#   z = %g + t*%g.\n"
                 "# }\n"
                 "#\n"
                 "# x-coord y-coord z-coord %s\n",
                 x_0,mx,
                 y_0,my,
                 z_0,mz,
                 fields_name[f]);
    for (p = 0; p < npoints; ++p)
    {
      /* doc test */
      if (!isfinite(interp_v[p]))
      {
        printf("%s[%s]{x=(%g,%g,%g)&X=(%g,%g,%g)} = %g\n",
                fields_name[f],
                grid->patch[pnt->patchn[p]]->name,
                pnt->x[p],pnt->y[p],pnt->z[p],
                pnt->X[p],pnt->Y[p],pnt->Z[p],interp_v[p]);
      }
      fprintf(file,"%f %f %f %f\n",
                   pnt->x[p],pnt->y[p],pnt->z[p],interp_v[p]);
    }
    Fclose(file);
    f++;
  }/* while(fields_name[f]) */
  free_2d_mem(pnt->f_index,grid->np);
  pnt->f_index = 0;
  Free(interp_v);
  
  /* if check det(metric) requested and if fields_name contain stem_g */
  if (stem_g && strstr(sfields_name,stem_g))
  {
    /* to avoid race condition between threads write all coeffs */
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *patch = grid->patch[p];
      char gname[MAX_STR];
      
      make_coeffs_3d(patch->fields[Ind(PrefixIt(gname,stem_g,"D0D0"))]);
      make_coeffs_3d(patch->fields[Ind(PrefixIt(gname,stem_g,"D0D1"))]);
      make_coeffs_3d(patch->fields[Ind(PrefixIt(gname,stem_g,"D0D2"))]);
      make_coeffs_3d(patch->fields[Ind(PrefixIt(gname,stem_g,"D1D1"))]);
      make_coeffs_3d(patch->fields[Ind(PrefixIt(gname,stem_g,"D1D2"))]);
      make_coeffs_3d(patch->fields[Ind(PrefixIt(gname,stem_g,"D2D2"))]);
    }
  
    /* interpolating each fields at the all given points */
    OpenMP_1d_Pragma(omp parallel for)
    for (p = 0; p < npoints; ++p)
    {
      Patch_T *patch  = grid->patch[pnt->patchn[p]];
      char gname[MAX_STR];
      double gxx,gyy,gzz,gxy,gxz,gyz,detg;
      
      Interpolation_T *interp_s = init_interpolation();
      interp_s->XYZ_dir_flag = 1;
      interp_s->X = pnt->X[p];
      interp_s->Y = pnt->Y[p];
      interp_s->Z = pnt->Z[p];
      
      interp_s->field = patch->fields[Ind(PrefixIt(gname,stem_g,"D0D0"))];
      plan_interpolation(interp_s);
      gxx = execute_interpolation(interp_s);
      
      interp_s->field = patch->fields[Ind(PrefixIt(gname,stem_g,"D0D1"))];
      plan_interpolation(interp_s);
      gxy = execute_interpolation(interp_s);
      
      interp_s->field = patch->fields[Ind(PrefixIt(gname,stem_g,"D0D2"))];
      plan_interpolation(interp_s);
      gxz = execute_interpolation(interp_s);
      
      interp_s->field = patch->fields[Ind(PrefixIt(gname,stem_g,"D1D1"))];
      plan_interpolation(interp_s);
      gyy = execute_interpolation(interp_s);
      
      interp_s->field = patch->fields[Ind(PrefixIt(gname,stem_g,"D1D2"))];
      plan_interpolation(interp_s);
      gyz = execute_interpolation(interp_s);
      
      interp_s->field = patch->fields[Ind(PrefixIt(gname,stem_g,"D2D2"))];
      plan_interpolation(interp_s);
      gzz = execute_interpolation(interp_s);
      
      detg=(2.*gxy*gxz*gyz + gxx*gyy*gzz -
              gzz*gxy*gxy  - gyy*gxz*gxz -
              gxx*gyz*gyz);

      if(detg <= SmallDet)
      {
        printf("det(%s_ij(%g,%g,%g)) = %g\n",
             stem_g,pnt->x[p], pnt->y[p], pnt->z[p],detg);
      }
      free_interpolation(interp_s);
    }
  }/* if(?) */
  
  free_2d(fields_name);
  Free(pnt->x);
  Free(pnt->y);
  Free(pnt->z);
  Free(pnt->X);
  Free(pnt->Y);
  Free(pnt->Z);
  Free(pnt->patchn);
  free_2d_mem(pnt->f_index,grid->np);
  pnt->f_index = 0;
  
  FUNC_TOC
}

/* ->: f(r) = f(r0)*exp(Att*(r-r0)/r0), for r < r0.
// conditions: f be C^0 continues across the surface. */
static double approx_expmr_C0(struct Demand_S *const demand)
{
 const double r0    = demand->r0;
 const double fr0   = demand->fr0;
 const double r     = demand->r;
 const double Att   = Attenuate_Factor0;

 return fr0*exp(Att*fabs(r-r0)/r0);
}

/* ->: f(r) = (a+b*r)*exp(Att*(r-r0)/r0), for r < r0.
// conditions: f be C^1 continues across the surface. */
static double approx_r_expmr_C1(struct Demand_S *const demand)
{
 const double r0    = demand->r0;
 const double fr0   = demand->fr0;
 const double dfr0  = demand->dfr0;
 const double r     = demand->r;
 const double Att   = Attenuate_Factor0;
 double a,b;

 a = fr0 + Att*fr0 - dfr0*r0;
 b = dfr0 - (Att*fr0)/r0;
 
 return (a+b*r)*exp(Att*fabs(r-r0)/r0);
}

