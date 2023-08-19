/*
// Alireza Rashti
// November 2020
*/

/* various utilities related to surface of NS stars. */

#include "star_NS_surface.h"

/* make extrapolation weaker or stronger */
static const double Attenuate_Factor0 = 0.5;/* came from experiment */ 
static const double Attenuate_Factor1 = 1E-2;/* came from experiment */ 

/* ->: EXIT_SUCCESS if succeeds, otherwise an error code
// extrapolate given fields_name outside star surface.
// mostly used for matter fields. 
// more method see param "extrapolate_matter_fields". */
int 
star_NS_extrapolate
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  )
{
  FUNC_TIC
  
  if (phys->ctype != NS)
   Error0("Wrong physics!");
  
  int ret = -1;
  
  printf(Pretty0"method = %s\n",method);
  
  /* initialize */
  struct Extrap_S *const extrap = extrap_init(phys,fields_name,method);
  
  /* call main function to extrapolate */
  ret = extrap->fmain(extrap);
  
  /* free */
  extrap_free(extrap);
  
  FUNC_TOC
  
  return ret;
}


/* initialize the extrap struct */
static struct Extrap_S* 
extrap_init
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method or instruction 
                          // to be used for extrapolating */
  )
{
  struct Extrap_S *const extrap = calloc(1,sizeof(*extrap));IsNull(extrap);
  Grid_T *const grid = phys->grid;
  
  extrap->phys = phys;
  extrap->grid = grid;
  sprintf(extrap->method,"%s",method);
  
  if (strcmp_i(method,"inverse_r_expmAr") ||
      strcmp_i(method,"inverse_r_expmr"))
  {
    Uint nf,npo,npi;
    
    nf = 0;/* number of fields */
    while(fields_name[nf]) ++nf;
    
    /* alloc */
    extrap->nf = nf;
    extrap->fld= calloc(nf,sizeof(*extrap->fld));IsNull(extrap->fld);
    
    /* collect names */
    collect_names(extrap,fields_name,nf);
    
    if (grid->kind == Grid_SplitCubedSpherical_NSNS ||
        grid->kind == Grid_SplitCubedSpherical_BHBH ||
        grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH  ||
        grid->kind == Grid_SplitCubedSpherical_SNS  ||
        grid->kind == Grid_CubedSpherical_BHNS      ||
        grid->kind == Grid_CubedSpherical_NSNS
       )
    {
     /* this function finds field values and its derivative
     // on the surface with the known values of the field. */ 
     extrap->fmain = fmain_f_df_ddf_CS;
    }
    else
     Error0(NO_OPTION);
    
    /* set the function extrapolate and approximates the values */
    if (strcmp_i(method,"inverse_r_expmAr"))
      extrap->extrap = approx_inverse_r_expmAr;
    else if (strcmp_i(method,"inverse_r_expmr"))
      extrap->extrap = approx_inverse_r_expmr;
    else
      Error0(NO_OPTION);
    
    extrap->patches_in  = collect_patches(phys->grid,Ftype("NS_OB"),&npi);
    extrap->patches_out = collect_patches(phys->grid,Ftype("NS_around"),&npo);
     
    extrap->npo = npo;
    extrap->npi = npi;
  }
  else if (strcmp_i(method,"poly2")             ||
           strcmp_i(method,"exp2")              ||
           strcmp_i(method,"inverse_r2")        ||
           strcmp_i(method,"inverse_r2_expmr")  ||
           strcmp_i(method,"inverse_r2_expmAr"))
  {
    Uint nf,npo,npi;
    
    nf = 0;/* number of fields */
    while(fields_name[nf]) ++nf;
    
    /* alloc */
    extrap->nf = nf;
    extrap->fld= calloc(nf,sizeof(*extrap->fld));IsNull(extrap->fld);
    
    /* collect names */
    collect_names(extrap,fields_name,nf);
    
    if (grid->kind == Grid_SplitCubedSpherical_NSNS ||
        grid->kind == Grid_SplitCubedSpherical_BHBH ||
        grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH  ||
        grid->kind == Grid_SplitCubedSpherical_SNS  ||
        grid->kind == Grid_CubedSpherical_BHNS      ||
        grid->kind == Grid_CubedSpherical_NSNS
       )
    {
     /* this function finds field values and its derivative
     // on the surface with the known values of the field. */ 
     extrap->fmain = fmain_f_df_ddf_CS;
    }
    else
     Error0(NO_OPTION);
    
    /* set the function extrapolate and approximates the values */
    if (strcmp_i(method,"poly2"))
      extrap->extrap = approx_poly2;
    else if (strcmp_i(method,"exp2"))
      extrap->extrap = approx_exp2;
    else if (strcmp_i(method,"inverse_r2"))
      extrap->extrap = approx_inverse_r2;
    else if (strcmp_i(method,"inverse_r2_expmr"))
      extrap->extrap = approx_inverse_r2_expmr;
    else if (strcmp_i(method,"inverse_r2_expmAr"))
      extrap->extrap = approx_inverse_r2_expmAr;
    else
      Error0(NO_OPTION);
    
    extrap->patches_in  = collect_patches(phys->grid,Ftype("NS_OB"),&npi);
    extrap->patches_out = collect_patches(phys->grid,Ftype("NS_around"),&npo);
     
    extrap->npo = npo;
    extrap->npi = npi;
    
    extrap->C2  = 1;/* we need second oreder derivatives */
  }
  else if (strcmp_i(method,"expmr"))
  {
    Uint nf,npo,npi;
    
    nf = 0;/* number of fields */
    while(fields_name[nf]) ++nf;
    
    /* alloc */
    extrap->nf = nf;
    extrap->fld= calloc(nf,sizeof(*extrap->fld));IsNull(extrap->fld);
    
    /* collect names */
    collect_names(extrap,fields_name,nf);
    
    if (grid->kind == Grid_SplitCubedSpherical_NSNS ||
        grid->kind == Grid_SplitCubedSpherical_BHBH ||
        grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH  ||
        grid->kind == Grid_SplitCubedSpherical_SNS  ||
        grid->kind == Grid_CubedSpherical_BHNS      ||
        grid->kind == Grid_CubedSpherical_NSNS
       )
    {
     /* this function finds field values and its derivative
     // on the surface with the known values of the field. */ 
     extrap->fmain = extrapolate_expmr_C0_CS;
    }
    else
     Error0(NO_OPTION);
    
    extrap->patches_in  = collect_patches(phys->grid,Ftype("NS_OB"),&npi);
    extrap->patches_out = collect_patches(phys->grid,Ftype("NS_around"),&npo);
     
    extrap->npo = npo;
    extrap->npi = npi;
  }
  else
    Error0(NO_OPTION);
  
  return extrap;
}

/* free extrap struct */
static void extrap_free(struct Extrap_S *const extrap)
{
  if(!extrap)
    return;
  
  free_2d_mem(extrap->fld,extrap->nf);
  Free(extrap->patches_out);
  Free(extrap->patches_in);
  free(extrap);
}

/* ->: EXIT_SUCESS if succeeds, otherwise an error code.
// this function finds field values and its derivative
// on the surface with the known values of the field and then
// using them to extrapolate.
// requirement:f, df, ddf and the grid must be cubed spherical type.
// note: if it could not find df and ddf, it makes them and after
// is done with them, removes them. */
static int fmain_f_df_ddf_CS(struct Extrap_S *const extrap)
{
  Physics_T *const phys = extrap->phys;
  const Uint npo = extrap->npo;
  const Uint npi = extrap->npi;
  const Uint nf  = extrap->nf;/* numebr of fields */
  const double NS_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  Uint p;
  
  /* update all coeffs to avoid race condition */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npi; p++)
  {
    Patch_T *patch = extrap->patches_in[p];
    Uint f = 0;

    /* make coeffs in  X and Y direction inside this patch */
    for (f = 0; f < nf; ++f)
    {
      int ii;
      
      /* add dfield if does not exist, field must exist already. */
      for (ii = 0; ii < 3; ++ii)
      {
        if (_Ind(extrap->fld[f]->df[ii]) < 0)
        {
          if(Verbose)
            printf(Pretty0"compute %s in %s\n",
                       extrap->fld[f]->df[ii],patch->name);
          extrap->fld[f]->did_add_df = 1;
          Field_T *df = add_field(extrap->fld[f]->df[ii],0,patch,NO);
          partial_derivative(df);
        }
      }
      /* add ddfield if does not exist, dfield must exist already. */
      if (extrap->C2)
      for (ii = 0; ii < 6; ++ii)
      {
        if (_Ind(extrap->fld[f]->ddf[ii]) < 0)
        {
          if(Verbose)
            printf(Pretty0"compute %s in %s\n",
                       extrap->fld[f]->ddf[ii],patch->name);
          extrap->fld[f]->did_add_ddf = 1;
          Field_T *ddf = add_field(extrap->fld[f]->ddf[ii],0,patch,NO);
          partial_derivative(ddf);
        }
      }
      /* Note: partial derivatives modify coeffs thus it is made here */
      /* populate coeffs */
      make_coeffs_2d(patch->fields[Ind(extrap->fld[f]->f)],0,1);
      
      for (ii = 0; ii < 3; ++ii)
        make_coeffs_2d(patch->fields[Ind(extrap->fld[f]->df[ii])],0,1);
      
      if (extrap->C2)  
      for (ii = 0; ii < 6; ++ii)
        make_coeffs_2d(patch->fields[Ind(extrap->fld[f]->ddf[ii])],0,1);
    }
  }
  
  /* populating f, df/dr, d^2f/dr^2 at each (th,ph) points */
  OpenMP_1d_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *opatch = extrap->patches_out[p];
    Uint f = 0;
    
    for (f = 0; f < nf; ++f)
    {
     struct Demand_S demand[1] = {0};
     
     Field_T *field = 
       opatch->fields[LookUpField_E(extrap->fld[f]->f,opatch)];
     empty_field(field);
     field->v = alloc_double(opatch->nn);
     
     for(Uint ijk = 0; ijk < opatch->nn; ++ijk)
     {
      Patch_T *patch = 0;/* patch inside NS to be used for f,df,ddf */
      double th = 0,ph = 0;
      double X[3] = {0}, N[3] = {0}, x[3] = {0};
      double KD[2]      = {0,1};
      double df_dx[3]   = {0};
      double ddf_ddx[6] = {0};
      double ddfddr = 0,dfdr = 0,fr0 = 0;
      double _ddfddr[3] = {0,0,0};
      double rSurf,rSurf3,r;
      Uint d1,d2;/* derivative */

      /* outside r,th,ph */
      x[0]= opatch->node[ijk]->x[0]-NS_center[0];
      x[1]= opatch->node[ijk]->x[1]-NS_center[1];
      x[2]= opatch->node[ijk]->x[2]-NS_center[2];
      r   = sqrt(Pow2(x[0])+Pow2(x[1])+Pow2(x[2]));
      th = acos(x[2]/r);
      ph = arctan(x[1],x[0]);
      
      /* normal vector */
      N[0]  = sin(th)*cos(ph);
      N[1]  = sin(th)*sin(ph);
      N[2]  = cos(th); 

      /* find X and patch inside surface */
      if (opatch->coordsys == CubedSpherical)
      {
        X[2] = 1.;
        find_XYZ_and_patch_of_theta_phi_CS
          (X,&patch,NS_center,th,ph,extrap->patches_in,npi);
      }
      else
        Error0(NO_OPTION);
      
      /* find x and r at surface */
      assert(x_of_X(x,X,patch));
      x[0] -= NS_center[0];
      x[1] -= NS_center[1];
      x[2] -= NS_center[2];
      rSurf  = sqrt(Pow2(x[0])+Pow2(x[1])+Pow2(x[2]));
      rSurf3 = rSurf*Pow2(rSurf);
      
      /* find f(r0),df(r0),ddf(r0): */
      /* 2d interpolate on the surface */
      Interpolation_T *interp_s = init_interpolation();
      interp_s->XY_dir_flag = 1;
      interp_s->X = X[0];
      interp_s->Y = X[1];
      
      if (patch->coordsys == CubedSpherical)
      {
        interp_s->K = patch->n[2]-1;
      }
      else 
        Error0(NO_OPTION);
       
      /* f value at r = r0 = rSurf*/
      interp_s->field = patch->fields[Ind(extrap->fld[f]->f)];
      plan_interpolation(interp_s);
      fr0 = execute_interpolation(interp_s);
       
      /* df/dx value */
      for (d1 = 0; d1 < 3; d1++)
      {
        interp_s->field = patch->fields[Ind(extrap->fld[f]->df[d1])];
        plan_interpolation(interp_s);
        df_dx[d1] = execute_interpolation(interp_s);
      }
      
      /* d^2f/dx^2 value */
      if (extrap->C2)
      for (d1 = 0; d1 < 3; d1++)
      {
        for (d2 = d1; d2 < 3; d2++)
        {
          interp_s->field = 
            patch->fields[Ind(extrap->fld[f]->ddf[IJsymm3(d1,d2)])];
          plan_interpolation(interp_s);
          ddf_ddx[IJsymm3(d1,d2)] = execute_interpolation(interp_s);
        }
      }
      
      /* free */
      free_interpolation(interp_s);
       
      /* df/dr */
      dfdr = (N[0]*df_dx[0]+N[1]*df_dx[1]+N[2]*df_dx[2]);
      
      ddfddr = 0;
      if (extrap->C2)
      {
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
      /* extrap */
      demand->r     = r;
      demand->r0    = rSurf;
      demand->fr0   = fr0;
      demand->dfr0  = dfdr;
      demand->ddfr0 = ddfddr;
      
      field->v[ijk] = extrap->extrap(demand);
     }
    }
  }
  
  /* remove automatically added fields */
  for (p = 0; p < npi; p++)
  {
    Patch_T *patch = extrap->patches_in[p];
    Uint f = 0;

    for (f = 0; f < nf; ++f)
    {
      int ii;
      
      if (extrap->fld[f]->did_add_df)
      for (ii = 0; ii < 3; ++ii)
      {
         Field_T *df = patch->fields[Ind(extrap->fld[f]->df[ii])];
         REMOVE_FIELD(df);
      }
      
      if (extrap->fld[f]->did_add_ddf)
      for (ii = 0; ii < 6; ++ii)
      {
        Field_T *ddf = patch->fields[Ind(extrap->fld[f]->ddf[ii])];
        REMOVE_FIELD(ddf);
      }
    }
  }
  
  return EXIT_SUCCESS;
}

/* collect names of the fields and their derivatives */
static void collect_names(struct Extrap_S *const extrap,
                          const char **const fields_name,
                          const Uint nf)
{
  const char *s = 0;
  Uint f,i,j;
    
  for (f = 0; f < nf; ++f)
  {
    extrap->fld[f] = calloc(1,sizeof(*extrap->fld[f]));IsNull(extrap->fld[f]);
    /* names of fields and its derivatives */
    sprintf(extrap->fld[f]->f,"%s",fields_name[f]);
    if (fields_name[f][0] == '_')/* => _dgamma */
    {
      s = extrap->fld[f]->f+1;
      /* if it is indexed */
      if (regex_search(".+_(U|D)[[:digit:]]",s))
      {
        for (i = 0; i < 3; ++i)
        {
          sprintf(extrap->fld[f]->df[i],"_d%sD%u",s,i);
          for (j = i; j < 3; ++j)
            sprintf(extrap->fld[f]->ddf[IJsymm3(i,j)],"_dd%sD%uD%u",s,i,j);
        }
      }
      else/* not indexed */
      {
        for (i = 0; i < 3; ++i)
        {
          sprintf(extrap->fld[f]->df[0],"_d%s_D%u",s,i);
          for (j = i; j < 3; ++j)
            sprintf(extrap->fld[f]->ddf[IJsymm3(i,j)],"_dd%s_D%uD%u",s,i,j);
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
          sprintf(extrap->fld[f]->df[i],"d%sD%u",s,i);
          for (j = i; j < 3; ++j)
            sprintf(extrap->fld[f]->ddf[IJsymm3(i,j)],"dd%sD%uD%u",s,i,j);
        }
      }
      else/* not indexed */
      {
        for (i = 0; i < 3; ++i)
        {
          sprintf(extrap->fld[f]->df[i],"d%s_D%u",s,i);
          for (j = i; j < 3; ++j)
            sprintf(extrap->fld[f]->ddf[IJsymm3(i,j)],"dd%s_D%uD%u",s,i,j);
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
      printf("fld[%u] = %s\n",f,extrap->fld[f]->f);
      printf("d[%s] = (%s,%s,%s)\n",
          extrap->fld[f]->f,extrap->fld[f]->df[0],
          extrap->fld[f]->df[1],extrap->fld[f]->df[2]);
      for (i = 0; i < 3; ++i)
      {
        for (j = i; j < 3; ++j)
        {
          printf("dd[%s](%u,%u) = %s\n",
            extrap->fld[f]->f,i,j,extrap->fld[f]->ddf[IJsymm3(i,j)]);
        }
      }
    }
  }
}

/* ->: f(r) = a+b*r+c*r^2.
// conditions: f be C^2 continues across the surface. */
static double approx_poly2(struct Demand_S *const demand)
{
 const double r0 = demand->r0;
 const double fr0 = demand->fr0;
 const double dfr0  = demand->dfr0;
 const double ddfr0 = demand->ddfr0;
 const double r     = demand->r;
 double a,b,c;
 
 a = fr0 + (r0*(-2*dfr0 + ddfr0*r0))/2.;
 b = dfr0 - ddfr0*r0;
 c = 0.5*ddfr0;
 
 return a+b*r+c*Pow2(r);
}

/* ->: f(r) = a+b/r+c/r^2.
// conditions: f be C^2 continues across the surface. */
static double approx_inverse_r2(struct Demand_S *const demand)
{
 const double r0 = demand->r0;
 const double fr0 = demand->fr0;
 const double dfr0  = demand->dfr0;
 const double ddfr0 = demand->ddfr0;
 const double r     = demand->r;
 double a,b,c;
 
 a = fr0 + (r0*(4*dfr0 + ddfr0*r0))/2.;
 b = -(Pow2(r0)*(3*dfr0 + ddfr0*r0));
 c = (Pow3(r0)*(2*dfr0 + ddfr0*r0))/2.;
 
 return a+b/r+c/Pow2(r);
}

/* ->: f(r) = (a+b/r+c/r^2)*exp(-Att*r/r0).
// conditions: f be C^2 continues across the surface. */
static double approx_inverse_r2_expmr(struct Demand_S *const demand)
{
 const double r0 = demand->r0;
 const double fr0 = demand->fr0;
 const double dfr0  = demand->dfr0;
 const double ddfr0 = demand->ddfr0;
 const double r     = demand->r;
 const double Att   = Attenuate_Factor1;
 double a,b,c;
 
 a = (pow(M_E,Att)*((2 + 4*Att + Pow2(Att))*fr0 + 
       r0*(2*(2 + Att)*dfr0 + ddfr0*r0)))/2.;

 b = -(pow(M_E,Att)*r0*(3*Att*fr0 + Pow2(Att)*fr0 + 3*dfr0*r0 + 
       2*Att*dfr0*r0 + ddfr0*Pow2(r0)));

 c = (pow(M_E,Att)*Pow2(r0)*(Pow2(Att)*fr0 + 
       r0*(2*dfr0 + ddfr0*r0) + 2*Att*(fr0 + dfr0*r0)))/2.;
 
 return (a+b/r+c/Pow2(r))*exp(-Att*r/r0);
}

/* ->: f(r) = (a+b/r+c/r^2)*exp(-Att*(r-r0)), for r >= r0.
// conditions: f be C^2 continues across the surface. */
static double approx_inverse_r2_expmAr(struct Demand_S *const demand)
{
 const double r0 = demand->r0;
 const double fr0 = demand->fr0;
 const double dfr0  = demand->dfr0;
 const double ddfr0 = demand->ddfr0;
 const double r     = demand->r;
 const double Att   = Attenuate_Factor0;
 double a,b,c;
 
 a = (fr0*(2 + 4*Att*r0 + Pow2(Att)*Pow2(r0)) + 
     r0*(ddfr0*r0 + 2*dfr0*(2 + Att*r0)))/2.;
 b = -(Pow2(r0)*(3*Att*fr0 + ddfr0*r0 + Pow2(Att)*fr0*r0 + 
       dfr0*(3 + 2*Att*r0)));
 c = (Pow3(r0)*(2*Att*fr0 + ddfr0*r0 + Pow2(Att)*fr0*r0 + 
       2*dfr0*(1 + Att*r0)))/2.;
 
 return (a+b/r+c/Pow2(r))*exp(-Att*fabs(r-r0));
}

/* ->: f(r) = (a+b/r)*exp(-Att*(r-r0)), for r >= r0.
// conditions: f be C^1 continues across the surface. */
static double approx_inverse_r_expmAr(struct Demand_S *const demand)
{
 const double r0    = demand->r0;
 const double fr0   = demand->fr0;
 const double dfr0  = demand->dfr0;
 const double r     = demand->r;
 const double Att   = Attenuate_Factor0;
 double a,b;

 a = fr0 + dfr0*r0 + Att*fr0*r0;
 b = -((dfr0 + Att*fr0)*Pow2(r0));
 
 return (a+b/r)*exp(-Att*fabs(r-r0));
}

/* ->: f(r) = (a+b/r)*exp(-Att*(r/r0)), for r >= r0.
// conditions: f be C^1 continues across the surface. */
static double approx_inverse_r_expmr(struct Demand_S *const demand)
{
 const double r0    = demand->r0;
 const double fr0   = demand->fr0;
 const double dfr0  = demand->dfr0;
 const double r     = demand->r;
 const double Att   = Attenuate_Factor1;
 double a,b;

 a = pow(M_E,Att)*(fr0 + Att*fr0 + dfr0*r0);
 b = -(pow(M_E,Att)*r0*(Att*fr0 + dfr0*r0));
 
 return (a+b/r)*exp(-Att*r/r0);
}


/* ->: f(r) = a+b*exp(c*r)
// conditions: f be C^2 continues across the surface. */
static double approx_exp2(struct Demand_S *const demand)
{
 const double r0 = demand->r0;
 const double fr0 = demand->fr0;
 const double dfr0  = demand->dfr0;
 const double ddfr0 = demand->ddfr0;
 const double r     = demand->r;
 double a,b,c;

 if (EQL(ddfr0,0.) || EQL(dfr0,0.))
  Error0("Division by zero!");
  
 a = -(pow(dfr0,2)/ddfr0) + fr0;
 
 b = pow(dfr0,2)/(ddfr0*exp(ddfr0*r0/dfr0));
 
 c = ddfr0/dfr0;
 
 return a+b*exp(c*r);
 
}

/* using bisect root finder to find NS surface for cubed spherical type.
// find the NS surface on Ylm points 
// i.e (theta,phi) collocations are = (Legendre,EquiSpaced),
// using the fact that at the surface enthalpy = 1.
// it fills also NS feature at phys->grid_char:
//
// grid_char->params[phys->igc]->obj    
// grid_char->params[phys->igc]->dir    
// grid_char->params[phys->igc]->relClm 
// grid_char->params[phys->igc]->imgClm 
// grid_char->params[phys->igc]->r_min  
// grid_char->params[phys->igc]->r_max  
// grid_char->params[phys->igc]->lmax   
// grid_char->params[phys->igc]->occupied = 1; */
static void find_NS_surface_Ylm_bisect_CS(Physics_T *const phys)
{
  FUNC_TIC
  
  printf(Pretty0"%s\n",phys->stype);
  fflush(stdout);
  
  Grid_Char_T *grid_char = phys->grid_char;
  struct NS_surface_RootFinder_S par[1] = {0};
  const Uint lmax   = (Uint)Geti("surface_Ylm_max_l");
  const Uint Ntheta = Ntheta_Ylm(lmax);
  const Uint Nphi   = Nphi_Ylm(lmax);
  const Uint Ntot   = Ntotal_Ylm(lmax);
  const double Residual     = sqrt(Getd("RootFinder_Tolerance"));
  const double max_h_L2_res = Getd("enthalpy_allowed_residual");
  const double NS_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  const double incr_a = 0.01;
  const double incr_b = 0.01;
  const Uint Nincr    = 500;
  Uint Np_h,Np_ob;/* number of patches below: */
  Patch_T **patches_h  = collect_patches(phys->grid,Ftype("NS,NS_around"),&Np_h);
  Patch_T **patches_ob = collect_patches(phys->grid,Ftype("NS_OB"),&Np_ob);
  double h_L2_res = 0;
  double theta,phi;
  double *Rnew_NS = 0;/* new R for NS */
  double Max_R_NS = 0;/* maximum radius of NS */
  double Min_R_NS = DBL_MAX;/* minimum radius of NS */
  double *h_res   = 0;/* residual of h */
  double X[3],x[3],N[3];
  int NS_surface_finder_work_flg = 1;/* whether surface finder worked or not */
  Uint i,j;
  Uint l,m;
  
  /* populate root finder */
  Root_Finder_T *root = init_root_finder(1);
  root->type      = "Bisect_Single";
  root->tolerance = Getd("RootFinder_Tolerance");
  root->MaxIter   = (Uint)Geti("RootFinder_Iteration");
  root->params    = par;
  root->f[0]      = NS_surface_enthalpy_root_finder_eq;
  root->verbose   = strstr_i(Gets("RootFinder_verbose"),"yes");
    
  plan_root_finder(root);
  
  /* parameters for root finder */
  par->root_finder = root;
  
  /* initialize tables */
  init_Legendre_root_function();
  
  Rnew_NS = alloc_double(Ntot);
  h_res   = alloc_double(Ntot);
  /* for each points of Ylm find the surface of NS */
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      
      Patch_T *patch = 0;
      double y2[3] = {0};
      double *dr,a,b,Fa,Fb;
      Uint iincr;
      
      /* find patch and X,Y,Z at NS surface in which theta and phi take place */
      X[2] = 1.;
      find_XYZ_and_patch_of_theta_phi_CS(X,&patch,NS_center,
                                         theta,phi,patches_ob,Np_ob);
      /* finding x */
      assert(x_of_X(x,X,patch));
      
      /* r^ = sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^ */
      N[0]  = sin(theta)*cos(phi);
      N[1]  = sin(theta)*sin(phi);
      N[2]  = cos(theta);
      y2[0] = x[0]-NS_center[0];
      y2[1] = x[1]-NS_center[1];
      y2[2] = x[2]-NS_center[2];  
      
      par->hpatches = patches_h;
      par->Nph      = Np_h;
      par->x0[0]    = x[0];
      par->x0[1]    = x[1];
      par->x0[2]    = x[2];
      par->N        = N;
      /* set [a,b] for bisect */
      a  = 0.;
      Fa = NS_surface_enthalpy_root_finder_eq(par,&a);
      /* since we expect at some point h > 1 inside the star */
      iincr = 0;
      while (LSS(Fa,0.) && iincr < Nincr)
      {
        a -= incr_a;
        Fa = NS_surface_enthalpy_root_finder_eq(par,&a);
        iincr++;
      }
      if (Fa < 0.)
      {
        printf(Pretty0"Star's radius shrunk too much in %s!\n",
               patch->name);
        root->interrupt = 1;
        break;
      }
      
      b  = 0.;
      Fb = NS_surface_enthalpy_root_finder_eq(par,&b);
      /* since we expect at some point h < 1 outside the star */
      iincr = 0;
      while(GRT(Fb,0.) && iincr < Nincr)
      {
        b += incr_b;
        Fb = NS_surface_enthalpy_root_finder_eq(par,&b);
        iincr++;
      }
      if (Fb > 0.)
      {
        printf(Pretty0"Star's radius expanded too much in %s!\n",
               patch->name);
        root->interrupt = 1;
        break;
      }
      root->a_bisect  = a;
      root->b_bisect  = b;
      dr = execute_root_finder(root);
      h_res[IJ_Ylm(i,j,Nphi)] = root->residual;
      
      /* if NS surface finder interrupted */
      if (root->interrupt)
      {
        break;
      }
      /* if root finder is not OK for some reason */
      if (GRT(root->residual,Residual))
      {
        printf(". Root finder for NS surface at %s:\n.. ",patch->name);
        print_root_finder_exit_status(root);
        printf(".. Residual = %g\n",root->residual);
        //NS_surface_finder_work_flg = 0;/* since the residual is large don't update */
      }
        
      /*  new coords of R respect to the center of NS */
      y2[0] += N[0]*dr[0];
      y2[1] += N[1]*dr[0];
      y2[2] += N[2]*dr[0];
      Rnew_NS[IJ_Ylm(i,j,Nphi)] = root_square(3,y2,0);
      
      /* find the max NS radius */
      if (Rnew_NS[IJ_Ylm(i,j,Nphi)] > Max_R_NS)
        Max_R_NS = Rnew_NS[IJ_Ylm(i,j,Nphi)];
      /* find the min NS radius */
      if (Rnew_NS[IJ_Ylm(i,j,Nphi)] < Min_R_NS)
        Min_R_NS = Rnew_NS[IJ_Ylm(i,j,Nphi)];
    }/* end of for (j = 0; j < Nphi; ++j) */
    
    /* if NS surface finder interrupted */
    if (root->interrupt)
      break;
  }/* end of for (i = 0; i < Ntheta; ++i) */
  
  /* if NS surface finder interrupted */
  if (root->interrupt)
  {
    /* track the number of fails! */
    Uint mfail = (Uint)Geti("surface_max_fail");
    Uint nfail = (Uint)Geti("surface_num_fail");
    Seti("surface_num_fail",(int)(++nfail));
    
    if (nfail > mfail)
    {
      Error0("NS surface finder failed too many times!");
    }
    
    /* these are crucial for the next grid */
    printf(Pretty0"NS surface finder was interrupted.\n");
    
    Seti("did_NS_surface_finder_work?",0);
  
    assert(!grid_char->params[phys->igc]->occupied);
    grid_char->params[phys->igc]->obj   = phys->stype;
    grid_char->params[phys->igc]->dir   = phys->spos;
    grid_char->params[phys->igc]->r_max = Getd("max_radius");
    grid_char->params[phys->igc]->r_min = Getd("min_radius");
    grid_char->params[phys->igc]->lmax  = lmax;
    grid_char->params[phys->igc]->occupied = 1;
    
    free(Rnew_NS);
    free(h_res);
    free_root_finder(root);
    Free(patches_h);
    Free(patches_ob);
  
    FUNC_TOC
    return;
  }
  
  h_L2_res = L2_norm(Ntot,h_res,0);
  if (h_L2_res > max_h_L2_res)
    NS_surface_finder_work_flg = 0;/* since the residual is large don't update */
  else
    NS_surface_finder_work_flg = 1;
    
  /* making radius of NS parameter at each patch using Ylm interpolation */
  double *realClm = alloc_ClmYlm(lmax);
  double *imagClm = alloc_ClmYlm(lmax);
  
  /* calculating coeffs */
  get_Ylm_coeffs(realClm,imagClm,Rnew_NS,Ntheta,Nphi,lmax);

  assert(!grid_char->params[phys->igc]->occupied);
  grid_char->params[phys->igc]->obj    = phys->stype;
  grid_char->params[phys->igc]->dir    = phys->spos;
  grid_char->params[phys->igc]->relClm = realClm;
  grid_char->params[phys->igc]->imgClm = imagClm;
  grid_char->params[phys->igc]->r_min  = Min_R_NS;
  grid_char->params[phys->igc]->r_max  = Max_R_NS;
  grid_char->params[phys->igc]->lmax   = lmax;
  grid_char->params[phys->igc]->occupied = 1;
    
  /* printing */
  printf(Pretty0"Max NS radius           = %e\n",Max_R_NS);
  printf(Pretty0"Min NS radius           = %e\n",Min_R_NS);
  printf(Pretty0"L2 norm of enthalpy     = %e\n",h_L2_res);
  printf(Pretty0"Mass shedding indicator = %e\n",
                      star_NS_mass_shedding_indicator(phys));

  Setd("enthalpy_L2_residual",h_L2_res);
  
  l = lmax;
  for (m = 0; m <= l; ++m)
  {
    Uint lm = lm2n(l,m);
    printf(Pretty0"Truncation error [Real(C[%u][%u])] = %e\n",l,m,realClm[lm]);
    printf(Pretty0"Truncation error [Imag(C[%u][%u])] = %e\n",l,m,imagClm[lm]);
  }
  
  /* if some day you wanna filter Clm's */
  if (0)
  {
    const double e = 0.1;
    for (l = 0; l <= lmax; ++l)
      for (m = 0; m <= l; ++m)
      {
        Uint lm = lm2n(l,m);
        realClm[lm] /= (1+e*Pow2(l)*Pow2(l+1));
        imagClm[lm] /= (1+e*Pow2(l)*Pow2(l+1));
      }
  }
  
  free(Rnew_NS);
  free(h_res);
  free_root_finder(root);
  
  Seti("did_NS_surface_finder_work?",NS_surface_finder_work_flg);
  
  Setd("max_radius",Max_R_NS);
  Setd("min_radius",Min_R_NS);
  
  Free(patches_h);
  Free(patches_ob);
  
  FUNC_TOC
}

/* find r such that f(h(r)) = h(r)-1 = 0.
// the root finder moving along the r and it seeks for h = 1. */
static double NS_surface_enthalpy_root_finder_eq(void *params,const double *const x)
{
  const struct NS_surface_RootFinder_S *const pars = params;
  const double dx        = x[0];
  const double *const x0 = pars->x0;
  const double *const N  = pars->N;
  const double y[3]      = {x0[0]+dx*N[0],x0[1]+dx*N[1],x0[2]+dx*N[2]};
  Patch_T *patch = 0;
  double X[3],h;
  int h_ind;
  
  /* find enthalpy at the (X,Y,Z) */
  patch = x_in_which_patch(y,pars->hpatches,pars->Nph);
  if (!patch)
  {
    /* try */
    patch = x_in_which_patch_force(y,pars->hpatches,pars->Nph,X);
    if (!patch)/* fails! */
    {
      Root_Finder_T *root_finder = pars->root_finder;
      root_finder->interrupt = 1;
      return DBL_MAX;
    }
  }
  else
    assert(X_of_x(X,y,patch));
    
  h_ind = Ind("enthalpy");
  if (!patch->fields[h_ind]->v)/* if there is no enthalpy defined in the patch */
    return -1.;
    
  Interpolation_T *interp_h = init_interpolation();
  interp_h->field = patch->fields[h_ind];
  interp_h->XYZ_dir_flag = 1;
  interp_h->X            = X[0];
  interp_h->Y            = X[1];
  interp_h->Z            = X[2];
  plan_interpolation(interp_h);
  h = execute_interpolation(interp_h);/* enthalpy */
  free_interpolation(interp_h);
  
  return h-1;
}

/* denthalpy(r)/dr for NS surface root finder */
static double NS_surface_denthalpy_dr_root_finder(void *params,const double *const x,const Uint dir)
{
  assert(dir == 0);
  const struct NS_surface_RootFinder_S *const pars = params;
  const double dx        = x[0];
  const double *const x0 = pars->x0;
  const double *const N  = pars->N;
  const double y[3]      = {x0[0]+dx*N[0],x0[1]+dx*N[1],x0[2]+dx*N[2]};
  Patch_T *patch = 0;
  double X[3],dh_dx,dh_dy,dh_dz;
  int dh_dx_ind,dh_dy_ind,dh_dz_ind;
  Interpolation_T *interp_dh_dx = 0,
                  *interp_dh_dy = 0,
                  *interp_dh_dz = 0;
  
    
    
  patch = x_in_which_patch(y,pars->hpatches,pars->Nph);
  if (!patch)
  {
    Root_Finder_T *root_finder = pars->root_finder;
    root_finder->interrupt = 1;
    return DBL_MAX;
  }
 
  /* find denthalpy/dr at the (X,Y,Z): */
  assert(X_of_x(X,y,patch));
  dh_dx_ind = _Ind("denthalpy_D0");
  if (dh_dx_ind < 0)/* if there is no enthalpy defined in the patch */
    return 1;
  interp_dh_dx = init_interpolation();
  interp_dh_dx->field = patch->fields[dh_dx_ind];
  interp_dh_dx->X = X[0];
  interp_dh_dx->Y = X[1];
  interp_dh_dx->Z = X[2];
  interp_dh_dx->XYZ_dir_flag = 1;
  plan_interpolation(interp_dh_dx);
  dh_dx = execute_interpolation(interp_dh_dx);
  free_interpolation(interp_dh_dx);
  
  dh_dy_ind = _Ind("denthalpy_D1");
  if (dh_dy_ind < 0)/* if there is no enthalpy defined in the patch */
    return 1;
  interp_dh_dy = init_interpolation();
  interp_dh_dy->field = patch->fields[dh_dy_ind];
  interp_dh_dy->X = X[0];
  interp_dh_dy->Y = X[1];
  interp_dh_dy->Z = X[2];
  interp_dh_dy->XYZ_dir_flag = 1;
  plan_interpolation(interp_dh_dy);
  dh_dy = execute_interpolation(interp_dh_dy);
  free_interpolation(interp_dh_dy);
  
  dh_dz_ind = _Ind("denthalpy_D2");
  if (dh_dz_ind < 0)/* if there is no enthalpy defined in the patch */
    return 1;
  interp_dh_dz = init_interpolation();
  interp_dh_dz->field = patch->fields[dh_dz_ind];
  interp_dh_dz->X = X[0];
  interp_dh_dz->Y = X[1];
  interp_dh_dz->Z = X[2];
  interp_dh_dz->XYZ_dir_flag = 1;
  plan_interpolation(interp_dh_dz);
  dh_dz = execute_interpolation(interp_dh_dz);
  free_interpolation(interp_dh_dz);
  
  /* Grad h . r^ = dh/dr */
  return N[0]*dh_dx+N[1]*dh_dy+N[2]*dh_dz;
}

/* calculating mass shedding indicator:
// Chi = {d(ln h)/dr|equator}/{d(ln h)/dr|pole} */
double star_NS_mass_shedding_indicator(Physics_T *const phys)
{
  if (phys->grid->kind == Grid_SplitCubedSpherical_SNS  ||
      phys->grid->kind == Grid_SplitCubedSpherical_NSNS ||
      phys->grid->kind == Grid_SplitCubedSpherical_BHNS
     )
  {
  
    const double NS_center[3] = {Getd("center_x"),
                                 Getd("center_y"),
                                 Getd("center_z")};
    Patch_T *patch    = 0;
    Patch_T **patches = 0;
    Interpolation_T *interp_s = init_interpolation();
    double dh_eq[3] = {0},dh_pole[3] = {0};/* h derivatives */
    double h_eq,h_pole;/* h values */
    double dr_dlnh_pole,dr_dlnh_eq;
    double N[3] = {0},X[3] = {0},x[3] = {0};
    double r,theta,phi;
    char regex[99] = {'\0'};
    const char *side = 0;
    Uint Np;
    
    /* opposite */
    if      (phys->pos == LEFT)   side = "right";
    else if (phys->pos == RIGHT)  side = "left";
    else if (phys->pos == CENTER) side = "left";// since assumed symm. doesn't matter which side
    else                          Error0(NO_OPTION);
    
    /* approx. equator value  (X,Y,Z)=(0,0,1) */
    sprintf(regex,".*%s_%s_%s.*",phys->spos,phys->stype,side);
    patches = collect_patches_regex(phys->grid,regex,&Np);
    X[0] = 0;
    X[1] = 0;
    X[2] = 1;
    patch = X_in_which_patch(X,patches,Np);
    Free(patches);
    assert(patch);
    x_of_X(x,X,patch);
    x[0] -= NS_center[0];
    x[1] -= NS_center[1];
    x[2] -= NS_center[2];
    r     = root_square(3,x,0);
    theta = acos(x[2]/r);
    phi   = arctan(x[1],x[0]);
    N[0]  = sin(theta)*cos(phi);
    N[1]  = sin(theta)*sin(phi);
    N[2]  = cos(theta);
    interp_s->X = 0;
    interp_s->Y = 0;
    interp_s->K = patch->n[2]-1;
    interp_s->XY_dir_flag = 1;
    
    /* derivatives */
    interp_s->field = patch->fields[Ind("denthalpy_D0")];
    plan_interpolation(interp_s);
    dh_eq[0] = execute_interpolation(interp_s);
    
    interp_s->field = patch->fields[Ind("denthalpy_D1")];
    plan_interpolation(interp_s);
    dh_eq[1] = execute_interpolation(interp_s);
    
    interp_s->field = patch->fields[Ind("denthalpy_D2")];
    plan_interpolation(interp_s);
    dh_eq[2] = execute_interpolation(interp_s);
    
    /* value */
    interp_s->field = patch->fields[Ind("enthalpy")];
    plan_interpolation(interp_s);
    h_eq = execute_interpolation(interp_s);
    
    dr_dlnh_eq = (N[0]*dh_eq[0]+N[1]*dh_eq[1]+N[2]*dh_eq[2])/h_eq;
    
    /* approx. north pole value (X,Y,Z)=(0,0,1) */
    sprintf(regex,".*%s_%s_up.*",phys->spos,phys->stype);
    patches = collect_patches_regex(phys->grid,regex,&Np);
    X[0] = 0;
    X[1] = 0;
    X[2] = 1;
    patch = X_in_which_patch(X,patches,Np);
    Free(patches);
    assert(patch);
    x_of_X(x,X,patch);
    x[0] -= NS_center[0];
    x[1] -= NS_center[1];
    x[2] -= NS_center[2];
    r     = root_square(3,x,0);
    theta = acos(x[2]/r);
    phi   = arctan(x[1],x[0]);
    N[0]  = sin(theta)*cos(phi);
    N[1]  = sin(theta)*sin(phi);
    N[2]  = cos(theta);
    interp_s->X = 0;
    interp_s->Y = 0;
    interp_s->K = patch->n[2]-1;
    interp_s->XY_dir_flag = 1;
    
    /* derivatives */
    interp_s->field = patch->fields[Ind("denthalpy_D0")];
    plan_interpolation(interp_s);
    dh_pole[0] = execute_interpolation(interp_s);
    
    interp_s->field = patch->fields[Ind("denthalpy_D1")];
    plan_interpolation(interp_s);
    dh_pole[1] = execute_interpolation(interp_s);
    
    interp_s->field = patch->fields[Ind("denthalpy_D2")];
    plan_interpolation(interp_s);
    dh_pole[2] = execute_interpolation(interp_s);
    
    /* value */
    interp_s->field = patch->fields[Ind("enthalpy")];
    plan_interpolation(interp_s);
    h_pole = execute_interpolation(interp_s);
    dr_dlnh_pole = (N[0]*dh_pole[0]+N[1]*dh_pole[1]+N[2]*dh_pole[2])/h_pole;
    
    free_interpolation(interp_s);
    
    return dr_dlnh_eq/dr_dlnh_pole;
  }
  else
    Error0(NO_OPTION);
   
  return 0; 
}

/* find NS surface */
int star_NS_find_star_surface(Physics_T *const phys)
{
  IF_sval("surface_type","topology_s2")
  {
    IF_sval("surface_finder","bisection")
      find_NS_surface_Ylm_bisect_CS(phys);
    else
      Error0(NO_OPTION);
  }
  else IF_sval("surface_type","perfect_s2")
  {
    find_NS_surface_perfect_s2(phys);
  }
  else
      Error0(NO_OPTION);
  
  return EXIT_SUCCESS;
}

/* find NS surface and then set grid characteristic for perfect S2 */
static void find_NS_surface_perfect_s2(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_Char_T *grid_char = phys->grid_char;
  const Uint lmax   = (Uint)Geti("surface_Ylm_max_l");
  const Uint Ntheta = Ntheta_Ylm(lmax);
  const Uint Nphi   = Nphi_Ylm(lmax);
  const Uint Ntot   = Ntotal_Ylm(lmax);
  const double R_NS = Getd("perfect_S2_radius");
  double *rns = alloc_double(Ntot);/* surface function r = r(th,ph). */
  double *reClm_rns = alloc_ClmYlm(lmax),
         *imClm_rns = alloc_ClmYlm(lmax);
  Uint ij;
  
  init_Legendre_root_function();
  for (ij = 0; ij < Ntot; ++ij)
  {
    rns[ij] = R_NS;
  }
  /* calculating coeffs */
  get_Ylm_coeffs(reClm_rns,imClm_rns,rns,Ntheta,Nphi,lmax);
  
  assert(!grid_char->params[phys->igc]->occupied);
  grid_char->params[phys->igc]->obj    = phys->stype;
  grid_char->params[phys->igc]->dir    = phys->spos;
  grid_char->params[phys->igc]->relClm = reClm_rns;
  grid_char->params[phys->igc]->imgClm = imClm_rns;
  grid_char->params[phys->igc]->r_min  = Getd("min_radius");
  grid_char->params[phys->igc]->r_max  = Getd("max_radius");
  grid_char->params[phys->igc]->lmax   = lmax;
  grid_char->params[phys->igc]->occupied = 1;
  Seti("did_NS_surface_finder_work?",1);
  
  Free(rns);
  
  FUNC_TOC
}

/* use TOV star (pefect S2) to start off a star parameters 
// and domain shape etc. */
void star_start_off_TOV(Physics_T *const phys)
{
  TOV_T *tov = TOV_init();
  tov->phys  = phys;
  tov->bar_m = Getd("baryonic_mass");
  tov = TOV_solution(tov);
  
  const double ns_R = tov->rbar[tov->N-1];
  
  Setd("perfect_S2_radius",ns_R);
  Setd("min_radius",ns_R);
  Setd("max_radius",ns_R);
  
  TOV_free(tov);
}

/* ->: success int.
// decrease given field(s) exponentially from the NS surface 
//  it is C^0 continues. */
static int extrapolate_expmr_C0_CS(struct Extrap_S *const extrap)
{
  Physics_T *const phys = extrap->phys;
  const Uint npo = extrap->npo;
  const Uint npi = extrap->npi;
  const Uint nf  = extrap->nf;/* numebr of fields */
  const double NS_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  const double att          = 0.1;/* exp(-att*(r2-r1)) */
  Uint p;

  /* update all coeffs to avoid race condition */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npi; p++)
  {
    Patch_T *patch = extrap->patches_in[p];
    Uint f = 0;

    /* make coeffs in  X and Y direction inside this patch */
    for (f = 0; f < nf; ++f)
    {
      /* must have the field */
      make_coeffs_2d(patch->fields[Ind(extrap->fld[f]->f)],0,1);
    }
  }
  
  /* populating f */
  OpenMP_1d_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *opatch = extrap->patches_out[p];
    Uint f = 0;
    
    for (f = 0; f < nf; ++f)
    {
      Field_T *field = 
        opatch->fields[LookUpField_E(extrap->fld[f]->f,opatch)];
      empty_field(field);
      field->v = alloc_double(opatch->nn);
      
      for(Uint ijk = 0; ijk < opatch->nn; ++ijk)
      {
        Patch_T *patch = 0;/* patch inside NS to be used for f0 */
        double th = 0,ph = 0;
        double X[3] = {0}, x[3] = {0};
        double fr0 = 0;
        double rSurf,r;

        /* outside r,th,ph */
        x[0]= opatch->node[ijk]->x[0]-NS_center[0];
        x[1]= opatch->node[ijk]->x[1]-NS_center[1];
        x[2]= opatch->node[ijk]->x[2]-NS_center[2];
        r   = sqrt(Pow2(x[0])+Pow2(x[1])+Pow2(x[2]));
        th = acos(x[2]/r);
        ph = arctan(x[1],x[0]);
        
        /* find X and patch inside surface */
        if (opatch->coordsys == CubedSpherical)
        {
          X[2] = 1.;
          find_XYZ_and_patch_of_theta_phi_CS
            (X,&patch,NS_center,th,ph,extrap->patches_in,npi);
        }
        else
          Error0(NO_OPTION);
        
        /* find x and r at surface */
        assert(x_of_X(x,X,patch));
        x[0] -= NS_center[0];
        x[1] -= NS_center[1];
        x[2] -= NS_center[2];
        rSurf  = sqrt(Pow2(x[0])+Pow2(x[1])+Pow2(x[2]));
        
        /* find f(r0) */
        /* 2d interpolate on the surface */
        Interpolation_T *interp_s = init_interpolation();
        interp_s->XY_dir_flag = 1;
        interp_s->X = X[0];
        interp_s->Y = X[1];
        
        if (patch->coordsys == CubedSpherical)
        {
          interp_s->K = patch->n[2]-1;
        }
        else 
          Error0(NO_OPTION);
         
        /* f value at r = r0 = rSurf*/
        interp_s->field = patch->fields[Ind(extrap->fld[f]->f)];
        plan_interpolation(interp_s);
        fr0 = execute_interpolation(interp_s);
        free_interpolation(interp_s);
         
        /* exp extrapolate */
        field->v[ijk] = fr0*exp(-att*(r-rSurf));
      }
    }
  }
  
  return EXIT_SUCCESS;
}

