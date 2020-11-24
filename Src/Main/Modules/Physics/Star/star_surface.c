/*
// Alireza Rashti
// November 2020
*/

/* various utilities related to surface of stars. */

#include "star_surface.h"

/* ->: EXIT_SUCCESS if succeeds, otherwise an error code
// extrapolate given fields_name outside star surface.
// mostly used for matter fields.
//
// method:
// ========
// note: r is coordinate distance to the center of patch covers the star
//
// poly2: 
// f(r) = a+b*r+c*r^2.
// conditions: f be C^2 continues across the surface.
//
// exp2:
// f(r) = a+b*exp(c*r).
// conditions: f be C^2 continues across the surface.
//
// */
int 
star_extrapolate
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  )
{
  FUNC_TIC
  printf(Pretty0"method = %s\n",method);
  
  if (phys->type != NS)
   Error0("Wrong physics!");
  
  int ret = -1;
  
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
  
  extrap->grid = grid;
  sprintf(extrap->method,"%s",method);
  
  if (strcmp_i(method,"poly2") ||
      strcmp_i(method,"exp2")
     )
  {
    unsigned nf,npo,npi;
    
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
        grid->kind == Grid_SplitCubedSpherical_BH   ||
        grid->kind == Grid_SplitCubedSpherical_NS   ||
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
    if (strcmp_i(method,"poly2"))      extrap->extrap = approx_poly2;
    else if (strcmp_i(method,"exp2"))  extrap->extrap = approx_exp2;
    else Error0(NO_OPTION);
    
    extrap->patches_in = collect_patches(phys->grid,Ftype("NS_OB"),&npi);
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
  _free(extrap->patches_out);
  _free(extrap->patches_in);
  free(extrap);
}

/* ->: EXIT_SUCESS if succeeds, otherwise an error code.
// this function finds field values and its derivative
// on the surface with the known values of the field. 
// requirement:f, df, ddf and the grid must be cubed spherical type. */
static int fmain_f_df_ddf_CS(struct Extrap_S *const extrap)
{
  const unsigned npo = extrap->npo;
  const unsigned npi = extrap->npi;
  const unsigned nf  = extrap->nf;/* numebr of fields */
  unsigned p;
  
  /* update all coeffs to avoid race condition */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npi; p++)
  {
    Patch_T *patch = extrap->patches_in[p];
    unsigned f = 0;

    /* make coeffs in  X and Y direction inside this patch */
    for (f = 0; f < nf; ++f)
    {
      int ii;
      
      /* must have the field */
      make_coeffs_2d(patch->pool[Ind(extrap->fld[f]->f)],0,1);
      
      /* must have dfield */
      for (ii = 0; ii < 3; ++ii)
      {
        int indxf = _Ind(extrap->fld[f]->df[ii]);
        if (indxf >= 0)
         make_coeffs_2d(patch->pool[indxf],0,1);
        else
        Error0("This kind of extrapolate depends on "
                "first and second order derivatives.\n"
                "If there in not any, please use or add another option.");
      }
      
      /* must have ddfield */
      for (ii = 0; ii < 6; ++ii)
      {
        int indxf = _Ind(extrap->fld[f]->ddf[ii]);
        if (indxf >= 0)
          make_coeffs_2d(patch->pool[indxf],0,1);
        else
         Error0("This kind of extrapolate depends on "
                "first and second order derivatives.\n"
                "If there in not any, please use or add another option.");
      }
      
      ++f;
    }
  }
  
  /* populating f, df/dr, d^2f/dr^2 at each (th,ph) points */
  OpenMP_1d_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = extrap->patches_out[p];
    unsigned f = 0;
    
    for (f = 0; f < nf; ++f)
    {
     struct Demand_S demand[1] = {0};
     
     Field_T *field = patch->pool[Ind(extrap->fld[f]->f)];
     empty_field(field);
     field->v = alloc_double(patch->nn);
     
     forall_ijk
     {
      Patch_T *patchp = 0;/* patch prime to be used for f,df,ddf */
      double th = 0,ph = 0;
      double X[3] = {0}, N[3] = {0}, x[3] = {0};
      double KD[2]      = {0,1};
      double df_dx[3]   = {0};
      double ddf_ddx[6] = {0};
      double ddfddr = 0,dfdr = 0,fr0 = 0;
      double _ddfddr[3] = {0,0,0};
      double rSurf,rSurf3,r;
      unsigned d1,d2;/* derivative */

      /* find th and ph and X */
      if (patch->coordsys == CubedSpherical)
      {
       find_theta_phi_of_XYZ_CS(&th,&ph,patch->node[ijk]->X,
             patch->CoordSysInfo->CubedSphericalCoord->side);
       /* find xp in patch_in */
       X[2] = 1.;
       find_XYZ_and_patch_of_theta_phi_CS(X,&patchp,th,ph,
                                          extrap->patches_in,npi);
      }
      else
       Error0(NO_OPTION);
      
      
      /* find r */
      assert(x_of_X(x,X,patch));
      x[0] -= patchp->c[0];
      x[1] -= patchp->c[1];
      x[2] -= patchp->c[2];
      r     = sqrt(Pow2(x[0])+Pow2(x[1])+Pow2(x[2]));
      
      /* normal vector */
      N[0]  = sin(th)*cos(ph);
      N[1]  = sin(th)*sin(ph);
      N[2]  = cos(th); 
      
      /* find rSurf */
      assert(x_of_X(x,X,patchp));
      x[0] -= patchp->c[0];
      x[1] -= patchp->c[1];
      x[2] -= patchp->c[2];
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
        interp_s->K = patchp->n[2]-1;
      }
      else 
       Error0(NO_OPTION);
       
      /* f value at r = r0 = rSurf*/
      interp_s->field = patchp->pool[Ind(extrap->fld[f]->f)];
      plan_interpolation(interp_s);
      fr0 = execute_interpolation(interp_s);
       
      /* df/dx value */
      for (d1 = 0; d1 < 3; d1++)
      {
        interp_s->field = patchp->pool[Ind(extrap->fld[f]->df[d1])];
        plan_interpolation(interp_s);
        df_dx[d1] = execute_interpolation(interp_s);
      }
      
      /* d^2f/dx^2 value */
      for (d1 = 0; d1 < 3; d1++)
      {
        for (d2 = d1; d2 < 3; d2++)
        {
          interp_s->field = 
            patchp->pool[Ind(extrap->fld[f]->ddf[IJsymm3(d1,d2)])];
          plan_interpolation(interp_s);
          ddf_ddx[IJsymm3(d1,d2)] = execute_interpolation(interp_s);
        }
      }
      
      /* free */
      free_interpolation(interp_s);
       
      /* df/dr */
      dfdr = (N[0]*df_dx[0]+N[1]*df_dx[1]+N[2]*df_dx[2]);
      
      ddfddr = 0;
      _ddfddr[0] = _ddfddr[1] = _ddfddr[2] = 0;
      /* d^2f/dr^2 */
      for (int _i = 0; _i < 3; ++_i)
      {
        for (int _j = 0; _j < 3; ++_j)
        {
          _ddfddr[_i] += (KD[_i==_j]/rSurf - x[_i]*x[_j]/rSurf3)*df_dx[_j];
          _ddfddr[_i] += N[_j]*ddf_ddx[IJsymm3(_i,_j)];
        }
        ddfddr += _ddfddr[_i]*N[_i]; \
      }
      
      /* extrap */
      demand->r     = r;
      demand->r0    = rSurf;
      demand->fr0   = fr0;
      demand->dfr0  = dfdr;
      demand->ddfr0 = ddfddr;
      
      field->v[ijk] = extrap->extrap(demand);
     }
     ++f;
    }
  }
  return EXIT_SUCCESS;
}

/* collect names of the fields and their derivatives */
static void collect_names(struct Extrap_S *const extrap,
                          const char **const fields_name,
                          const unsigned nf)
{
  const char *s = 0;
  unsigned f,i,j;
    
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

