/*
// Alireza Rashti
// July 2020
*/

#include "bbn_s2_induced_metric.h"

/* inclusion map S2->M and computing AKV vector from the derivatives.
// here we use (theta,phi) coords on sub-manifold S2 and 
// (x,y,z) = r(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) on M */
void
bbn_compute_AKV_from_z
  (
  Grid_T *const grid/* grid */,
  const double *const akv/* akv scalar values */,
  const char *const dakv_D0/* d/dx akv name */,
  const char *const dakv_D1/* d/dy akv name */,
  const char *const dakv_D2/* d/dz akv name */,
  const char *const type/* NS or BH */,
  const unsigned Ntheta/* number of points in theta direction */,
  const unsigned Nphi/* number of points in theta direction */,
  const unsigned lmax/* l max in Ylm, if asked for spherical harmonic */,
  const int expansion_type/* 1 double fourier, 0: spherical harmonic */
  )
{
  FUNC_TIC
  
  unsigned (*surface_patch)(const Patch_T *const patch) = 0;
  unsigned i,j,k,p;
  
  if (strcmp_i(type,"BH"))
    surface_patch = IsItHorizonPatch;
  else if (strcmp_i(type,"NS"))
    surface_patch = IsItNSSurface;
  else
    Error0("No such type.");
 
  if (expansion_type == 1)
  {
    double *realC,*imagC;
    
    /* find FT coeffs */
    r2cft_2d_coeffs_S2(akv,Ntheta,Nphi,&realC,&imagC,1);
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      if (!surface_patch(patch))
        continue;
      
      const unsigned *N = patch->n;
      Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
      const double *X;
      double theta = 0,phi = 0,iz;
      ADD_AND_ALLOC_FIELD(__z_scalar_TEMP)
      DECLARE_FIELD(__z_scalar_TEMP)
      
      /* populate z scalar in 3d */
      for (i = 0; i < N[0]; ++i)
      {
        for (j = 0; j < N[1]; ++j)
        {
          X = patch->node[L(N,i,j,0)]->X;
          find_theta_phi_of_XYZ_CS(&theta,&phi,X,side);
          
          iz = r2cft_2d_interpolation_S2(realC,imagC,Ntheta,Nphi,theta,phi);
          
          for (k = 0; k < N[2]; ++k)
            __z_scalar_TEMP->v[L(N,i,j,k)] = iz;
        }
      }/* for (i = 0; i < N[0]; ++i) */
      
      /* compute derivatives */
      Field_T *dAKV_D0 = patch->pool[Ind(dakv_D0)];
      Field_T *dAKV_D1 = patch->pool[Ind(dakv_D1)];
      Field_T *dAKV_D2 = patch->pool[Ind(dakv_D2)];
      empty_field(dAKV_D0);
      empty_field(dAKV_D1);
      empty_field(dAKV_D2);
      dAKV_D0->v = Partial_Derivative(__z_scalar_TEMP,"x");
      dAKV_D1->v = Partial_Derivative(__z_scalar_TEMP,"y");
      dAKV_D2->v = Partial_Derivative(__z_scalar_TEMP,"z");
      REMOVE_FIELD(__z_scalar_TEMP)
    }
    free(realC);
    free(imagC);
  }/* double Fourier */
  else if (expansion_type == 0)/* Ylm expansion */
  {
    double *realClm = alloc_ClmYlm(lmax);
    double *imagClm = alloc_ClmYlm(lmax);
      
    /* find Ylm coeffs */
    get_Ylm_coeffs(realClm,imagClm,akv,Ntheta,Nphi,lmax);
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      if (!surface_patch(patch))
        continue;
      
      const unsigned *N = patch->n;
      Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
      const double *X;
      double theta = 0,phi = 0,iz;
      ADD_AND_ALLOC_FIELD(__z_scalar_TEMP)
      DECLARE_FIELD(__z_scalar_TEMP)
      
      /* populate z scalar in 3d */
      for (i = 0; i < N[0]; ++i)
      {
        for (j = 0; j < N[1]; ++j)
        {
          X = patch->node[L(N,i,j,0)]->X;
          find_theta_phi_of_XYZ_CS(&theta,&phi,X,side);
          
          iz = interpolation_Ylm(realClm,imagClm,lmax,theta,phi);
          
          for (k = 0; k < N[2]; ++k)
            __z_scalar_TEMP->v[L(N,i,j,k)] = iz;
        }
      }/* for (i = 0; i < N[0]; ++i) */
      
      /* compute derivatives */
      Field_T *dAKV_D0 = patch->pool[Ind(dakv_D0)];
      Field_T *dAKV_D1 = patch->pool[Ind(dakv_D1)];
      Field_T *dAKV_D2 = patch->pool[Ind(dakv_D2)];
      empty_field(dAKV_D0);
      empty_field(dAKV_D1);
      empty_field(dAKV_D2);
      dAKV_D0->v = Partial_Derivative(__z_scalar_TEMP,"x");
      dAKV_D1->v = Partial_Derivative(__z_scalar_TEMP,"y");
      dAKV_D2->v = Partial_Derivative(__z_scalar_TEMP,"z");
      REMOVE_FIELD(__z_scalar_TEMP)
    }
    free(realClm);
    free(imagClm);
  }/* Ylm expansion */
  else
    Error0(NO_OPTION);
    
  FUNC_TOC
}

/* computing the induced metric on S2 (cubedspherical,CTS). */
void
bbn_compute_induced_metric_on_S2_CS_CTS
  (
  Grid_T *const grid/* grid */,
  const char *const type/* NS or BH */,
  const unsigned Ntheta/* number of points in theta direction */,
  const unsigned Nphi/* number of points in phi direction */,
  const unsigned lmax/* l max in Ylm */,
  double **const ph_D0D0/* induced h00  pointer */,
  double **const ph_D0D1/* induced h01  pointer */,
  double **const ph_D1D1/* induced h11  pointer */,
  const int expansion_type/* 1 double fourier, 0: spherical harmonic */
  )
{

  if (expansion_type == 1)
  {
   compute_induced_metric_on_S2_CS_FT_CTS
     (grid,type,Ntheta,Nphi,ph_D0D0,ph_D0D1,ph_D1D1);
  }
  else if (expansion_type == 0)
  {
   compute_induced_metric_on_S2_CS_Ylm_CTS
     (grid,type,lmax,ph_D0D0,ph_D0D1,ph_D1D1);
  }
  else
    Error0(NO_OPTION);

}

/* computing the induced metric on S2 (cubedspherical,Ylm,CTS).
// here we use (theta,phi) coords on sub-manifold S2 and 
// (x,y,z) = r(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) 
// coords for the manifold. 
// also assuming CTS method is used; thus metic gamma = psi^4*_gamma. 
// note: it allocates memory for the induced metric.
// note: it uses Ylm collocation points */
static void 
compute_induced_metric_on_S2_CS_Ylm_CTS
  (
  Grid_T *const grid,
  const char *const type,/* NS or BH */
  const unsigned lmax,/* l max in Ylm */
  double **const ph_D0D0,/* induced h00  pointer */
  double **const ph_D0D1,/* induced h01  pointer */
  double **const ph_D1D1 /* induced h11  pointer */
  )
{
  FUNC_TIC
  const unsigned N = Pow2(2*lmax+1);
  double *const h_D0D0 = alloc_double(N);
  double *const h_D0D1 = alloc_double(N);
  double *const h_D1D1 = alloc_double(N);
  const unsigned Ntheta= 2*lmax+1;
  const unsigned Nphi  = 2*lmax+1;
  double theta,phi;
  unsigned i,j;
  int type_flg = -1;/* 1 for NS , 0 for BH */
  
  if (strcmp_i(type,"BH"))
    type_flg = 0;
  else if (strcmp_i(type,"NS"))
    type_flg = 1;
  else
    Error0("Bad argument: no such type.");
    
  /* initialize tables */
  init_Legendre_root_function();
  
  /* for each points of Ylm find manifold metric g_D?D? */
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      
      unsigned ij = IJ(i,j);
      Patch_T *patch = 0;
      double ipsi,r,x[3],X[3],
            i_gamma_D0D0,i_gamma_D0D1,i_gamma_D0D2,
            i_gamma_D1D1,i_gamma_D1D2,i_gamma_D2D2,
            gamma_D0D0,gamma_D0D1,gamma_D0D2,
            gamma_D1D1,gamma_D1D2,gamma_D2D2;
      /* find patch and X,Y,Z on the surface in which theta and phi take place */
      find_XYZ_and_patch_of_theta_phi_CS(X,&patch,theta,phi,grid,type);

      /* finding x */
      x_of_X(x,X,patch);
      x[0] -= patch->c[0];
      x[1] -= patch->c[1];
      x[2] -= patch->c[2];
      r     = root_square(3,x,0);
      
      /* find value at the (X,Y,Z) */
      INTERPOLATE_macro(psi);
      INTERPOLATE_macro(_gamma_D0D0);
      INTERPOLATE_macro(_gamma_D0D1);
      INTERPOLATE_macro(_gamma_D0D2);
      INTERPOLATE_macro(_gamma_D1D1);
      INTERPOLATE_macro(_gamma_D1D2);
      INTERPOLATE_macro(_gamma_D2D2);
      
      /* gamma = psi^4 * _gamma */
      double ipsi4 = pow(ipsi,4);
      gamma_D0D0 = ipsi4*i_gamma_D0D0;
      gamma_D0D1 = ipsi4*i_gamma_D0D1;
      gamma_D0D2 = ipsi4*i_gamma_D0D2;
      gamma_D1D1 = ipsi4*i_gamma_D1D1;
      gamma_D1D2 = ipsi4*i_gamma_D1D2;
      gamma_D2D2 = ipsi4*i_gamma_D2D2;
      
      bbn_populate_2d_induced_metric_S2_theta_phi(
         &h_D0D0[ij],&h_D0D1[ij],&h_D1D1[ij],
         &gamma_D0D0,&gamma_D0D1,&gamma_D0D2,
         &gamma_D1D1,&gamma_D1D2,&gamma_D2D2,
         r,theta,phi);
      
    }/* end of for (j = 0; j < Nphi; ++j) */
  }/* end for (i = 0; i < Ntheta; ++i) */

  *ph_D0D0 = h_D0D0;
  *ph_D0D1 = h_D0D1;
  *ph_D1D1 = h_D1D1;
  
  FUNC_TOC
}

/* computing the induced metric on S2 (cubedspherical,FT,CTS).
// here we use (theta,phi) coords on sub-manifold S2 and 
// (x,y,z) = r(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) 
// coords for the manifold. 
// also assuming CTS method is used; thus metic gamma = psi^4*_gamma. 
// note: it allocates memory for the induced metric.
// note: this using double Fourier on S2 collocation points. */
static void 
compute_induced_metric_on_S2_CS_FT_CTS
  (
  Grid_T *const grid,
  const char *const type,/* NS or BH */
  const unsigned Ntheta,/* # of collocation points in theta direction */
  const unsigned Nphi,/* # of collocation points in phi direction */
  double **const ph_D0D0,/* induced h00  pointer */
  double **const ph_D0D1,/* induced h01  pointer */
  double **const ph_D1D1 /* induced h11  pointer */
  )
{
  FUNC_TIC
  const unsigned N = Ntheta*Nphi;
  double *const h_D0D0 = alloc_double(N);
  double *const h_D0D1 = alloc_double(N);
  double *const h_D1D1 = alloc_double(N);
  double theta,phi;
  unsigned i,j;
  int type_flg = -1;/* 1 for NS , 0 for BH */
  
  if (strcmp_i(type,"BH"))
    type_flg = 0;
  else if (strcmp_i(type,"NS"))
    type_flg = 1;
  else
    Error0("Bad argument: no such type.");
    
  
  /* for each points find manifold metric g_D?D? */
  for (i = 0; i < Ntheta; ++i)
  {
    theta = i*M_PI/(Ntheta-1);
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      
      unsigned ij = IJ(i,j);
      Patch_T *patch = 0;
      double ipsi,r,x[3],X[3],
            i_gamma_D0D0,i_gamma_D0D1,i_gamma_D0D2,
            i_gamma_D1D1,i_gamma_D1D2,i_gamma_D2D2,
            gamma_D0D0,gamma_D0D1,gamma_D0D2,
            gamma_D1D1,gamma_D1D2,gamma_D2D2;
      /* find patch and X,Y,Z on the surface in which theta and phi take place */
      find_XYZ_and_patch_of_theta_phi_CS(X,&patch,theta,phi,grid,type);

      /* finding x */
      x_of_X(x,X,patch);
      x[0] -= patch->c[0];
      x[1] -= patch->c[1];
      x[2] -= patch->c[2];
      r     = root_square(3,x,0);
      
      /* find value at the (X,Y,Z) */
      INTERPOLATE_macro(psi);
      INTERPOLATE_macro(_gamma_D0D0);
      INTERPOLATE_macro(_gamma_D0D1);
      INTERPOLATE_macro(_gamma_D0D2);
      INTERPOLATE_macro(_gamma_D1D1);
      INTERPOLATE_macro(_gamma_D1D2);
      INTERPOLATE_macro(_gamma_D2D2);
      
      /* gamma = psi^4 * _gamma */
      double ipsi4 = pow(ipsi,4);
      gamma_D0D0 = ipsi4*i_gamma_D0D0;
      gamma_D0D1 = ipsi4*i_gamma_D0D1;
      gamma_D0D2 = ipsi4*i_gamma_D0D2;
      gamma_D1D1 = ipsi4*i_gamma_D1D1;
      gamma_D1D2 = ipsi4*i_gamma_D1D2;
      gamma_D2D2 = ipsi4*i_gamma_D2D2;
      
      bbn_populate_2d_induced_metric_S2_theta_phi(
         &h_D0D0[ij],&h_D0D1[ij],&h_D1D1[ij],
         &gamma_D0D0,&gamma_D0D1,&gamma_D0D2,
         &gamma_D1D1,&gamma_D1D2,&gamma_D2D2,
         r,theta,phi);
      
    }/* end of for (j = 0; j < Nphi; ++j) */
  }/* end for (i = 0; i < Ntheta; ++i) */

  *ph_D0D0 = h_D0D0;
  *ph_D0D1 = h_D0D1;
  *ph_D1D1 = h_D1D1;
  
  FUNC_TOC  
}


/* given theta, phi on the surface and type (BH or NS)
// it finds the corresponding patch and X,Y,Z coordinate. */
static void 
find_XYZ_and_patch_of_theta_phi_CS
  (
  double *const X,/* fill X coords */
  Patch_T **const ppatch,/* fill patch */
  const double theta,/* (theta,phi) */
  const double phi,/* (theta,phi) */
  Grid_T *const grid,/* in this grid */
  const char *const type/* [BH,NS] */
  )
{
  const double tan_phi    = tan(phi);
  const double cos_theta  = cos(theta);
  const double tan_phi2   = Pow2(tan_phi);
  const double cos_theta2 = Pow2(cos_theta);
  Flag_T found_flg = NO;
  unsigned (*surface_patch)(const Patch_T *const patch) = 0;
  unsigned p;
  
  
  if (strcmp_i(type,"BH"))
  {
    surface_patch = IsItHorizonPatch;
    X[2] = 0;/* since we are on BH surface from BH surrounding side */
  }
  else if (strcmp_i(type,"NS"))
  {
    surface_patch = IsItNSSurface;
    X[2] = 1;/* since we are on NS surface */
  }
  else
    Error0("No such type.");
  
  /* check all of patches in which (x,y,z) and 
  // (X,Y,Z) and (theta,phi) are consistent */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!surface_patch(patch))
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
    else if (EQL(theta2,0) && side == UP)
    {
      found_flg = YES;
      *ppatch = patch;
      break;
    }
    else if (EQL(theta2,M_PI) && side == DOWN)
    {
      found_flg = YES;
      *ppatch = patch;
      break;
    }
  }
  if (found_flg == NO)
    Error0("(X,Y,Z) or patch could not be found.\n");
}

/* given (X,Y,Z) in the specified slice of NS/BH in cubed spherical coords
// it finds the associated polar and azimuthal angels on the surface  */
static void find_theta_phi_of_XYZ_CS(double *const theta,double *const phi,const double *const X,const Flag_T side)
{
  const double a = X[0];
  const double b = X[1];
  const double d = sqrt(1+Pow2(a)+Pow2(b));
  
  switch (side)
  {
    case UP:
      *phi   = arctan(b,a);
      *theta = acos(1/d);
    break;
    case DOWN:
      *phi   = arctan(a,b);
      *theta = acos(-1/d);
    break;
    case LEFT:
      *phi   = arctan(-1,a);
      *theta = acos(b/d);
    break;
    case RIGHT:
      *phi   = arctan(1,b);
      *theta = acos(a/d);
    break;
    case BACK:
      *phi   = arctan(b,-1);
      *theta = acos(a/d);
    break;
    case FRONT:
      *phi   = arctan(a,1);
      *theta = acos(b/d);
    break;
    default:
      Error0(NO_OPTION);
  }
  
}

/* test induced metric h algorithm.
// it tests both NS and BH for a concrete example:
// 3-d => ds^2 = dx^2+dy^2+dz^2+dx*dy+dx*dz+dy*dz.
// NOTE: it changes the values of _gamma and psi. */
void bbn_test_induced_metric_algorithm(Grid_T *const grid)
{
  FUNC_TIC
  
  unsigned lmax   = UINT_MAX;
  unsigned Ntheta = (unsigned)Pgeti("akv_n_theta");
  unsigned Nphi   = (unsigned)Pgeti("akv_n_phi");

  if (Pcmps("akv_expansion","spherical_harmonic"))
  {
    lmax = (unsigned)Pgeti("akv_lmax");/* lmax in Ylm */
  }

  const char *type     = 0;
  double *h_D0D0=0,*h_D0D1=0,*h_D1D1=0;/* induced metric */
  double x[3],X[3],r,theta,phi;
  unsigned p,i,j,ij;
  int status;/* 0 success, 1 failed */
  int type_flg;/* NS = 1, BH = 0 */
  
  /* change psi and _gamma to sphere */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn,ijk;
    nn = patch->nn;
    
    REALLOC_v_WRITE_v(_gamma_D0D2)
    REALLOC_v_WRITE_v(_gamma_D0D0)
    REALLOC_v_WRITE_v(_gamma_D0D1)
    REALLOC_v_WRITE_v(_gamma_D1D2)
    REALLOC_v_WRITE_v(_gamma_D1D1)
    REALLOC_v_WRITE_v(_gamma_D2D2)
    REALLOC_v_WRITE_v(psi)
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      _gamma_D0D0[ijk] = 1;
      _gamma_D0D1[ijk] = 0.5;
      _gamma_D0D2[ijk] = 0.5;
      _gamma_D1D2[ijk] = 0.5;
      _gamma_D1D1[ijk] = 1;
      _gamma_D2D2[ijk] = 1;
      psi[ijk]         = 0.7;/* arbitrary */
    }
  }/* end of FOR_ALL_PATCHES */
  
  /* test NS */
  status   = 0;
  type     = "NS";
  type_flg = 1;
  
  if (Pcmps("akv_expansion","spherical_harmonic"))
  {
    compute_induced_metric_on_S2_CS_Ylm_CTS
      (grid,type,lmax,&h_D0D0,&h_D0D1,&h_D1D1);
    for (i = 0; i < Ntheta; ++i)
    {
      theta = acos(-Legendre_root_function(i,Ntheta));
      for (j = 0; j < Nphi; ++j)
      {
        phi = j*2*M_PI/Nphi;
        ij  = IJ(i,j);
        Patch_T *patch = 0;
        find_XYZ_and_patch_of_theta_phi_CS(X,&patch,theta,phi,grid,type);
        /* finding x */
        x_of_X(x,X,patch);
        x[0] -= patch->c[0];
        x[1] -= patch->c[1];
        x[2] -= patch->c[2];
        r     = root_square(3,x,0);
        
        double ipsi;
        INTERPOLATE_macro(psi)
        double ipsi4 = pow(ipsi,4);
        
        double h00 = 
          Pow2(r)*ipsi4*(Power(Cos(phi),2)*Power(Cos(theta),2) + 
          Power(Cos(theta),2)*Power(Sin(phi),2) + 
          Cos(phi)*Cos(theta)*(Cos(theta)*Sin(phi) - Sin(theta)) - 
          Cos(theta)*Sin(phi)*Sin(theta) + Power(Sin(theta),2));
          
        double h01 = 
          0.5*Pow2(r)*ipsi4*(Cos(phi) - Sin(phi))*Sin(theta)*
          (Cos(phi)*Cos(theta) + Cos(theta)*Sin(phi) - Sin(theta));
        
        double h11 =
          Pow2(r)*ipsi4*(1 - Cos(phi)*Sin(phi))*Power(Sin(theta),2);
          
        if (!EQL(h_D0D0[ij],h00))
        {
          printf("dh00 = %g\n",h_D0D0[ij]-h00);
          status = 1;
        }
        if (!EQL(h_D0D1[ij],h01)) 
        {
          printf("dh01 = %g\n",h_D0D1[ij]-h01);
          status = 1;
        }
        if (!EQL(h_D1D1[ij],h11))
        {
          printf("dh11 = %g\n",h_D1D1[ij]-h11);
          status = 1;
        }
      }
    }
  }/* if (Pcmps("akv_expansion","spherical_harmonic")) */
  else if (Pcmps("akv_expansion","double_fourier"))
  {
    compute_induced_metric_on_S2_CS_FT_CTS
      (grid,type,Ntheta,Nphi,&h_D0D0,&h_D0D1,&h_D1D1);
    for (i = 0; i < Ntheta; ++i)
    {
      theta = i*M_PI/(Ntheta-1);
      for (j = 0; j < Nphi; ++j)
      {
        phi = j*2*M_PI/Nphi;
        ij  = IJ(i,j);
        Patch_T *patch = 0;
        find_XYZ_and_patch_of_theta_phi_CS(X,&patch,theta,phi,grid,type);
        /* finding x */
        x_of_X(x,X,patch);
        x[0] -= patch->c[0];
        x[1] -= patch->c[1];
        x[2] -= patch->c[2];
        r     = root_square(3,x,0);
        
        double ipsi;
        INTERPOLATE_macro(psi)
        double ipsi4 = pow(ipsi,4);
        
        double h00 = 
          Pow2(r)*ipsi4*(Power(Cos(phi),2)*Power(Cos(theta),2) + 
          Power(Cos(theta),2)*Power(Sin(phi),2) + 
          Cos(phi)*Cos(theta)*(Cos(theta)*Sin(phi) - Sin(theta)) - 
          Cos(theta)*Sin(phi)*Sin(theta) + Power(Sin(theta),2));
          
        double h01 = 
          0.5*Pow2(r)*ipsi4*(Cos(phi) - Sin(phi))*Sin(theta)*
          (Cos(phi)*Cos(theta) + Cos(theta)*Sin(phi) - Sin(theta));
        
        double h11 =
          Pow2(r)*ipsi4*(1 - Cos(phi)*Sin(phi))*Power(Sin(theta),2);
          
        if (!EQL(h_D0D0[ij],h00))
        {
          printf("dh00 = %g\n",h_D0D0[ij]-h00);
          status = 1;
        }
        if (!EQL(h_D0D1[ij],h01)) 
        {
          printf("dh01 = %g\n",h_D0D1[ij]-h01);
          status = 1;
        }
        if (!EQL(h_D1D1[ij],h11))
        {
          printf("dh11 = %g\n",h_D1D1[ij]-h11);
          status = 1;
        }
      }
    }
  }/* else if (Pcmps("akv_expansion","double_fourier")) */
  else
    Error0(NO_OPTION);
  
  
  printf("~> Testing induced metric algorithm %s side:",type);
  if (status)
    printf(" => [FAILED] <=\n");
  else
    printf(" => [PASSED] <=\n");
  
  free(h_D0D0); h_D0D0 = 0;
  free(h_D0D1); h_D0D1 = 0;
  free(h_D1D1); h_D1D1 = 0;
  
  /* test BH */
  status   = 0;
  type     = "BH";
  type_flg = 0;
  
  if (Pcmps("akv_expansion","spherical_harmonic"))
  {
    compute_induced_metric_on_S2_CS_Ylm_CTS
      (grid,type,lmax,&h_D0D0,&h_D0D1,&h_D1D1);
    for (i = 0; i < Ntheta; ++i)
    {
      theta = acos(-Legendre_root_function(i,Ntheta));
      for (j = 0; j < Nphi; ++j)
      {
        phi = j*2*M_PI/Nphi;
        ij  = IJ(i,j);
        Patch_T *patch = 0;
        find_XYZ_and_patch_of_theta_phi_CS(X,&patch,theta,phi,grid,type);
        /* finding x */
        x_of_X(x,X,patch);
        x[0] -= patch->c[0];
        x[1] -= patch->c[1];
        x[2] -= patch->c[2];
        r     = root_square(3,x,0);
        
        double ipsi;
        INTERPOLATE_macro(psi)
        double ipsi4 = pow(ipsi,4);
        
        double h00 = 
          Pow2(r)*ipsi4*(Power(Cos(phi),2)*Power(Cos(theta),2) + 
          Power(Cos(theta),2)*Power(Sin(phi),2) + 
          Cos(phi)*Cos(theta)*(Cos(theta)*Sin(phi) - Sin(theta)) - 
          Cos(theta)*Sin(phi)*Sin(theta) + Power(Sin(theta),2));
          
        double h01 = 
          0.5*Pow2(r)*ipsi4*(Cos(phi) - Sin(phi))*Sin(theta)*
          (Cos(phi)*Cos(theta) + Cos(theta)*Sin(phi) - Sin(theta));
        
        double h11 =
          Pow2(r)*ipsi4*(1 - Cos(phi)*Sin(phi))*Power(Sin(theta),2);
          
        if (!EQL(h_D0D0[ij],h00))
        {
          printf("dh00 = %g\n",h_D0D0[ij]-h00);
          status = 1;
        }
        if (!EQL(h_D0D1[ij],h01)) 
        {
          printf("dh01 = %g\n",h_D0D1[ij]-h01);
          status = 1;
        }
        if (!EQL(h_D1D1[ij],h11))
        {
          printf("dh11 = %g\n",h_D1D1[ij]-h11);
          status = 1;
        }
      }
    }
  }/* if (Pcmps("akv_expansion","spherical_harmonic")) */
  else if (Pcmps("akv_expansion","double_fourier"))
  {
    compute_induced_metric_on_S2_CS_FT_CTS
      (grid,type,Ntheta,Nphi,&h_D0D0,&h_D0D1,&h_D1D1);
    for (i = 0; i < Ntheta; ++i)
    {
      theta = i*M_PI/(Ntheta-1);
      for (j = 0; j < Nphi; ++j)
      {
        phi = j*2*M_PI/Nphi;
        ij  = IJ(i,j);
        Patch_T *patch = 0;
        find_XYZ_and_patch_of_theta_phi_CS(X,&patch,theta,phi,grid,type);
        /* finding x */
        x_of_X(x,X,patch);
        x[0] -= patch->c[0];
        x[1] -= patch->c[1];
        x[2] -= patch->c[2];
        r     = root_square(3,x,0);
        
        double ipsi;
        INTERPOLATE_macro(psi)
        double ipsi4 = pow(ipsi,4);
        
        double h00 = 
          Pow2(r)*ipsi4*(Power(Cos(phi),2)*Power(Cos(theta),2) + 
          Power(Cos(theta),2)*Power(Sin(phi),2) + 
          Cos(phi)*Cos(theta)*(Cos(theta)*Sin(phi) - Sin(theta)) - 
          Cos(theta)*Sin(phi)*Sin(theta) + Power(Sin(theta),2));
          
        double h01 = 
          0.5*Pow2(r)*ipsi4*(Cos(phi) - Sin(phi))*Sin(theta)*
          (Cos(phi)*Cos(theta) + Cos(theta)*Sin(phi) - Sin(theta));
        
        double h11 =
          Pow2(r)*ipsi4*(1 - Cos(phi)*Sin(phi))*Power(Sin(theta),2);
          
        if (!EQL(h_D0D0[ij],h00))
        {
          printf("dh00 = %g\n",h_D0D0[ij]-h00);
          status = 1;
        }
        if (!EQL(h_D0D1[ij],h01)) 
        {
          printf("dh01 = %g\n",h_D0D1[ij]-h01);
          status = 1;
        }
        if (!EQL(h_D1D1[ij],h11))
        {
          printf("dh11 = %g\n",h_D1D1[ij]-h11);
          status = 1;
        }
      }
    }
  }/* else if (Pcmps("akv_expansion","double_fourier")) */
  else
    Error0(NO_OPTION);
    
  printf("~> Testing induced metric algorithm %s side:",type);
  if (status)
    printf(" => [FAILED] <=\n");
  else
    printf(" => [PASSED] <=\n");
  
  free(h_D0D0); h_D0D0 = 0;
  free(h_D0D1); h_D0D1 = 0;
  free(h_D1D1); h_D1D1 = 0;
  FUNC_TOC
}


