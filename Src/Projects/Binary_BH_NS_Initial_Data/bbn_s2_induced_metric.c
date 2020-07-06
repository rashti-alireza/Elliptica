
/* computing the induced metric on S2 (cubedspherical,Ylm,CTS).
// here we use (theta,phi) coords on sub-manifold S2 and 
// (x,y,z) = r(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) 
// coords for the manifold. 
// also assuming CTS method is used thus metic gamma = psi*_gamma. */
void 
bbn_compute_induced_metric_on_S2_CS_Ylm_CTS
  (
  Grid_T *const grid,
  const char *const type,/* NS or BH */
  const unsigned lmax,/* l max in Ylm */
  double *const h_D0D0,/* induced h00  */
  double *const h_D0D1,/* induced h01  */
  double *const h_D1D1 /* induced h11  */
  )
{
  double theta,phi;
  unsigned Ntheta,Nphi;
  
  /* initialize tables */
  init_Legendre_root_function();
  Ntheta  = Nphi = 2*lmax+1;
  
  /* for each points of Ylm find manifold metric g_D?D? */
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      
      Patch_T *patch = 0;
      /* find patch and X,Y,Z on the surface in which theta and phi take place */
      find_XYZ_and_patch_of_theta_phi_CS(X,&patch,theta,phi,grid,type);

      /* find enthalpy at the (X,Y,Z) */
      Interpolation_T *interp_psi = init_interpolation();
      interp_h->field = patch->pool[Ind("psi")];
      interp_h->XY_dir_flag  = 1;
      interp_h->X            = X[0];
      interp_h->Y            = X[1];
      interp_h->K            = patch->n[2]-1;
      plan_interpolation(interp_h);
      h = execute_interpolation(interp_h);/* enthalpy */

      assert(h > 0);

      free_interpolation(interp_h);

      
      
    }/* end of for (j = 0; j < Nphi; ++j) */
  }/* end for (i = 0; i < Ntheta; ++i) */

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
  
  X[2] = 0;/* since we are on BH surface from BH surrounding side */
  
  if (strcmp_i(type,"BH"))
    surface_patch = IsItHorizonPatch;
  else if (strcmp_i(type,"NS"))
    surface_patch = IsItNSSurface;
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
  }
  if (found_flg == NO)
    Error0("(X,Y,Z) or patch could not be found.\n");
}
