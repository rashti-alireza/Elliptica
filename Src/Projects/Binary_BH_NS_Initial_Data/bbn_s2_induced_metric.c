
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
      /* find patch and X,Y,Z at NS surface in which theta and phi take place */
      find_XYZ_and_patch_of_theta_phi_NS_CS(X,&patch,theta,phi,grid);

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