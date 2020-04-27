/*
// Alireza Rashti
// June 2019
*/

#include "integration.h"

/* synopsis:
// =========
//
// ** filling the integration struct **
// Integration_T *I = init_integration();
//
// *** for example if you want to use Spectral method: ***
// VOLUME integration:
// *** for integrating {fdV} over a whole patch in physical domain: ***
// Note: dV is made automatically using the given metric, however,
// the integrand f(x) must be provided fully by the user.
// I->type = "Integral{f(x)dV},Spectral"; # dV is sqrt(det(g)) dxdydz
// I->Spectral->f = field;# this is the Field f(x)
//
// # filling the metric (assumed symmetry):
// I->g00 = gamma_D0D0;
// I->g01 = gamma_D0D1;
// I->g02 = gamma_D0D2;
// I->g11 = gamma_D1D1;
// I->g12 = gamma_D1D2;
// I->g22 = gamma_D2D2;

// SURFACE integration:
// *** for integrating {fdS} over a whole physical hypersurface: ***
// Note: dS is made automatically using the given metric, however,
// the integrand f(x) must be provided fully by the user.
// for example: if we have intergral f(i)*dS(-i) one needs to
// do the contraction of f^i and normal vector on the surface
// BUT the element dS is made automatically out of the given metric.
//
// I->type = "Integral{f(x)dS},Spectral"; # dS is sqrt(det(g_{ij} dx^i/dY dx^j/dZ))dYdZ for hypersurface X = const.
//                                        # dS is sqrt(det(g_{ij} dx^i/dX dx^j/dZ))dXdZ for hypersurface Y = const.
//                                        # dS is sqrt(det(g_{ij} dx^i/dX dx^j/dY))dXdY for hypersurface Z = const.
// I->Spectral->f = field;# this is the Field f(x)
// # selecting the hypersurface for example:
// I->Spectral->X_surface = 1; # for hypersurface X = const.
// I->Spectral->I         = i; # the index shows  X = const. plane
// 
// I->Spectral->Y_surface = 1; # for hypersurface Y = const.
// I->Spectral->J         = j; # the index shows  Y = const. plane
//
// I->Spectral->Z_surface = 1; # for hypersurface Z = const.
// I->Spectral->K         = k; # the index shows  Z = const. plane
//
// # filling the metric (assumed symmetry):
// I->g00 = gamma_D0D0;
// I->g01 = gamma_D0D1;
// I->g02 = gamma_D0D2;
// I->g11 = gamma_D1D1;
// I->g12 = gamma_D1D2;
// I->g22 = gamma_D2D2;
//
// *** for example if you want to use Composite Simpson's Rule 1D: ***
// I->type = "Composite Simpson's Rule 1D"
// I->Composite_Simpson_1D->a = 1.5;// e.g
// I->Composite_Simpson_1D->b = M_PI;// e.g
// I->Composite_Simpson_1D->n = 11;// e.g
// I->Composite_Simpson_1D->f = array;
//
// *** for example if you want to use Gaussian Quadrature Chebyshev Extrema: ***
// I->type = "Gaussian Quadrature Chebyshev Extrema"
// I->GQ_ChebyshevExtrema->f = array;
// I->GQ_ChebyshevExtrema->n = 10;// the dimension of array
//
// *** for example if you want to use Gaussian Quadrature Lobatto method: ***
// I->type = "Gaussian Quadrature Lobatto"
// I->GQ_Lobatto->f = array;
// I->GQ_Lobatto->n = 10;// the dimension of array
//
// *** for example if you want to use Gaussian Quadrature Legendre method: ***
// I->type = "Gaussian Quadrature Legendre"
// I->GQ_Legendre->f = array;
// I->GQ_Legendre->n = 10;// the dimension of array
//
// having populated the I structure, now we have: 
// ** planning the appropriate function for integration **
// plan_integration(I);
//
// ** evaluating integration **
// double integral = execute_integration(I);
//
// ** freeing **
// free_integration(I);
*/


/* initializing an Integration_T struct with calloc.
// ->return value: a pristine struct */
Integration_T *init_integration(void)
{
  Integration_T *I = calloc(1,sizeof(*I));
  pointerEr(I);
  
  return I;
}

/* ->return value: integral */
double execute_integration(Integration_T *const I)
{
  return I->integration_func(I);
}

/* given the information it decides how to perform the integral */
void plan_integration(Integration_T *const I)
{
  Patch_T *patch;
  Coord_T coordsys;
  
  unsigned i;
  
  if (strcmp_i(I->type,"Composite Simpson's Rule 1D"))
  {
    I->integration_func = Composite_Simpson_1D;
  }
  else if (strcmp_i(I->type,"Gaussian Quadrature Chebyshev Extrema"))
  {
    I->integration_func = GaussQuadrature_ChebyshevExtrema;
  }
  else if (strcmp_i(I->type,"Gaussian Quadrature Lobatto"))
  {
    I->integration_func = GaussQuadrature_Lobatto;
  }
  else if (strcmp_i(I->type,"Gaussian Quadrature Legendre"))
  {
    I->integration_func = GaussQuadrature_Legendre;
  }
  else if (strcmp_i(I->type,"Integral{f(x)dV},Spectral"))/* means over the whole physical volume covered by the patch */
  {
    patch = I->Spectral->f->patch;
    coordsys = patch->coordsys;
    
    if (coordsys == Cartesian || coordsys == CubedSpherical)
    {
      for (i = 0; i < 3; ++i)
      {
        if (patch->basis[i]       == Chebyshev_Tn_BASIS &&
            patch->collocation[i] == Chebyshev_Extrema    )
        {
          I->integration_func = f_xyz_dV_Cheb_Ext_Spec;
        }
        else
          Error0(INCOMPLETE_FUNC);
      }
    }
    else
      Error0(NO_OPTION);

  }
  else if (strcmp_i(I->type,"Integral{f(x)dS},Spectral"))/* means over the whole specified slice in physical area in the patch */
  {
    patch = I->Spectral->f->patch;
    coordsys = patch->coordsys;
    unsigned count = 0;
    
    if (I->Spectral->X_surface) count++;
    if (I->Spectral->Y_surface) count++;
    if (I->Spectral->Z_surface) count++;
    if (count > 1)
      Error0("At 'Integral{f(x)dS},Spectral' more than one hypersurface was flagged.\n");
    
    if (coordsys == Cartesian || coordsys == CubedSpherical)
    {
      for (i = 0; i < 3; ++i)
      {
        if (patch->basis[i]       == Chebyshev_Tn_BASIS &&
            patch->collocation[i] == Chebyshev_Extrema    )
        {
          I->integration_func = f_xyz_dS_Cheb_Ext_Spec;
        }
        else
          Error0(INCOMPLETE_FUNC);
      }
    }
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
}

/* free the integral struct */
void free_integration(Integration_T *I)
{
  if (!I)
    return;
    
  free(I);
}

/* Composite Simpson's Rule in 1D
// ->return value: the result of inegral */
static double Composite_Simpson_1D(Integration_T *const I)
{
  if (I->Composite_Simpson_1D->n % 2 != 1)
    Error0("Composite Simpson's Rule requires odd number of points.\n");
    
  const double h = (I->Composite_Simpson_1D->b-I->Composite_Simpson_1D->a)/(I->Composite_Simpson_1D->n-1);
  const double *const f = I->Composite_Simpson_1D->f;
  double i0 = 0,i1 = 0,i2 = 0;
  unsigned j;
  
  I->err = fabs(I->Composite_Simpson_1D->b-I->Composite_Simpson_1D->a)/180*pow(h,4)*100;
  
  i0 = f[0]+f[I->Composite_Simpson_1D->n-1];
  for (j = 1; j < I->Composite_Simpson_1D->n-1; ++j)
  {
    if (j%2 == 0)
      i2 += f[j];
    else
      i1 += f[j];
  }
  
  return h*(i0+2*i2+4*i1)/3;
}

/* performing the integral \integral_{-1}^{1} dx f(x)/sqrt(1-x^2) 
// using Guassian Quadrature method.
// note: to get a good result for this method, you need many points
//       unless f(x) is a polynomial which works best.
// note: the collocation points for f(x) are Chebyshev Extrema
// note: the integration is from -1 to 1 then the expected order of f(x)
//       is f(-1) = f[0] and f(1) = f[n-1].
// ->return value: \integral_{-1}^{1} dx f(x)/sqrt(1-x^2) */
static double GaussQuadrature_ChebyshevExtrema(Integration_T *const I)
{
  double i0 = 0;
  const double *const f = I->GQ_ChebyshevExtrema->f;
  const unsigned n      = I->GQ_ChebyshevExtrema->n;
  const double   w      = M_PI/(n-1);
  double err = M_PI;
  unsigned i;
  
  err /= Factorial(2*(int)n);
  err *= L_inf(n,f);/* approximately */
  err /= pow(2,2*n-1);
  I-> err = err;/* note: this error is valid only for polynomial */
  
  for (i = 1; i <= n-2; ++i)
    i0 += f[i];
  i0 += (f[0]+f[n-1])/2;
  i0 *= w;
  
  return i0;
}

/* performing the integral \integral_{-1}^{1} dx f(x)
// using Guassian Quadrature Lobatto's method.
// note: the collocation points for f(x) are zeros of d(Legendre(x,n-1))/dx
// note: the integration is from -1 to 1 then the expected order of f(x)
//       is f(-1) = f[0] and f(1) = f[n-1].
// ->return value: \integral_{-1}^{1} f(x)dx */
static double GaussQuadrature_Lobatto(Integration_T *const I)
{
  double i0 = 0;
  const double *const f = I->GQ_Lobatto->f;
  const unsigned n      = I->GQ_Lobatto->n;
  const int ni          = (int)n;
  double (*w)(const double x, const unsigned n) = Lobatto_weight_function;
  double err;
  unsigned i;
  
  /* initializing root tables */
  init_Lobatto_root_function();
  /* initialzing associated_legendre needed for Lobatto_weight_function */
  init_associated_legendre();
  
  /* trying to tame err */
  err = 1./Factorial(2*ni-2);
  err *= n*pow(n-1,3);
  err *= pow(2,2*n-1)/Factorial(2*ni-2);
  err /= Factorial(2*ni-2);
  err *= Factorial(ni-2);
  err *= Factorial(ni-2);
  err /= (2*n-1);
  err *= Factorial(ni-2);
  err *= Factorial(ni-2);
  err *= L_inf(n,f);/* approximately */
  if (fabs(err) > 10.)/* if n gets big, we have overflow, so put it 0, which is not valid */
    err = 0.;
  I-> err = err;/* note: this error is valid only for polynomial */

  for (i = 1; i <= n-2; ++i)
    i0 += w(Lobatto_root_function(i-1,n-1),n)*f[i];
  i0 += 2*(f[0]+f[n-1])/(n*(n-1));
  
  return i0;
}

/* performing the integral \integral_{-1}^{1} dx f(x)
// using Guass Legendre quadrature method.
// note: the collocation points for f(x) are zeros of Legendre(n,x) 
//       for array f of dimension n.
// note: the integration is from -1 to 1 then the expected order of f(x)
//       is f(-1) = f[0] and f(1) = f[n-1].
// ->return value: \integral_{-1}^{1} f(x)dx */
static double GaussQuadrature_Legendre(Integration_T *const I)
{
  double i0 = 0;
  const double *const f = I->GQ_Legendre->f;
  const unsigned n      = I->GQ_Legendre->n;
  const int ni          = (int)n;
  double (*w)(const double x, const unsigned n) = Legendre_weight_function;
  double err;
  unsigned i;
  
  /* initializing root tables */
  init_Legendre_root_function();
  
  /* initialzing dPn/dx needed for Legendre_weight_function */
  init_dLegendre_dx();
  
  /* trying to tame err */
  err = 1./Factorial(2*ni);
  err *= Factorial(ni);
  err = pow(err,3);
  err /= (2*n+1);
  err *= Factorial(ni);
  err *= pow(2,2*n+1);
  err *= L_inf(n,f);/* approximately */
  
  if (fabs(err) > 10.)/* if n gets big, we have overflow, so put it 0, which is not valid */
    err = 0.;
  I-> err = err;/* note: this error is valid only for polynomial */

  for (i = 0; i < n; ++i)
    i0 += w(Legendre_root_function(i,n),n)*f[i];
  
  return i0;
}

/* ->return value: integrating over the whole physical volume of the patch
// which has collocation points of Chebyshev extrema and
// bases of Chebyshev, USING CARTESIAN COORDINATES.
// Thread Safe */
static double f_xyz_dV_Cheb_Ext_Spec(Integration_T *const I)
{
  const Field_T *const f = I->Spectral->f;
  Patch_T patch = make_temp_patch(f->patch);
  Field_T *F    = add_field("integrand","(3dim)",&patch,YES);
  const double *const fv = f->v;
  double *const Fv = F->v;
  const double *Fc;
  const unsigned nn       = patch.nn;
  const unsigned *const n = patch.n;
  double sum = 0.,g;
  unsigned ijk,i,j,k;
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double g00 = I->g00[ijk];
    double g01 = I->g01[ijk];
    double g02 = I->g02[ijk];
    double g11 = I->g11[ijk];
    double g12 = I->g12[ijk];
    double g22 = I->g22[ijk];
    g = g00*g11*g22 - g00*pow(g12, 2) - pow(g01, 2)*g22 + 
        2*g01*g02*g12 - pow(g02, 2)*g11;/* determinant */
    Fv[ijk] = fv[ijk]*sqrt(g)*J_xyzN0N1N2(&patch,ijk);
  }
  
  Fc = make_coeffs_3d(F);
  
  /* note: if i,j or k is odd Int_ChebTn is 0, so it only loops over even numbers */
  for (i = 0; i < n[0]; i += 2)
    for (j = 0; j < n[1]; j += 2)
      for (k = 0; k < n[2]; k += 2)
        sum += Fc[L(n,i,j,k)]*
               Int_ChebTn_OPTM(i,n[0])*
               Int_ChebTn_OPTM(j,n[1])*
               Int_ChebTn_OPTM(k,n[2]);
        
  sum *= -1.;/* - sign is for integration */
  
  remove_field(F);
  free_temp_patch(&patch);
  
  return sum;
}

/* ->return value: integrating over the whole specified physical 
// hypersurface of the patch that has collocation points of Chebyshev extrema and
// bases of Chebyshev, USING CARTESIAN COORDINATES.
// Thread Safe */
static double f_xyz_dS_Cheb_Ext_Spec(Integration_T *const I)
{
  const Field_T *const f = I->Spectral->f;
  Patch_T patch = make_temp_patch(f->patch);
  Field_T *F    = add_field("integrand","(3dim)",&patch,YES);
  const double *const fv = f->v;
  double *const Fv = F->v;
  const double *Fc;
  const unsigned *const n = patch.n;
  double sum = 0.;
  unsigned ijk,i,j,k;
  
  /* populate integrand at X hypresurface */
  if (I->Spectral->X_surface)
  {
    i = I->Spectral->I;
    for (j = 0; j < n[1]; ++j)
      for (k = 0; k < n[2]; ++k)
      {
        ijk     = L(n,i,j,k);
        Fv[ijk] = fv[ijk]*sqrt(det_h_xyzN1N2_Cheb_Ext(&patch,I,ijk));
      }
    
    Fc = make_coeffs_2d(F,1,2);
    /* note: if i,j or k is odd Int_ChebTn is 0, so it only loops over even numbers */
    for (j = 0; j < n[1]; j += 2)
      for (k = 0; k < n[2]; k += 2)
        sum += Fc[L(n,i,j,k)]*
               Int_ChebTn_OPTM(j,n[1])*
               Int_ChebTn_OPTM(k,n[2]);
    
  }
  /* populate integrand at Y hypresurface */
  else if (I->Spectral->Y_surface)
  {
    j = I->Spectral->J;
    for (i = 0; i < n[0]; ++i)
      for (k = 0; k < n[2]; ++k)
      {
        ijk     = L(n,i,j,k);
        Fv[ijk] = fv[ijk]*sqrt(det_h_xyzN0N2_Cheb_Ext(&patch,I,ijk));
      }
    
    Fc = make_coeffs_2d(F,0,2);
    /* note: if i,j or k is odd Int_ChebTn is 0, so it only loops over even numbers */
    for (i = 0; i < n[0]; i += 2)
      for (k = 0; k < n[2]; k += 2)
        sum += Fc[L(n,i,j,k)]*
               Int_ChebTn_OPTM(i,n[0])*
               Int_ChebTn_OPTM(k,n[2]);
  }
  /* populate integrand at Z hypresurface */
  else if(I->Spectral->Z_surface)
  {
    k = I->Spectral->K;
    for (i = 0; i < n[0]; ++i)
      for (j = 0; j < n[1]; ++j)
      {
        ijk     = L(n,i,j,k);
        Fv[ijk] = fv[ijk]*sqrt(det_h_xyzN0N1_Cheb_Ext(&patch,I,ijk));
      }
    
    Fc = make_coeffs_2d(F,0,1);
    /* note: if i,j or k is odd Int_ChebTn is 0, so it only loops over even numbers */
    for (i = 0; i < n[0]; i += 2)
      for (j = 0; j < n[1]; j += 2)
        sum += Fc[L(n,i,j,k)]*
               Int_ChebTn_OPTM(i,n[0])*
               Int_ChebTn_OPTM(j,n[1]);
  }
  else
    Error0(NO_OPTION);
  
  remove_field(F);
  free_temp_patch(&patch);
  
  return sum;
}

/* taking the integral of Chebyshev T(n,x) appears in Chebyshev expansion
// from {xi,xf}, namaly: f = c_{0}+c_{N-1}*T(N-1,x) + 2*sum c_{n}*T(n,x)
// note: n is the coeffs number in c_{n} and N is total number of coeffs.
// -> return value: integral_{xi}^{xf} T(n,x)dx */
double Integrate_ChebTn(const unsigned n,const double xi,const double xf)
{
  if (fabs(xi) > 1. || fabs(xf) > 1.)
    Error0("Bad argument for Int_ChebTn function.\n");
    
  const double i = acos(xi);
  const double f = acos(xf);
  double Int;
  
  if (n == 1)
    Int = (cos(2*f) - cos(2*i))/4.;
  else
    Int = (-cos(f)*cos(f*n) + cos(i)*cos(i*n) -
           n*sin(f)*sin(f*n) + n*sin(i)*sin(i*n))/
           (Pow2(n)-1.);
     
  return Int;
}

/* OPTIMIZED version of the integral of Chebyshev T(n,x) appears in Chebyshev expansion
// from {-1,1}, namaly: f = c_{0}+c_{N-1}*T(N-1,x) + sum c_{n}*T(n,x).
// in this version, we assume n is even, so it won't check it every time.
// note: n is the coeffs number in c_{n} and N is total number of coeffs.
// -> return value: integral_{-1}^{1} T(n,x)dx */
static double Int_ChebTn_OPTM(const unsigned n,const unsigned N)
{
  if (n == 0)
    return 2.;
  else if (n == N-1)
    return -2./(Pow2(n)-1.);
  else
    return -4./(Pow2(n)-1.);
  
  return DBL_MAX; 
}

/* -> return value: det (d(x,y,z)/d(N0,N1,N2)) , 
// note: det (d(x,y,z)/d(N0,N1,N2)) = 1/(det (d(N0,N1,N2)/d(x,y,z)))
// we use, this method, calculating d(N0,N1,N2)/d(x,y,z) is faster
// than d(x,y,z)/d(N0,N1,N2). */
static double J_xyzN0N1N2(Patch_T *const patch,const unsigned ijk)
{
  const double a00 = dq2_dq1(patch,_N0_,_x_,ijk);
  const double a01 = dq2_dq1(patch,_N0_,_y_,ijk);
  const double a02 = dq2_dq1(patch,_N0_,_z_,ijk);
  const double a10 = dq2_dq1(patch,_N1_,_x_,ijk);
  const double a11 = dq2_dq1(patch,_N1_,_y_,ijk);
  const double a12 = dq2_dq1(patch,_N1_,_z_,ijk);
  const double a20 = dq2_dq1(patch,_N2_,_x_,ijk);
  const double a21 = dq2_dq1(patch,_N2_,_y_,ijk);
  const double a22 = dq2_dq1(patch,_N2_,_z_,ijk);
  
  return 1./(a00 *a11 *a22  - a00 *a12 *a21  -
             a01 *a10 *a22  + a01 *a12 *a20  +
             a02 *a10 *a21  - a02 *a11 *a20  );
  
}
