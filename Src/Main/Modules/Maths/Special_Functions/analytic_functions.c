/*
// Alireza Rashti
// July 2018
*/

#include "analytic_functions.h"
#define CONST -30.2340
/* these used for compatibility with CForm produced by mathematica */
#define E M_E
#define xM x_(i) /* x Macro */
#define yM y_(i) /* y Macro */
#define zM z_(i) /* z Macro */
#define Cos(a) cos(a)
#define Sin(a) sin(a)
#define Cosh(a) cosh(a)
#define Sinh(a) sinh(a)
#define Log(a) log(a)
#define Power(a,b) pow(a,b)
#define Sqrt(a) sqrt(a)

/* the followings are some analytic functions used for various purposes.
// ->return value: a pointer to value of function on whole patch.
*/

/* f: patch -> const */
double *c_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
      f[i] = CONST;
      
  return f;
}

/* f: patch -> x */
double *x_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = x_(i);
    
  return f;
}

/* f: patch -> y */
double *y_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = y_(i);
  return f;
}

/* f: patch -> z */
double *z_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = z_(i);
  return f;
}

/* f: patch -> x^3-y^3+z^3-x^2*y*z+z^2*y^2*x */
double *poly5_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = Power(xM,3) - Power(yM,3) - Power(xM,2)*yM*zM + xM*Power(yM,2)*Power(zM,2) + 
   Power(zM,3);
  return f;
}

/* f: (x,y,z) -> x^3-y^3+z^3-x^2*y*z+z^2*y^2*x */
double poly5_f_point(const double x,const double y,const double z)
{
  double f;
  f = Power(x,3) - Power(y,3) - Power(x,2)*y*z + x*Power(y,2)*Power(z,2) + 
      Power(z,3);
        
  return f;
}

/* f: patch -> r = sqrt(x^2+y^2+z^2) */
double *r_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
      f[i] = sqrt(pow(x_(i),2)+pow(y_(i),2)+pow(z_(i),2));
  return f;
}

/* f: patch -> 1/sqrt(1+r^2) */
double *inv_rP1_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
      f[i] = 1/sqrt(pow(x_(i),2)+pow(y_(i),2)+pow(z_(i),2)+1);
  return f;
}

/* f: patch -> cos(x) */
double *cosx_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
      f[i] = cos(x_(i));
  return f;
}

/* f: patch -> cos(y) */
double *cosy_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = cos(y_(i));
  return f;
}
/* f: patch -> cos(z) */
double *cosz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = cos(z_(i));
  return f;
}

/* f: patch -> sin(x) */
double *sinx_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = sin(x_(i));
  return f;
}

/* f: patch -> sin(y) */
double *siny_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = sin(y_(i));
  return f;
}

/* f: patch -> sin(z) */
double *sinz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = sin(z_(i));
  return f;
}

/* f: patch -> cos(xyz) */
double *cosxyz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = cos(x_(i)*y_(i)*z_(i));
  return f;
}

/* f: patch -> cos^4(xyz) */
double *cos4xyz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = pow(cos(x_(i)*y_(i)*z_(i)),4);
  return f;
}

/* f: patch -> cos^5(xyz) */
double *cos5xyz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = pow(cos(x_(i)*y_(i)*z_(i)),5);
  return f;
}

/* f: patch -> sin(xyz) */
double *sinxyz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = sin(x_(i)*y_(i)*z_(i));
  return f;
}

/* f: patch -> sin^3(xyz) */
double *sin3xyz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = pow(sin(x_(i)*y_(i)*z_(i)),3);
  return f;
}

/* f: patch -> cosh(xyz) */
double *coshxyz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = cosh(x_(i)*y_(i)*z_(i));
  return f;
}


/* f: patch -> sinh(xyz) */
double *sinhxyz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = sinh(x_(i)*y_(i)*z_(i));
  return f;
}


/* f: patch -> tanh(xyz) */
double *tanhxyz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = tanh(x_(i)*y_(i)*z_(i));
  return f;
}

/* f: patch -> log(1+(xyz)^2) */
double *logxyz_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = log(1+pow(x_(i)*y_(i)*z_(i),2));
  return f;
}


/* f: patch -> x^2*cos(x*y^3*z) +y^3*z^4*x*sin(x^2*z) */
double *mix1_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = pow(x_(i),2)*cos(x_(i)*pow(y_(i),3)*z_(i))+
         pow(y_(i),3)*pow(z_(i),4)*x_(i)*sin(pow(x_(i),2)*z_(i));
  return f;
}

/* f: patch -> log(1+z^2)*cosh{sqrt(x^2+y^2+z^2)+sin(exp(x*y*z))} */
double *mix2_f(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = log(1+pow(z_(i),2))*cosh(
      sqrt(pow(x_(i),2)+pow(y_(i),2)+pow(z_(i),2))+
        sin(exp(x_(i)*y_(i)*z_(i)))
                      );
    
  return f;
}

/* the following is the derivative of functions above, shown by _* */

/* f: c -> dc/dx */
double *c_f_x(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: c -> dc/dy */
double *c_f_y(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}


/* f: c -> dc/dz */
double *c_f_z(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}


/* f: x -> dx/dx */
double *x_f_x(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 1;
  return f;
}

/* f: x -> dx/dy */
double *x_f_y(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: x -> dx/dz */
double *x_f_z(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: x -> d2x/dxx */
double *x_f_xx(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: y -> dy/dx */
double *y_f_x(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: y -> dy/dy */
double *y_f_y(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 1;
  return f;
}

/* f: y -> dy/dz */
double *y_f_z(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: y -> d2y/dyy */
double *y_f_yy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: z -> dz/dx */
double *z_f_x(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: z -> dz/dy */
double *z_f_y(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: z -> dz/dz */
double *z_f_z(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 1;
  return f;
}

/* f: z -> d2z/dzz */
double *z_f_zz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: patch -> d(1/sqrt(1+r))/dx */
double *inv_rP1_f_x(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -(xM/Power(1 + Power(xM,2) + Power(yM,2) + Power(zM,2),1.5));
  return f;
}

/* f: patch -> d(1/sqrt(1+r))/dy */
double *inv_rP1_f_y(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -(yM/Power(1 + Power(xM,2) + Power(yM,2) + Power(zM,2),1.5));
  return f;
}

/* f: patch -> d(1/sqrt(1+r))/dz */
double *inv_rP1_f_z(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -(zM/Power(1 + Power(xM,2) + Power(yM,2) + Power(zM,2),1.5));
  return f;
}

/* f: patch -> d(1/sqrt(1+r))/dxx */
double *inv_rP1_f_xx(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (-1 + 2*Power(xM,2) - Power(yM,2) - Power(zM,2))/
   Power(1 + Power(xM,2) + Power(yM,2) + Power(zM,2),2.5);
  return f;
}

/* f: patch -> d(1/sqrt(1+r))/dyy */
double *inv_rP1_f_yy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (-1 - Power(xM,2) + 2*Power(yM,2) - Power(zM,2))/
   Power(1 + Power(xM,2) + Power(yM,2) + Power(zM,2),2.5);
  return f;
}

/* f: patch -> d(1/sqrt(1+r))/dzz */
double *inv_rP1_f_zz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (-1 - Power(xM,2) - Power(yM,2) + 2*Power(zM,2))/
   Power(1 + Power(xM,2) + Power(yM,2) + Power(zM,2),2.5);
  return f;
}

/* f: patch -> d(1/sqrt(1+r))/dxy */
double *inv_rP1_f_xy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (3*xM*yM)/Power(1 + Power(xM,2) + Power(yM,2) + Power(zM,2),2.5);
  return f;
}

/* f: patch -> d(1/sqrt(1+r))/dxz */
double *inv_rP1_f_xz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (3*xM*zM)/Power(1 + Power(xM,2) + Power(yM,2) + Power(zM,2),2.5);
  return f;
}

/* f: patch -> d(1/sqrt(1+r))/dyz */
double *inv_rP1_f_yz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (3*yM*zM)/Power(1 + Power(xM,2) + Power(yM,2) + Power(zM,2),2.5);
  return f;
}

/* f: patch -> d(1/sqrt(1+r))/dxyz */
double *inv_rP1_f_xyz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (-15*xM*yM*zM)/Power(1 + Power(xM,2) + Power(yM,2) + Power(zM,2),3.5);
  return f;
}

/* f: r -> dr/dx */
double *r_f_x(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = xM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2));
  return f;
}

/* f: r -> dr/dy */
double *r_f_y(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = yM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2));
  return f;
}

/* f: r -> dr/dz */
double *r_f_z(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = zM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2));
  return f;
}

/* f: r -> dr/dxx */
double *r_f_xx(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (Power(yM,2) + Power(zM,2))/
   Power(Power(xM,2) + Power(yM,2) + Power(zM,2),1.5);
  return f;
}

/* f: r -> dr/dyy */
double *r_f_yy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (Power(xM,2) + Power(zM,2))/
   Power(Power(xM,2) + Power(yM,2) + Power(zM,2),1.5);
  return f;
}

/* f: r -> dr/dzz */
double *r_f_zz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (Power(xM,2) + Power(yM,2))/
   Power(Power(xM,2) + Power(yM,2) + Power(zM,2),1.5);
  return f;
}

/* f: r -> dr/dxy */
double *r_f_xy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -((xM*yM)/Power(Power(xM,2) + Power(yM,2) + Power(zM,2),1.5));
  return f;
}

/* f: r -> dr/dxz */
double *r_f_xz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -((xM*zM)/Power(Power(xM,2) + Power(yM,2) + Power(zM,2),1.5));
  return f;
}

/* f: r -> dr/dyz */
double *r_f_yz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -((yM*zM)/Power(Power(xM,2) + Power(yM,2) + Power(zM,2),1.5));
  return f;
}

/* f: r -> dr/dxyz */
double *r_f_xyz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (3*xM*yM*zM)/Power(Power(xM,2) + Power(yM,2) + Power(zM,2),2.5);
  return f;
}

/* f: sin(xyz)->d(sin(xyz))/dx */
double *sinxyz_f_x(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = yM*zM*Cos(xM*yM*zM);
  return f;
}

/* f: sin(xyz)->d(sin(xyz))/dy */
double *sinxyz_f_y(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = xM*zM*Cos(xM*yM*zM);
  return f;
}

/* f: sin(xyz)->d(sin(xyz))/dz */
double *sinxyz_f_z(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = xM*yM*Cos(xM*yM*zM);
  return f;
}

/* f: sin(xyz)->d(sin(xyz))/dxx */
double *sinxyz_f_xx(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -(Power(yM,2)*Power(zM,2)*Sin(xM*yM*zM));
  return f;
}

/* f: sin(xyz)->d(sin(xyz))/dyy */
double *sinxyz_f_yy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -(Power(xM,2)*Power(zM,2)*Sin(xM*yM*zM));
  return f;
}

/* f: sin(xyz)->d(sin(xyz))/dzz */
double *sinxyz_f_zz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -(Power(xM,2)*Power(yM,2)*Sin(xM*yM*zM));
  return f;
}

/* f: sin(xyz)->d(sin(xyz))/dxy */
double *sinxyz_f_xy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = zM*(Cos(xM*yM*zM) - xM*yM*zM*Sin(xM*yM*zM));
  return f;
}

/* f: sin(xyz)->d(sin(xyz))/dxz */
double *sinxyz_f_xz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = yM*(Cos(xM*yM*zM) - xM*yM*zM*Sin(xM*yM*zM));
  return f;
}

/* f: sin(xyz)->d(sin(xyz))/dyz */
double *sinxyz_f_yz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = xM*(Cos(xM*yM*zM) - xM*yM*zM*Sin(xM*yM*zM));
  return f;
}

/* f: sin(xyz)->d(sin(xyz))/dxyz*/
double *sinxyz_f_xyz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (1 - Power(xM,2)*Power(yM,2)*Power(zM,2))*Cos(xM*yM*zM) - 
   3*xM*yM*zM*Sin(xM*yM*zM);
  return f;
}

/* f: patch -> d(poly5)_dx */
double *poly5_f_x(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 3*Power(xM,2) - 2*xM*yM*zM + Power(yM,2)*Power(zM,2);
  return f;
}

/* f: patch -> d(poly5)_dy */
double *poly5_f_y(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -3*Power(yM,2) - Power(xM,2)*zM + 2*xM*yM*Power(zM,2);
  return f;
}

/* f: patch -> d(poly5)_dz */
double *poly5_f_z(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -(Power(xM,2)*yM) + 2*xM*Power(yM,2)*zM + 3*Power(zM,2);
  return f;
}

/* f: patch -> d(poly5)_dxx */
double *poly5_f_xx(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 6*xM - 2*yM*zM;
  return f;
}

/* f: patch -> d(poly5)_dyy */
double *poly5_f_yy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -6*yM + 2*xM*Power(zM,2);
  return f;
}

/* f: patch -> d(poly5)_dzz */
double *poly5_f_zz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 2*xM*Power(yM,2) + 6*zM;
  return f;
}

/* f: patch -> d(poly5)_dxy */
double *poly5_f_xy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 2*zM*(-xM + yM*zM);
  return f;
}

/* f: patch -> d(poly5)_dxz */
double *poly5_f_xz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 2*yM*(-xM + yM*zM);
  return f;
}

/* f: patch -> d(poly5)_dyz */
double *poly5_f_yz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -(xM*(xM - 4*yM*zM));
  return f;
}

/* f: patch -> d(poly5)_dxyz */
double *poly5_f_xyz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -2*xM + 4*yM*zM;
  return f;
}

/* f: dmix2/dx */
double *mix2_f_x(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (xM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
     Power(E,xM*yM*zM)*yM*zM*Cos(Power(E,xM*yM*zM)))*
   Log(1 + Power(zM,2))*Sinh(Sqrt(Power(xM,2) + Power(yM,2) + 
     Power(zM,2)) + Sin(Power(E,xM*yM*zM)));
  return f;
}

/* f: dmix2/dy */
double *mix2_f_y(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (yM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
     Power(E,xM*yM*zM)*xM*zM*Cos(Power(E,xM*yM*zM)))*
   Log(1 + Power(zM,2))*Sinh(Sqrt(Power(xM,2) + Power(yM,2) + 
     Power(zM,2)) + Sin(Power(E,xM*yM*zM)));
  return f;
}


/* f: dmix2/dz */
double *mix2_f_z(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (2*zM*Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))))/(1 + Power(zM,2)) + 
   (zM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*yM*Cos(Power(E,xM*yM*zM)))*
    Log(1 + Power(zM,2))*Sinh(Sqrt(Power(xM,2) + Power(yM,2) + 
      Power(zM,2)) + Sin(Power(E,xM*yM*zM)));
  return f;
}

/* f: dmix2/dxx */
double *mix2_f_xx(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = Log(1 + Power(zM,2))*(Power(xM/
       Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Power(E,xM*yM*zM)*yM*zM*Cos(Power(E,xM*yM*zM)),2)*
    Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))) + 
     ((Power(yM,2) + Power(zM,2))/
       Power(Power(xM,2) + Power(yM,2) + Power(zM,2),1.5) + 
      Power(E,xM*yM*zM)*Power(yM,2)*Power(zM,2)*
       Cos(Power(E,xM*yM*zM)) - 
      Power(E,2*xM*yM*zM)*Power(yM,2)*Power(zM,2)*
       Sin(Power(E,xM*yM*zM)))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))));
    
  return f;
}


/* f: dmix2/dyy */
double *mix2_f_yy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = Log(1 + Power(zM,2))*(Power(yM/
       Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Power(E,xM*yM*zM)*xM*zM*Cos(Power(E,xM*yM*zM)),2)*
    Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))) + 
     ((Power(xM,2) + Power(zM,2))/
       Power(Power(xM,2) + Power(yM,2) + Power(zM,2),1.5) + 
      Power(E,xM*yM*zM)*Power(xM,2)*Power(zM,2)*
       Cos(Power(E,xM*yM*zM)) - 
      Power(E,2*xM*yM*zM)*Power(xM,2)*Power(zM,2)*
       Sin(Power(E,xM*yM*zM)))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))));
  return f;
}

/* f: dmix2/dzz */
double *mix2_f_zz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (-4*Power(zM,2)*Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))))/Power(1 + Power(zM,2),2) + 
   (2*Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))))/(1 + Power(zM,2)) + 
   Power(zM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*yM*Cos(Power(E,xM*yM*zM)),2)*
    Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)))*Log(1 + Power(zM,2)) + 
   (4*zM*(zM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Power(E,xM*yM*zM)*xM*yM*Cos(Power(E,xM*yM*zM)))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))))/(1 + Power(zM,2)) + 
   Log(1 + Power(zM,2))*((Power(xM,2) + Power(yM,2))/
     Power(Power(xM,2) + Power(yM,2) + Power(zM,2),1.5) + 
    Power(E,xM*yM*zM)*Power(xM,2)*Power(yM,2)*
     Cos(Power(E,xM*yM*zM)) - 
    Power(E,2*xM*yM*zM)*Power(xM,2)*Power(yM,2)*
     Sin(Power(E,xM*yM*zM)))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)));
    
  return f;
}


/* f: dmix2/dxy */
double *mix2_f_xy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (zM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*yM*Cos(Power(E,xM*yM*zM)))*
    (yM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*zM*Cos(Power(E,xM*yM*zM)))*
    Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)))*Log(1 + Power(zM,2)) + 
   (2*zM*(yM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Power(E,xM*yM*zM)*xM*zM*Cos(Power(E,xM*yM*zM)))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))))/(1 + Power(zM,2)) + 
   Log(1 + Power(zM,2))*(Power(E,xM*yM*zM)*xM*(1 + xM*yM*zM)*
     Cos(Power(E,xM*yM*zM)) + 
    yM*zM*(-Power(Power(xM,2) + Power(yM,2) + Power(zM,2),-1.5) - 
       Power(E,2*xM*yM*zM)*Power(xM,2)*Sin(Power(E,xM*yM*zM))))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)));
    
  return f;
}

/* f: dmix2/dxz */
double *mix2_f_xz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (zM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*yM*Cos(Power(E,xM*yM*zM)))*
    (xM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*yM*zM*Cos(Power(E,xM*yM*zM)))*
    Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)))*Log(1 + Power(zM,2)) + 
   (2*zM*(xM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Power(E,xM*yM*zM)*yM*zM*Cos(Power(E,xM*yM*zM)))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))))/(1 + Power(zM,2)) + 
   Log(1 + Power(zM,2))*(-((xM*zM)/
       Power(Power(xM,2) + Power(yM,2) + Power(zM,2),1.5)) + 
    Power(E,xM*yM*zM)*yM*Cos(Power(E,xM*yM*zM)) + 
    Power(E,xM*yM*zM)*xM*Power(yM,2)*zM*Cos(Power(E,xM*yM*zM)) - 
    Power(E,2*xM*yM*zM)*xM*Power(yM,2)*zM*Sin(Power(E,xM*yM*zM)))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)));
    
  return f;
}

/* f: dmix2/dyz */
double *mix2_f_yz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (zM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*yM*Cos(Power(E,xM*yM*zM)))*
    (yM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*zM*Cos(Power(E,xM*yM*zM)))*
    Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)))*Log(1 + Power(zM,2)) + 
   (2*zM*(yM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Power(E,xM*yM*zM)*xM*zM*Cos(Power(E,xM*yM*zM)))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))))/(1 + Power(zM,2)) + 
   Log(1 + Power(zM,2))*(Power(E,xM*yM*zM)*xM*(1 + xM*yM*zM)*
     Cos(Power(E,xM*yM*zM)) + 
    yM*zM*(-Power(Power(xM,2) + Power(yM,2) + Power(zM,2),-1.5) - 
       Power(E,2*xM*yM*zM)*Power(xM,2)*Sin(Power(E,xM*yM*zM))))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)));
  return f;
}

/* f: dmix2/dxyz */
double *mix2_f_xyz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = (2*zM*(yM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Power(E,xM*yM*zM)*xM*zM*Cos(Power(E,xM*yM*zM)))*
    (xM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Power(E,xM*yM*zM)*yM*zM*Cos(Power(E,xM*yM*zM)))*
    Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))))/(1 + Power(zM,2)) + 
   (xM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*yM*zM*Cos(Power(E,xM*yM*zM)))*
    Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)))*Log(1 + Power(zM,2))*
    (Power(E,xM*yM*zM)*xM*(1 + xM*yM*zM)*Cos(Power(E,xM*yM*zM)) + 
    yM*zM*(-Power(Power(xM,2) + Power(yM,2) + Power(zM,2),-1.5) - 
       Power(E,2*xM*yM*zM)*Power(xM,2)*Sin(Power(E,xM*yM*zM)))) + 
   (yM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*zM*Cos(Power(E,xM*yM*zM)))*
    Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)))*Log(1 + Power(zM,2))*
    (Power(E,xM*yM*zM)*yM*(1 + xM*yM*zM)*Cos(Power(E,xM*yM*zM)) + 
    xM*zM*(-Power(Power(xM,2) + Power(yM,2) + Power(zM,2),-1.5) - 
       Power(E,2*xM*yM*zM)*Power(yM,2)*Sin(Power(E,xM*yM*zM)))) + 
   (zM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*yM*Cos(Power(E,xM*yM*zM)))*
    Cosh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)))*Log(1 + Power(zM,2))*
    (Power(E,xM*yM*zM)*zM*(1 + xM*yM*zM)*Cos(Power(E,xM*yM*zM)) + 
    xM*yM*(-Power(Power(xM,2) + Power(yM,2) + Power(zM,2),-1.5) - 
       Power(E,2*xM*yM*zM)*Power(zM,2)*Sin(Power(E,xM*yM*zM)))) + 
   (zM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*yM*Cos(Power(E,xM*yM*zM)))*
    (yM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*xM*zM*Cos(Power(E,xM*yM*zM)))*
    (xM/Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Power(E,xM*yM*zM)*yM*zM*Cos(Power(E,xM*yM*zM)))*
    Log(1 + Power(zM,2))*Sinh(Sqrt(Power(xM,2) + Power(yM,2) + 
      Power(zM,2)) + Sin(Power(E,xM*yM*zM))) + 
   (2*zM*(Power(E,xM*yM*zM)*zM*(1 + xM*yM*zM)*
       Cos(Power(E,xM*yM*zM)) + 
      xM*yM*(-Power(Power(xM,2) + Power(yM,2) + Power(zM,2),-1.5) - 
         Power(E,2*xM*yM*zM)*Power(zM,2)*Sin(Power(E,xM*yM*zM))))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
      Sin(Power(E,xM*yM*zM))))/(1 + Power(zM,2)) + 
   Log(1 + Power(zM,2))*(-(Power(E,xM*yM*zM)*
       (-1 - 3*xM*yM*zM + (-1 + Power(E,2*xM*yM*zM))*Power(xM,2)*
        Power(yM,2)*Power(zM,2))*Cos(Power(E,xM*yM*zM))) + 
    3*xM*yM*zM*(Power(Power(xM,2) + Power(yM,2) + Power(zM,2),
        -2.5) - Power(E,2*xM*yM*zM)*(1 + xM*yM*zM)*
        Sin(Power(E,xM*yM*zM))))*
    Sinh(Sqrt(Power(xM,2) + Power(yM,2) + Power(zM,2)) + 
    Sin(Power(E,xM*yM*zM)));
        
  return f;
}

/* f: patch -> dsin(x)/dx */
double *sinx_f_x(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = cos(xM);
  return f;
}

/* f: patch -> dsin(x)/dy */
double *sinx_f_y(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: patch -> dsin(x)/dz */
double *sinx_f_z(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: patch -> dsin(x)/dxx */
double *sinx_f_xx(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = -sin(xM);
  return f;
}

/* f: patch -> dsin(x)/dyy */
double *sinx_f_yy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: patch -> dsin(x)/dzz */
double *sinx_f_zz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: patch -> dsin(x)/dxy */
double *sinx_f_xy(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: patch -> dsin(x)/dxz */
double *sinx_f_xz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: patch -> dsin(x)/dyz */
double *sinx_f_yz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

/* f: patch -> dsin(x)/dxyz */
double *sinx_f_xyz(Patch_T *const patch)
{
  double *f = alloc_double(patch->nn);
  Uint i;
  
  FOR_ALL_POINTS(i,patch)
    f[i] = 0;
  return f;
}

