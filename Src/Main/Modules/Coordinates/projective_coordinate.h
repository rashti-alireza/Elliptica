#include "coordinate_shared_lib.h"
#include "maths_calculus_lib.h"
#include "maths_approximation_lib.h"

#define Power(a,b) pow(a,b)
#define Sqrt(a) sqrt(a)
#define Power2(a) SQR(a)
#define Power3(a) (a)*SQR(a)
#define Sqrt2 1.4142135623730951
#define dZ_du(u,w) dY_du(u,w)
#define dZ_dw(u,w) dY_dv(u,w)
#define dX_dw(u,w) dX_dv(u,w)


enum enum_dA_da
{
  da_dx = 0,
  da_dy,
  da_dz,
  db_dx,
  db_dy,
  db_dz,
  dc_dx,
  dc_dy,
  dc_dz,
  dA_da_UNDEFINED
};

void make_nodes_ProjectiveHemisphereUp_coord(Patch_T *const patch);
void make_nodes_ProjectiveHemisphereDown_coord(Patch_T *const patch);
void make_nodes_StereographicSphereLeft_coord(Patch_T *const patch);
void make_nodes_StereographicSphereRight_coord(Patch_T *const patch);
void make_JacobianT_ProjectiveHemisphere_coord(Patch_T *const patch);
double JT_ProjectiveHemisphere(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double dN0_dx_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X);
double dN0_dy_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X);
double dN0_dz_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X);
double dN1_dx_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X);
double dN1_dy_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X);
double dN1_dz_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X);
double dN2_dx_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X);
double dN2_dy_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X);
double dN2_dz_ProjectiveHemisphere_patch(Patch_T *const patch,const double *const X);
double dN0_dx_StereographicSphereRight_patch(Patch_T *const patch,const double *const X);
double dN0_dy_StereographicSphereRight_patch(Patch_T *const patch,const double *const X);
double dN0_dz_StereographicSphereRight_patch(Patch_T *const patch,const double *const X);
double dN1_dx_StereographicSphereRight_patch(Patch_T *const patch,const double *const X);
double dN1_dy_StereographicSphereRight_patch(Patch_T *const patch,const double *const X);
double dN1_dz_StereographicSphereRight_patch(Patch_T *const patch,const double *const X);
double dN2_dx_StereographicSphereRight_patch(Patch_T *const patch,const double *const X);
double dN2_dy_StereographicSphereRight_patch(Patch_T *const patch,const double *const X);
double dN2_dz_StereographicSphereRight_patch(Patch_T *const patch,const double *const X);
double dN0_dx_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X);
double dN0_dy_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X);
double dN0_dz_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X);
double dN1_dx_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X);
double dN1_dy_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X);
double dN1_dz_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X);
double dN2_dx_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X);
double dN2_dy_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X);
double dN2_dz_StereographicSphereLeft_patch(Patch_T *const patch,const double *const X);
static double dNi_dxj_ProjectiveHemisphere(Patch_T *const patch, const Dd_T Ni, const Dd_T xj,const double *const X);
static double dNi_dxj_StereographicSphereLeft(Patch_T *const patch, const Dd_T Ni, const Dd_T xj,const double *const X);
static double dNi_dxj_StereographicSphereRight(Patch_T *const patch, const Dd_T Ni, const Dd_T xj,const double *const X);
static double dN_dX(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e);
enum enum_dA_da get_dA_da(const Dd_T q2_e, const Dd_T q1_e);
static void R1_R2_derivative(Patch_T *const patch);
int x_of_X(double *const x,const double *const X,const Patch_T *const patch);
double interpolation_2d_PH(Field_T *const R, const Patch_T *const patch,const double *const X);
void make_JacobianT_StereographicSphereLeft_coord(Patch_T *const patch);
void make_JacobianT_StereographicSphereRight_coord(Patch_T *const patch);
double JT_StereographicSphere_Left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_StereographicSphere_Right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
static double dX_du(const double u, const double v);
static double dX_dv(const double u, const double v);
static double dY_du(const double u, const double v);
static double dY_dv(const double u, const double v);
double dq2_dq1(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);



