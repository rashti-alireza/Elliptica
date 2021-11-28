#include "core_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "manifold_lib.h"
#include "error_handling_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "maths_linear_algebra_lib.h"
#include "fields_lib.h"
#include "maths_equation_solvings_lib.h"

#define _MAX_STR_ (400)

/* basic trig math */
#define J__Cos(a) cos((a))
#define J__Sin(a) sin((a))
#define J__Csc(a) (1./sin((a)))
#define J__Cot(a) (1./tan((a)))

/* 2*M_PI */
#define J__2M_PI (6.283185307179586)

/* Kronecker Delta */
#define J__KD(i,j)  ( (i)==(j) ? 1. : 0.)

/* eta_i in Elliptica's paper, NOTE: n = patch->n[?] - 1 */
#define J__eta(i,n) ( (i) == 0 || (i) == (n) ? 1. : 2. )

/* (-1)^n */
#define J__sign(n) (((n)%2) ? -1. : 1.)

/* quick theta: 
// NOTE1: assuming Chebyshev Extrema points.
// NOTE2: assuming patch is defined. */
#define J__theta(i,X_axis) ( (i)*M_PI/((patch->n[(X_axis)])-1.) )

/* normalization, NOTE: n = patch->n */
#define J__norm(n) ( 0.5/((n)-1.) )

/* dX/dx */
#define J__dX_dx(patch,ijk,dX_axis,dx_axis) \
  ( (patch)->JacobianT->dX_dx[(dX_axis)][(dx_axis)][(ijk)] )


/* d2X/dxdy */
#define J__d2X_dxdy(patch,ijk,dX_axis,dxdy_axis) \
  ( (patch)->JacobianT->d2X_dxdy[(dX_axis)][(dxdy_axis)][(ijk)] )

/* dN/dX */
#define J__dN_dX(patch,dX_axis) \
  ( (patch)->JacobianT->dN_dX[(dX_axis)] )


/* sum_{n=0}^{N} cos(n lambda) = 
// 0.5 + 0.5*( sin( (N+0.5)*(lambda) ) ) / ( sin( 0.5*(lambda) ) ),
// N0 = N+0.5. */
#define J__sum_0_N_cos_nlambda(N,N0,lambda) \
  ( \
    0.5 + 0.5*( sin( (N0)*(lambda) ) ) / ( sin( 0.5*(lambda) ) ) \
  )

/* d/dlambda sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define J__d_dlambda_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J__2M_PI) ?\
    (0.0):\
    ( \
      J__Csc(0.5*(lambda))*(2.*(N0)*J__Cos((lambda)*(N0)) - \
      J__Cot(0.5*(lambda))*J__Sin((lambda)*(N0)))\
    )*0.25\
  )

/* d^2/dlambda^2 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define J__d2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J__2M_PI) ?\
    -(Pow3(N)/3.+Pow2(N)/2.+N/6.)/* simplified, don't forget - sign! */ :\
    (\
      J__Csc(0.5*(lambda))*(-4.*(N0)*J__Cos((lambda)*(N0))*J__Cot(0.5*(lambda)) + \
      (-1. - 4.*Pow2(N0) + 2.*Pow2(J__Csc(0.5*(lambda))))*J__Sin((lambda)*(N0)))\
    )*0.125\
  )

/* d^3/dlambda^3 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define J__d3_dlambda3_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J__2M_PI) ?\
    (0.0):\
    (\
      pow(J__Csc((lambda)/2.,3))*(2*(N0)*\
        (9 - 4*Pow2((N0)) + (3 + 4*Pow2((N0)))*J__Cos((lambda)))*\
        J__Cos((lambda)*(N0)) - (11 - 12*Pow2((N0)) + J__Cos((lambda)) + \
          12*Pow2((N0))*J__Cos((lambda)))*J__Cot((lambda)/2.)*J__Sin((lambda)*(N0)))\
    )/32.\
  )

/* d^4/dlambda^4 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define J__d4_dlambda4_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J__2M_PI) ?\
    (Pow4(N)*(N/5.+0.5)+Pow3(N)/3.-N/30.)/* simplified */:\
    (\
      pow(J__Csc((lambda)/2.),5)*(-16*(N0)*\
        (11 - 4*Pow2((N0)) + J__Cos((lambda)) + \
          4*Pow2((N0))*J__Cos((lambda)))*J__Cos((lambda)*(N0))*J__Sin((lambda)) + \
       (115 - 120*Pow2((N0)) + 48*Pow4((N0)) + \
          (76 + 96*Pow2((N0)) - 64*Pow4((N0)))*J__Cos((lambda)) + \
          (1 + 24*Pow2((N0)) + 16*Pow4((N0)))*J__Cos(2*(lambda)))*\
        J__Sin((lambda)*(N0)))\
    )/256.\
  )


/* -> eta_j{d/dX(df/du)=d/dX (2*sum_0^N (Tn(Xj) Tn(X)) -1 -(-1)^j *T_{N}(X))},
// NOTE: X = cos(th), N = patch->n[?]-1. */
#define J__d2f_dudX(thi,thj,N,j) \
   (J__eta(j,N)*( d_dXi_2xsum_0_N_Tnj_Tni(thi,thj,N) - J__sign(j)*dT_dx((int)(N),cos(thi)) ))

/* -> d2/dX^2(df/du)=d2/dX^2 (2*sum_0^N (Tn(Xj) Tn(X)) -1 -(-1)^j *T_{N}(X)),
// NOTE: X = cos(th)), N = patch->n[?]-1. */
#define J__d3f_dudXdX(thi,thj,N,j) \
   (J__eta(j,N)*( d2_dXi2_2xsum_0_N_Tnj_Tni(thi,thj,N) - J__sign(j)*d2T_dx2((int)(N),cos(thi)) ))

/* normalization * coords jacobian * J__d2f_dudX */
#define J__d2f_dudx(patch, dx_axis, X_axis, ijk, lmn, qi,qj) \
  ( J__norm(patch->n[X_axis])*J__dX_dx(patch,ijk,X_axis,dx_axis)*J__dN_dX(patch,X_axis)*\
    J__d2f_dudX(J__theta(qi,X_axis),J__theta(qj,X_axis),patch->n[X_axis]-1,qj) )

/* normalization * coords jacobian * J__d3f_dudXdX */
#define J__d3f_dudxdy(patch, dx_axis, dy_axis, dxdy_axis,X_axis, ijk, lmn, qi,qj) \
  ( \
    J__norm(patch->n[X_axis])*J__dN_dX(patch,X_axis)*\
    ( \
      J__d2X_dxdy(patch,ijk,X_axis,dxdy_axis)*\
      J__d2f_dudX(J__theta(qi,X_axis),J__theta(qj,X_axis),patch->n[X_axis]-1,qj) +     \
      J__dX_dx(patch,ijk,X_axis,dx_axis)*J__dX_dx(patch,ijk,X_axis,dy_axis)*J__dN_dX(patch,X_axis)* \
      J__d3f_dudXdX(J__theta(qi,X_axis),J__theta(qj,X_axis),patch->n[X_axis]-1,qj)     \
    )\
  )


/* Jacobian type */
typedef enum JTYPE_E
{
  T_x/* df_x/df */,
  T_xx/* df_xx/df */,
  T_xy/* df_xy/df */,
  T_xz/* df_xz/df */,
  T_y/* df_y/df */,
  T_yy/* df_yy/df */,
  T_yz/* df_yz/df */,
  T_z/* df_z/df */,
  T_zz/* df_zz/df */,
  T_UNDEF
}JType_E;

/* Jacobain for equation */
typedef void Js_Jacobian_eq_F(double **const J,Patch_T *const patch,JType_E jt_e);

void free_patch_SolMan_jacobian(Patch_T *const patch);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
void make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void obsolete_fill_jacobian_spectral_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type);
double read_matrix_entry_ccs(Matrix_T *const m, const long r,const long c);
fdInterp_dfs_T *get_dInterp_df(const Patch_T *const patch,const SubFace_T *const sf,const char *const dir);
static JType_E str2JType_E(const char *const str);
static char *interpret_type(const char *const type);
static void JType_E2str(const JType_E e,char *const str);
static void make_jacobian_spectral_method(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_spectral_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_spectral_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir);
static void make_jacobian_direct_method(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_direct_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_direct_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir);
static double ChebExtrema_1point(const Uint n, const Uint p);
static double dc_df(const Uint n,const Uint i,const Uint l);
static double dT_dx(const int n,const double x);
static void read_1st_and_2nd_deriv(const JType_E deriv_dir,JType_E *const deriv_1st,JType_E *const deriv_2nd);
static void JType_E2Dd_T(const JType_E jt_e, Dd_T *const q_dir);
static void write_J_in_disk_ccs(void);
static double J_sizeMb_ccs(const Matrix_T *const m);
fJs_T *get_j_reader(const Matrix_T *const m);
static void coarse_grain_Ap_ccs_matrix(Matrix_T *const m,const int Nslice) __attribute__((unused)) ;
static double dInterp_x_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_y_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_z_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_x_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_y_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_z_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_x_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_y_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_z_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_x_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_y_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_z_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);

static double
d2_dXi2_2xsum_0_N_Tnj_Tni(double thi/* X_i = cos(theta_i) */,
                          double thj/* X_i = cos(theta_i) */,
                          Uint N/* the sum upper limit */);
static double
d_dXi_2xsum_0_N_Tnj_Tni(double thi/* X_i = cos(theta_i) */,
                        double thj/* X_i = cos(theta_i) */,
                        Uint N/* the sum upper limit */);

double
  d2f_dxdu_spectral_Jacobian_analytic(Patch_T *const patch,
                                      const Uint dx_axis, 
                                      const Uint ijk,const Uint lmn);
  
double
  d3f_dxdydu_spectral_Jacobian_analytic(Patch_T *const patch,
                                        const int dxdy_axis,
                                        const Uint ijk,const Uint lmn);


